# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""Boltz-2 complex prediction + ipSAE scoring for 322 RBX1 designs.

Mirrors an internal complex re-scoring stack (the Boltz-2
half) but pinned to the RBX1 target and the rbx1_gem_paper repo layout:

* **Target**: RBX1 (UniProt P62877, 108 aa) — chain ``A``, MSA computed via
  ``--use_msa_server``.
* **Binder**: each row in ``data/designs.csv`` — chain ``B``, ``msa: empty``.

Outputs live on a Modal Volume ``rbx1-rerun-results``:

* ``boltz2/{pb_id}.json`` — IPSAE + native Boltz confidence
* ``structures/boltz2/{pb_id}.cif`` — predicted complex
* ``raw_data/boltz2/{pb_id}.npz`` — raw PAE matrix + scalars

Usage (Modal CLI auth via ``modal token set …`` or env vars)::

    # Fire-and-forget on Modal; survives local disconnect with --detach
    cd <repo_root>
    modal run --detach scripts/modal/modal_boltz2_rbx1.py

    # Download finished results to disk (idempotent, skip already-on-disk)
    modal run scripts/modal/modal_boltz2_rbx1.py --download

    # Tune GPU / concurrency
    GPU=A100 CONCURRENCY=10 modal run --detach scripts/modal/modal_boltz2_rbx1.py
"""

from __future__ import annotations

import contextlib
import json
import os
from pathlib import Path

import modal

GPU = os.environ.get("GPU", "A100")
TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 60))
CONCURRENCY = int(os.environ.get("CONCURRENCY", 10))

# Modal app / volume names — overridable so external users can run with
# their own workspace naming without editing source.
APP_NAME = os.environ.get("MODAL_APP_NAME", "rbx1-boltz2")
RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "rbx1-rerun-results")
BOLTZ_VOLUME_NAME = os.environ.get("MODAL_BOLTZ_VOLUME", "boltz-models")

# RBX1 — UniProt P62877, 108 aa. Matches data/target/rbx1.fasta exactly.
RBX1_TARGET_SEQ = (
    "MAAAMDVDTPSGTNSGAGKKRFEVKKWNAVALWAWDIVVDNCAICRNHIMDLCIECQANQ"
    "ASATSEECTVAWGVCNHAFHFHCISRWLKTRQVCPLDNREWEFQKYGH"
)

PAE_CUTOFF = 15.0
DIST_CUTOFF = 8.0  # DunbrackLab pDockQ default

# ---------------------------------------------------------------------------
# IPSAE / pDockQ / LIS scoring — inlined from rescoring stack so the Modal
# image is self-contained. DO NOT EDIT to diverge from the upstream copy.
# ---------------------------------------------------------------------------


def _calc_d0(n_res: float, min_value: float = 1.0) -> float:
    """d0 normalization from Yang & Skolnick, PROTEINS 57:702-710 (2004)."""
    n_res = max(27.0, float(n_res))
    return max(min_value, 1.24 * (n_res - 15.0) ** (1.0 / 3.0) - 1.8)


def _ptm_func(pae_values, d0: float):
    return 1.0 / (1.0 + (pae_values / d0) ** 2.0)


def compute_ipsae(
    pae_matrix,
    structure_path: str,
    target_len: int,
    binder_len: int,
    target_chain: str = "A",
    binder_chain: str = "B",
    pae_cutoff: float = PAE_CUTOFF,
    dist_cutoff: float = DIST_CUTOFF,
) -> dict[str, float]:
    """Compute ipSAE, pDockQ, pDockQ2, LIS from PAE + CIF/PDB."""
    import numpy as np

    try:
        import gemmi
    except ImportError:
        return _empty_metrics()

    pae_matrix = np.asarray(pae_matrix, dtype=np.float64)
    total_len = target_len + binder_len

    try:
        st = gemmi.read_structure(str(structure_path))
    except Exception as e:
        print(f"Failed to read structure {structure_path}: {e}")
        return _empty_metrics()

    model = st[0]
    chain_names = [c.name for c in model]
    if target_chain not in chain_names or binder_chain not in chain_names:
        print(f"Chains {target_chain}/{binder_chain} not in {chain_names}")
        return _empty_metrics()

    def _extract(chain_id: str) -> list[list[float]]:
        out: list[list[float]] = []
        for res in model[chain_id]:
            atom = res.find_atom("CB", "*") if res.name != "GLY" else res.find_atom("CA", "*")
            if not atom:
                atom = res.find_atom("CA", "*")
            if atom:
                out.append([atom.pos.x, atom.pos.y, atom.pos.z])
        return out

    target_coords = _extract(target_chain)
    binder_coords = _extract(binder_chain)
    struct_target_len = len(target_coords)
    struct_binder_len = len(binder_coords)
    struct_total = struct_target_len + struct_binder_len

    if pae_matrix.shape[0] != struct_total:
        if pae_matrix.shape == (total_len, total_len) and total_len != struct_total:
            idx = list(range(struct_target_len)) + list(
                range(target_len, target_len + struct_binder_len)
            )
            if max(idx) < pae_matrix.shape[0]:
                pae_matrix = pae_matrix[np.ix_(idx, idx)]
            else:
                print(
                    f"Cannot reconcile PAE {pae_matrix.shape} with structure "
                    f"({struct_target_len}+{struct_binder_len})"
                )
                return _empty_metrics()
        else:
            print(f"PAE shape {pae_matrix.shape} != struct ({struct_total}) or seq ({total_len})")
            return _empty_metrics()

    target_len = struct_target_len
    binder_len = struct_binder_len
    total_len = struct_total

    cb_coords = np.array(target_coords + binder_coords)
    distances = np.sqrt(((cb_coords[:, None, :] - cb_coords[None, :, :]) ** 2).sum(axis=2))
    target_mask = np.arange(total_len) < target_len
    binder_mask = ~target_mask

    pae_bt = pae_matrix[binder_mask][:, target_mask]
    pae_tb = pae_matrix[target_mask][:, binder_mask]

    n0chn = total_len
    d0chn = _calc_d0(n0chn)
    valid_bt = pae_bt < pae_cutoff
    valid_tb = pae_tb < pae_cutoff

    ipsae_d0chn_bt = float(_ptm_func(pae_bt[valid_bt], d0chn).mean()) if valid_bt.any() else 0.0
    ipsae_d0chn_tb = float(_ptm_func(pae_tb[valid_tb], d0chn).mean()) if valid_tb.any() else 0.0

    n0dom_bt = int(valid_bt.sum())
    ipsae_d0dom_bt = (
        float(_ptm_func(pae_bt[valid_bt], _calc_d0(n0dom_bt)).mean()) if n0dom_bt > 0 else 0.0
    )
    n0dom_tb = int(valid_tb.sum())
    ipsae_d0dom_tb = (
        float(_ptm_func(pae_tb[valid_tb], _calc_d0(n0dom_tb)).mean()) if n0dom_tb > 0 else 0.0
    )

    def _d0res_scores(pae_block, n_rows: int) -> float:
        vals: list[float] = []
        for i in range(n_rows):
            row = pae_block[i]
            good = row[row < pae_cutoff]
            if good.size > 0:
                vals.append(float(_ptm_func(good, _calc_d0(good.size)).mean()))
        return max(vals) if vals else 0.0

    ipsae_d0res_bt = _d0res_scores(pae_bt, binder_len)
    ipsae_d0res_tb = _d0res_scores(pae_tb, target_len)

    iptm_d0chn_bt = float(_ptm_func(pae_bt, d0chn).mean())
    iptm_d0chn_tb = float(_ptm_func(pae_tb, d0chn).mean())
    iptm_af_bt = float(_ptm_func(pae_bt, 10.0).mean())
    iptm_af_tb = float(_ptm_func(pae_tb, 10.0).mean())

    dist_bt = distances[binder_mask][:, target_mask]
    interface_mask = dist_bt <= dist_cutoff
    n_interface = int(interface_mask.sum())

    pdockq = pdockq2 = 0.0
    if n_interface > 0:
        binder_iface = np.where(interface_mask.any(axis=1))[0] + target_len
        target_iface = np.where(interface_mask.any(axis=0))[0]
        iface_idx = np.concatenate([binder_iface, target_iface])

        plddt_list: list[float] = []
        for chain_id in [target_chain, binder_chain]:
            for res in model[chain_id]:
                ca = res.find_atom("CA", "*")
                if ca:
                    plddt_list.append(ca.b_iso)
        plddt = np.array(plddt_list)
        if len(plddt) == total_len:
            mean_iface_plddt = float(plddt[iface_idx].mean())
            x = mean_iface_plddt * np.log10(n_interface)
            pdockq = 0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018
            ptm_iface = float(_ptm_func(pae_bt[interface_mask], 10.0).mean())
            x2 = mean_iface_plddt * ptm_iface
            pdockq2 = 1.31 / (1 + np.exp(-0.075 * (x2 - 84.733))) + 0.005

    pae_lis = pae_bt[pae_bt < 12.0]
    lis = float(((12.0 - pae_lis) / 12.0).mean()) if pae_lis.size > 0 else 0.0

    return {
        "ipsae_d0res_min": min(ipsae_d0res_bt, ipsae_d0res_tb),
        "ipsae_d0res_max": max(ipsae_d0res_bt, ipsae_d0res_tb),
        "ipsae_d0chn_min": min(ipsae_d0chn_bt, ipsae_d0chn_tb),
        "ipsae_d0chn_max": max(ipsae_d0chn_bt, ipsae_d0chn_tb),
        "ipsae_d0dom_min": min(ipsae_d0dom_bt, ipsae_d0dom_tb),
        "ipsae_d0dom_max": max(ipsae_d0dom_bt, ipsae_d0dom_tb),
        "iptm_d0chn_min": min(iptm_d0chn_bt, iptm_d0chn_tb),
        "iptm_d0chn_max": max(iptm_d0chn_bt, iptm_d0chn_tb),
        "iptm_af_min": min(iptm_af_bt, iptm_af_tb),
        "iptm_af_max": max(iptm_af_bt, iptm_af_tb),
        "pdockq": pdockq,
        "pdockq2": pdockq2,
        "lis": lis,
        "n_interface": n_interface,
    }


def _empty_metrics() -> dict[str, float]:
    keys = [
        "iptm",
        "ptm",
        "mean_plddt",
        "ipsae_d0res_min",
        "ipsae_d0res_max",
        "ipsae_d0chn_min",
        "ipsae_d0chn_max",
        "ipsae_d0dom_min",
        "ipsae_d0dom_max",
        "iptm_d0chn_min",
        "iptm_d0chn_max",
        "iptm_af_min",
        "iptm_af_max",
        "pdockq",
        "pdockq2",
        "lis",
        "n_interface",
    ]
    return dict.fromkeys(keys, 0.0)


# ---------------------------------------------------------------------------
# Modal image — Boltz-2 with ColabFold's MSA server for the target chain
# ---------------------------------------------------------------------------

BOLTZ_MODEL_VOLUME = modal.Volume.from_name(BOLTZ_VOLUME_NAME, create_if_missing=True)
BOLTZ_CACHE_DIR = f"/{BOLTZ_VOLUME_NAME}"


def _download_boltz_models() -> None:
    from boltz.main import download_boltz1, download_boltz2

    cache = Path(BOLTZ_CACHE_DIR)
    if not (cache / "boltz1_conf.ckpt").exists():
        download_boltz1(cache)
    if not (cache / "boltz2_conf.ckpt").exists():
        download_boltz2(cache)


boltz_image = (
    modal.Image.debian_slim(python_version="3.11")
    .micromamba()
    .apt_install("wget", "git", "gcc", "g++", "build-essential")
    .pip_install(
        "colabfold[alphafold-minus-jax]@git+https://github.com/sokrypton/ColabFold@acc0bf772f22feb7f887ad132b7313ff415c8a9f"
    )
    .micromamba_install("kalign2=2.04", "hhsuite=3.3.0", channels=["conda-forge", "bioconda"])
    .run_commands(
        'pip install --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html',
        gpu="a100",
    )
    .run_commands("python -m colabfold.download")
    .pip_install(
        "boltz==2.2.1",
        "pyyaml",
        "gemmi",
        "cuequivariance-torch",
        "cuequivariance-ops-torch-cu12",
    )
    .run_function(
        _download_boltz_models,
        gpu="a10g",
        volumes={BOLTZ_CACHE_DIR: BOLTZ_MODEL_VOLUME},
    )
)

app = modal.App(APP_NAME)

RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"
PREDICTOR = "boltz2"


# ---------------------------------------------------------------------------
# Boltz-2 prediction + IPSAE
# ---------------------------------------------------------------------------


@app.function(
    image=boltz_image,
    gpu=GPU,
    timeout=TIMEOUT_MIN * 60,
    max_containers=CONCURRENCY,
    volumes={BOLTZ_CACHE_DIR: BOLTZ_MODEL_VOLUME, RESULTS_DIR: RESULTS_VOLUME},
)
def predict_boltz2(pb_id: str, binder_seq: str, target_seq: str) -> dict:
    """Run Boltz-2 on (target||binder), compute IPSAE from PAE + CIF."""
    safe = _sanitize(pb_id)
    result_path = Path(RESULTS_DIR) / PREDICTOR / f"{safe}.json"
    if result_path.exists():
        try:
            cached = json.loads(result_path.read_text())
            if cached.get("status") == "ok":
                print(f"  {pb_id}: cached, skipping")
                return cached
        except Exception:
            pass

    try:
        import subprocess

        import numpy as np
        import requests.adapters
        import urllib3
        import yaml

        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        _orig_send = requests.adapters.HTTPAdapter.send
        requests.adapters.HTTPAdapter.send = lambda self, request, **kw: _orig_send(
            self, request, **{**kw, "verify": False}
        )

        target_len = len(target_seq)
        binder_len = len(binder_seq)

        in_dir = Path("/tmp/in_boltz")
        out_dir = Path("/tmp/out_boltz")
        in_dir.mkdir(parents=True, exist_ok=True)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Target gets MSA via --use_msa_server. Binder is explicitly empty
        # ("msa: empty") — matches ProteinTyper's `msa_mode: target_only`.
        yaml_data = {
            "version": 1,
            "sequences": [
                {"protein": {"id": "A", "sequence": target_seq}},
                {"protein": {"id": "B", "sequence": binder_seq, "msa": "empty"}},
            ],
        }
        yaml_path = in_dir / f"{safe}.yaml"
        yaml_path.write_text(yaml.dump(yaml_data, sort_keys=False))

        proc = subprocess.run(
            f'boltz predict "{yaml_path}"'
            f' --out_dir "{out_dir}"'
            f' --cache "{BOLTZ_CACHE_DIR}"'
            f" --use_msa_server --seed 42 --write_full_pae",
            shell=True,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            err = (proc.stderr or proc.stdout or "unknown error")[-500:]
            raise RuntimeError(f"boltz exit {proc.returncode}: {err}")

        cif_files = list(out_dir.glob("**/predictions/**/*.cif"))
        npz_files = list(out_dir.glob("**/predictions/**/pae_*.npz"))
        conf_files = list(out_dir.glob("**/predictions/**/confidence_*.json"))

        if not cif_files or not npz_files:
            result = {
                "pb_id": pb_id,
                "predictor": PREDICTOR,
                "status": "failed",
                **_empty_metrics(),
            }
            _save_result(result)
            return result

        cif_path = cif_files[0]
        pae_data = np.load(str(npz_files[0]))
        pae_matrix = pae_data["pae"]
        if pae_matrix.ndim == 3:
            pae_matrix = pae_matrix[0]

        native_iptm = native_ptm = native_plddt = 0.0
        if conf_files:
            try:
                conf_data = json.loads(conf_files[0].read_text())
                native_iptm = float(conf_data.get("iptm", 0.0))
                native_ptm = float(conf_data.get("ptm", 0.0))
                if "plddt" in conf_data:
                    p = conf_data["plddt"]
                    native_plddt = float(np.mean(p)) if isinstance(p, list) else float(p)
                elif "complex_plddt" in conf_data:
                    native_plddt = float(conf_data["complex_plddt"])
            except Exception as e:
                print(f"  warn: parse conf for {pb_id}: {e}")

        scores = compute_ipsae(
            pae_matrix,
            str(cif_path),
            target_len,
            binder_len,
            target_chain="A",
            binder_chain="B",
        )
        scores["iptm"] = native_iptm
        scores["ptm"] = native_ptm
        scores["mean_plddt"] = native_plddt
        result = {"pb_id": pb_id, "predictor": PREDICTOR, "status": "ok", **scores}
        _save_result(result)
        _save_raw(pb_id, pae_matrix, native_iptm, native_ptm, native_plddt, target_len, binder_len)
        _save_structure(cif_path, pb_id)
        return result

    except Exception as e:
        import traceback

        print(f"Boltz-2 failed for {pb_id}: {traceback.format_exc()}")
        result = {
            "pb_id": pb_id,
            "predictor": PREDICTOR,
            "status": f"error: {e}",
            **_empty_metrics(),
        }
        _save_result(result)
        return result


# ---------------------------------------------------------------------------
# Volume helpers
# ---------------------------------------------------------------------------


def _sanitize(name: str) -> str:
    return (
        name.replace("/", "_SLASH_")
        .replace("\\", "_BSLASH_")
        .replace("|", "_PIPE_")
        .replace(" ", "_")
        .replace(",", "_COMMA_")
    )


def _save_result(result: dict) -> None:
    out_dir = Path(RESULTS_DIR) / PREDICTOR
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{_sanitize(result['pb_id'])}.json").write_text(json.dumps(result))
    RESULTS_VOLUME.commit()


def _save_structure(local_path, pb_id: str) -> None:
    import shutil

    struct_dir = Path(RESULTS_DIR) / "structures" / PREDICTOR
    struct_dir.mkdir(parents=True, exist_ok=True)
    dest = struct_dir / f"{_sanitize(pb_id)}{Path(local_path).suffix}"
    shutil.copy2(str(local_path), str(dest))
    RESULTS_VOLUME.commit()


def _save_raw(
    pb_id: str,
    pae_matrix,
    native_iptm: float,
    native_ptm: float,
    native_plddt: float,
    target_len: int,
    binder_len: int,
) -> None:
    import numpy as np

    raw_dir = Path(RESULTS_DIR) / "raw_data" / PREDICTOR
    raw_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        str(raw_dir / f"{_sanitize(pb_id)}.npz"),
        pae=np.asarray(pae_matrix, dtype=np.float32),
        target_len=np.int32(target_len),
        binder_len=np.int32(binder_len),
        native_iptm=np.float32(native_iptm),
        native_ptm=np.float32(native_ptm),
        native_plddt=np.float32(native_plddt),
    )
    RESULTS_VOLUME.commit()


# ---------------------------------------------------------------------------
# Remote orchestrator + local entrypoint
# ---------------------------------------------------------------------------

orchestrator_image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas")


@app.function(
    image=orchestrator_image,
    timeout=24 * 3600,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def run_batch(pb_ids: list[str], seqs: list[str], target_seq: str) -> None:
    import pandas as pd

    RESULTS_VOLUME.reload()
    pred_dir = Path(RESULTS_DIR) / PREDICTOR
    completed: set[str] = set()
    if pred_dir.exists():
        for f in pred_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                d = json.loads(f.read_text())
                if d.get("status") == "ok":
                    completed.add(d["pb_id"])

    pending = [(p, s) for p, s in zip(pb_ids, seqs, strict=True) if p not in completed]
    print(f"Boltz-2: {len(completed)} done, {len(pending)} pending")

    if pending:
        for done, result in enumerate(
            predict_boltz2.map(
                [p[0] for p in pending],
                [p[1] for p in pending],
                [target_seq] * len(pending),
                return_exceptions=True,
            ),
            start=1,
        ):
            if isinstance(result, Exception):
                print(f"  exception: {result}")
            elif done % 25 == 0:
                print(f"  {done}/{len(pending)} done")

    RESULTS_VOLUME.reload()
    rows = []
    if pred_dir.exists():
        for f in pred_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                rows.append(json.loads(f.read_text()))
    df = pd.DataFrame(rows)
    csv_path = Path(RESULTS_DIR) / f"{PREDICTOR}_summary.csv"
    df.to_csv(csv_path, index=False)
    RESULTS_VOLUME.commit()
    ok = int((df["status"] == "ok").sum()) if "status" in df.columns else 0
    print(f"Boltz-2: {ok}/{len(df)} succeeded.")


@app.local_entrypoint()
def main(
    designs_csv: str = "./data/designs.csv",
    target_fasta: str = "./data/target/rbx1.fasta",
    limit: int | None = None,
    download: bool = False,
    retry_failed: bool = False,
) -> None:
    """Trigger a batch on Modal, or pull results to local disk.

    With ``--detach`` the batch survives client disconnects (Modal best
    practice for multi-hour jobs).
    """
    import pandas as pd

    if download:
        _download_to_local(designs_csv)
        return

    if retry_failed:
        _clear_failed()
        return

    target_seq = _read_fasta(Path(target_fasta))
    if target_seq != RBX1_TARGET_SEQ:
        raise RuntimeError(
            f"data/target/rbx1.fasta sequence does not match RBX1_TARGET_SEQ in this script "
            f"(file={len(target_seq)} aa, script={len(RBX1_TARGET_SEQ)} aa)"
        )

    df = pd.read_csv(designs_csv)
    df = df[df["sequence"].notna() & (df["sequence"].str.len() > 0)].copy()
    if limit:
        df = df.head(limit)

    pb_ids = df["pb_id"].tolist()
    seqs = df["sequence"].tolist()
    print(f"Triggering Boltz-2 batch: {len(pb_ids)} designs (target=RBX1, 108 aa)")
    # Use .remote() + --detach. .spawn() only schedules on `modal deploy`'d
    # apps; under `modal run` (ephemeral), spawn returns immediately but the
    # call never actually runs. The `--detach` flag preserves the LAST
    # .remote() call across local-CLI disconnect.
    run_batch.remote(pb_ids, seqs, target_seq)
    print("Done. Pull results with `--download`.")


def _read_fasta(path: Path) -> str:
    """Return the first sequence in the FASTA file (no header), upper-case."""
    seq_parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">") or not line.strip():
            if seq_parts:
                break
            continue
        seq_parts.append(line.strip())
    return "".join(seq_parts).upper()


def _download_to_local(designs_csv: str) -> None:
    """Pull all per-pb_id JSON + CIF files into ./data/{metrics,structures}/."""
    import pandas as pd

    df = pd.read_csv(designs_csv)
    expected_ids = set(df["pb_id"].astype(str).tolist())

    json_out = Path("./data/metrics") / PREDICTOR
    cif_out = Path("./data/structures") / PREDICTOR
    json_out.mkdir(parents=True, exist_ok=True)
    cif_out.mkdir(parents=True, exist_ok=True)

    n_json = n_cif = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(PREDICTOR):
            if entry.path.endswith(".json"):
                payload = b"".join(RESULTS_VOLUME.read_file(entry.path))
                data = json.loads(payload)
                pb_id = data.get("pb_id") or Path(entry.path).stem
                if pb_id not in expected_ids:
                    continue
                (json_out / f"{pb_id}.json").write_bytes(payload)
                n_json += 1
    except Exception as e:
        print(f"  no metric JSONs to download: {e}")

    try:
        for entry in RESULTS_VOLUME.iterdir(f"structures/{PREDICTOR}"):
            if entry.path.endswith(".cif"):
                pb_id = Path(entry.path).stem
                if pb_id not in expected_ids:
                    continue
                payload = b"".join(RESULTS_VOLUME.read_file(entry.path))
                (cif_out / f"{pb_id}.cif").write_bytes(payload)
                n_cif += 1
    except Exception as e:
        print(f"  no CIFs to download: {e}")

    print(f"Downloaded {n_json} JSON + {n_cif} CIF for {PREDICTOR}.")


def _clear_failed() -> None:
    cleared = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(PREDICTOR):
            if not entry.path.endswith(".json"):
                continue
            data = json.loads(b"".join(RESULTS_VOLUME.read_file(entry.path)))
            status = str(data.get("status", ""))
            if status.startswith("error") or status.startswith("failed"):
                RESULTS_VOLUME.remove_file(entry.path)
                cleared += 1
    except Exception:
        pass
    print(f"  cleared {cleared} failed results")
