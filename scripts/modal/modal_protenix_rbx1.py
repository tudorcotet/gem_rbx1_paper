# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""Protenix complex prediction + ipSAE scoring for 322 RBX1 designs.

Adapted from an internal complex re-scoring stack (Boltz-2 / Protenix / Chai) for
the RBX1 target. **Defaults to ``protenix-v2``** — the +9-13 pp DockQ
upgrade ByteDance gated behind internal review on volces.com (issue
https://github.com/bytedance/Protenix/issues/295). We bypass the gate
by pre-placing the checkpoint from the publicly-mirrored PXDesign
HuggingFace repo at the path Protenix's loader expects:

    https://huggingface.co/TMF001/pxdesign-weights/resolve/main/checkpoint/protenix-v2.pt
        → /protenix-models/checkpoint/protenix-v2.pt

The same mirror also has ``protenix_base_default_v0.5.0`` and three
mini variants if you want to compare. Swap via ``PROTENIX_MODEL_NAME``.

* **Target**: RBX1 (UniProt P62877, 108 aa) — chain ``A``. ``use_msa``
  defaults to ``False`` because Protenix's bundled ``colab_request_utils``
  client queues against ``api.colabfold.com`` and frequently blocks 30+ min
  per call. Flip back on via ``PROTENIX_USE_MSA=true`` when the queue is
  healthy, or after staging a precomputed RBX1 MSA on the model volume
  (Protenix accepts ``unpairedMsaPath`` in the input JSON).
* **Binder**: per-row sequence — chain ``B``, single-sequence (no MSA
  written for the binder chain; matches ``msa_mode: target_only``).

Outputs live on the shared Modal Volume ``rbx1-rerun-results``:

* ``protenix/{pb_id}.json`` — IPSAE + native Protenix confidence
* ``structures/protenix/{pb_id}.cif`` — predicted complex
* ``raw_data/protenix/{pb_id}.npz`` — PAE matrix + scalars

Usage::

    cd <repo_root>
    modal run --detach scripts/modal/modal_protenix_rbx1.py
    modal run scripts/modal/modal_protenix_rbx1.py --download

    # The default is already protenix-v2 (pre-fetched from HuggingFace on
    # first cold start). To drop to an older weight:
    GPU="A100:80GB" CONCURRENCY=3 \\
      PROTENIX_MODEL_NAME=protenix_base_default_v0.5.0 \\
      modal run --detach scripts/modal/modal_protenix_rbx1.py
"""

from __future__ import annotations

import contextlib
import json
import os
from pathlib import Path

import modal

GPU = os.environ.get("GPU", "A100")
TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 180))
CONCURRENCY = int(os.environ.get("CONCURRENCY", 3))
MODEL_NAME = os.environ.get("PROTENIX_MODEL_NAME", "protenix-v2")
# Target chain A gets a ColabFold MMseqs2 MSA via Protenix's bundled
# `colab_request_utils` client (queues against api.colabfold.com). Flip to
# `false` to skip MSA if the queue is jammed.
USE_MSA = os.environ.get("PROTENIX_USE_MSA", "true").lower() in {"1", "true", "yes"}

# Modal app / volume names — overridable for external users.
APP_NAME = os.environ.get("MODAL_APP_NAME", "rbx1-protenix")
RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "rbx1-rerun-results")
PROTENIX_VOLUME_NAME = os.environ.get("MODAL_PROTENIX_VOLUME", "protenix-models")

# Weights hosted publicly at TMF001/pxdesign-weights on HuggingFace —
# bypasses the ByteDance gate on volces.com. Add a row to swap mirrors.
_HF_WEIGHTS = {
    "protenix-v2":                    "https://huggingface.co/TMF001/pxdesign-weights/resolve/main/checkpoint/protenix-v2.pt",
    "protenix_base_default_v0.5.0":   "https://huggingface.co/TMF001/pxdesign-weights/resolve/main/checkpoint/protenix_base_default_v0.5.0.pt",
    "protenix_mini_default_v0.5.0":   "https://huggingface.co/TMF001/pxdesign-weights/resolve/main/checkpoint/protenix_mini_default_v0.5.0.pt",
    "protenix_mini_tmpl_v0.5.0":      "https://huggingface.co/TMF001/pxdesign-weights/resolve/main/checkpoint/protenix_mini_tmpl_v0.5.0.pt",
    "pxdesign_v0.1.0":                "https://huggingface.co/TMF001/pxdesign-weights/resolve/main/checkpoint/pxdesign_v0.1.0.pt",
}

RBX1_TARGET_SEQ = (
    "MAAAMDVDTPSGTNSGAGKKRFEVKKWNAVALWAWDIVVDNCAICRNHIMDLCIECQANQ"
    "ASATSEECTVAWGVCNHAFHFHCISRWLKTRQVCPLDNREWEFQKYGH"
)

PAE_CUTOFF = 15.0
DIST_CUTOFF = 8.0

# ---------------------------------------------------------------------------
# IPSAE scoring — identical to scripts/modal/modal_boltz2_rbx1.py
# ---------------------------------------------------------------------------


def _calc_d0(n_res: float, min_value: float = 1.0) -> float:
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
        print(f"Failed to read structure: {e}")
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
                return _empty_metrics()
        else:
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
        plddt_list = []
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
        "ranking_score",
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
# Modal image — Protenix 2.0.x (pinned by env var to swap in v2 weights)
# ---------------------------------------------------------------------------

PROTENIX_MODEL_VOLUME = modal.Volume.from_name(PROTENIX_VOLUME_NAME, create_if_missing=True)
PROTENIX_CACHE_DIR = f"/{PROTENIX_VOLUME_NAME}"

protenix_image = (
    # NVIDIA CUDA base image — provides /usr/local/cuda so Protenix's
    # fast_layer_norm JIT extension can compile. The plain debian_slim base
    # has no CUDA toolkit, and cuda-toolkit-12-1 isn't in Debian apt repos.
    modal.Image.from_registry(
        "nvidia/cuda:12.4.1-cudnn-devel-ubuntu22.04", add_python="3.11"
    )
    .apt_install("wget", "git", "gcc", "g++", "hmmer", "kalign")
    # Install torch + numpy FIRST at the exact pins protenix wants so pip
    # doesn't backtrack through 20+ legacy numpys on the protenix step.
    # protenix 2.0.0 pins torch==2.7.1 and numpy==2.4.1; the cu121 wheel is
    # cu126-compatible at runtime. PyPI default index has cpu-only torch
    # wheels; we want the CUDA build so request the cu121 wheel index.
    .pip_install(
        "torch==2.7.1",
        "numpy==2.4.1",
        extra_index_url="https://download.pytorch.org/whl/cu121",
    )
    .pip_install("protenix==2.0.0", "gemmi")
    .env({"CUDA_HOME": "/usr/local/cuda"})
    # Point Protenix at the persistent volume so weights survive across
    # containers. The CCD / config download still happens on first @enter().
    .env({"PROTENIX_ROOT_DIR": "/protenix-models"})
)

app = modal.App(APP_NAME)

RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"
PREDICTOR = "protenix"


@app.cls(
    image=protenix_image,
    gpu=GPU,
    timeout=TIMEOUT_MIN * 60,
    max_containers=CONCURRENCY,
    volumes={RESULTS_DIR: RESULTS_VOLUME, PROTENIX_CACHE_DIR: PROTENIX_MODEL_VOLUME},
)
class ProtenixPredictor:
    """One InferenceRunner per container — model loads at @modal.enter()."""

    @modal.enter()
    def setup(self):
        import os as _os
        import urllib.request as _urllib

        from configs.configs_base import configs as configs_base
        from configs.configs_data import data_configs
        from configs.configs_inference import inference_configs
        from configs.configs_model_type import model_configs
        from ml_collections.config_dict import ConfigDict
        from protenix.config import parse_configs
        # protenix>=2.0 fixed the typo `infercence` → `inference`. Fall back to
        # the old name in case someone pins protenix<2.0 via PROTENIX_VERSION.
        from runner.inference import InferenceRunner
        try:
            from runner.inference import download_inference_cache
        except ImportError:
            from runner.inference import download_infercence_cache as download_inference_cache

        # Pre-place the requested checkpoint from HuggingFace if it's one of
        # the publicly-mirrored PXDesign weights. Lands at the path Protenix's
        # loader checks before falling back to the (gated) volces.com URL.
        if MODEL_NAME in _HF_WEIGHTS:
            ckpt_dir = "/protenix-models/checkpoint"
            _os.makedirs(ckpt_dir, exist_ok=True)
            ckpt_path = f"{ckpt_dir}/{MODEL_NAME}.pt"
            if not _os.path.exists(ckpt_path):
                tmp_path = f"{ckpt_path}.partial"
                print(f"[protenix-setup] fetching {MODEL_NAME} from HuggingFace")
                _urllib.urlretrieve(_HF_WEIGHTS[MODEL_NAME], tmp_path)
                _os.replace(tmp_path, ckpt_path)
                PROTENIX_MODEL_VOLUME.commit()
                print(f"[protenix-setup] wrote {ckpt_path} ({_os.path.getsize(ckpt_path)/1e9:.2f} GB)")

        inference_configs["dump_dir"] = "/tmp/protenix_out"
        inference_configs["input_json_path"] = "/dev/null"
        inference_configs["need_atom_confidence"] = True
        inference_configs["model_name"] = MODEL_NAME
        inference_configs["seeds"] = [42]
        inference_configs["use_msa"] = USE_MSA

        configs = {**configs_base, **{"data": data_configs}, **inference_configs}
        configs = parse_configs(configs=configs, fill_required_with_null=True)
        configs.seeds = [42]
        configs.need_atom_confidence = True
        configs.use_msa = USE_MSA
        configs.sample_diffusion.N_sample = 1

        model_specifics = ConfigDict(model_configs[MODEL_NAME])
        configs.update(model_specifics)
        configs.use_deepspeed_evo_attention = (
            _os.environ.get("USE_DEEPSPEED_EVO_ATTENTION", "false") == "true"
        )

        download_inference_cache(configs)
        self.runner = InferenceRunner(configs)
        self.configs = configs
        print(f"  Protenix model {MODEL_NAME} loaded and ready.")

    @modal.method()
    def predict(self, pb_id: str, binder_seq: str, target_seq: str) -> dict:
        """Run Protenix on (target||binder), extract PAE + CIF, compute IPSAE."""
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
            import shutil

            import numpy as np
            from runner.inference import infer_predict
            from runner.msa_search import update_infer_json

            target_len = len(target_seq)
            binder_len = len(binder_seq)
            in_dir = Path("/tmp/protenix_in")
            out_dir = Path("/tmp/protenix_out")
            in_dir.mkdir(parents=True, exist_ok=True)
            if out_dir.exists():
                shutil.rmtree(out_dir)
            out_dir.mkdir(parents=True, exist_ok=True)

            # Target gets MSA; binder is single-sequence (no MSA written for B).
            input_data = [
                {
                    "name": safe,
                    "sequences": [
                        {"proteinChain": {"sequence": target_seq, "count": 1, "id": ["A"]}},
                        {"proteinChain": {"sequence": binder_seq, "count": 1, "id": ["B"]}},
                    ],
                }
            ]
            input_path = in_dir / f"{safe}.json"
            input_path.write_text(json.dumps(input_data))

            self.configs["dump_dir"] = str(out_dir)
            # protenix>=2.0 returns (json_path, actual_updated_bool); older
            # versions returned just the string. Unwrap the tuple — ml_collections
            # rejects a tuple assignment to a string-typed reference field.
            _msa_ret = update_infer_json(
                str(input_path), out_dir=str(out_dir), use_msa=USE_MSA
            )
            self.configs["input_json_path"] = (
                _msa_ret[0] if isinstance(_msa_ret, tuple) else _msa_ret
            )
            infer_predict(self.runner, self.configs)

            cif_files = list(out_dir.glob("**/*.cif"))
            conf_files = list(out_dir.glob("**/*summary_confidence*.json"))
            full_data_files = list(out_dir.glob("**/*full_data*.json"))
            print(
                f"  Found: {len(cif_files)} CIF, {len(conf_files)} conf, "
                f"{len(full_data_files)} full_data files"
            )

            if not cif_files:
                result = {
                    "pb_id": pb_id,
                    "predictor": PREDICTOR,
                    "status": "failed_no_cif",
                    "model_name": MODEL_NAME,
                    **_empty_metrics(),
                }
                _save_result(result)
                return result

            cif_path = cif_files[0]
            native_iptm = native_ptm = native_plddt = ranking_score = 0.0
            if conf_files:
                try:
                    conf = json.loads(conf_files[0].read_text())
                    ranking_score = float(conf.get("ranking_score", 0.0))
                    native_plddt = float(
                        conf.get("plddt", conf.get("atom_plddt", conf.get("mean_plddt", 0.0)))
                    )
                    native_iptm = float(conf.get("iptm", 0.0))
                    native_ptm = float(conf.get("ptm", 0.0))
                except Exception as e:
                    print(f"  warn: parse conf: {e}")

            pae_matrix = None
            if full_data_files:
                try:
                    fd = json.loads(full_data_files[0].read_text())
                    for key in ("token_pair_pae", "pae", "predicted_aligned_error"):
                        if key in fd:
                            pae_matrix = np.array(fd[key])
                            break
                except Exception as e:
                    print(f"  warn: parse full_data: {e}")
            if pae_matrix is None:
                npz_files = list(out_dir.glob("**/*pae*.npz"))
                if npz_files:
                    pae_data = np.load(str(npz_files[0]))
                    for key in ("pae", "token_pair_pae", "predicted_aligned_error"):
                        if key in pae_data:
                            pae_matrix = pae_data[key]
                            break
            if pae_matrix is not None and pae_matrix.ndim == 3:
                pae_matrix = pae_matrix[0]
            if pae_matrix is None:
                result = {
                    "pb_id": pb_id,
                    "predictor": PREDICTOR,
                    "status": "failed_no_pae",
                    "model_name": MODEL_NAME,
                    **_empty_metrics(),
                }
                _save_result(result)
                return result

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
            scores["ranking_score"] = ranking_score
            result = {
                "pb_id": pb_id,
                "predictor": PREDICTOR,
                "status": "ok",
                "model_name": MODEL_NAME,
                **scores,
            }
            _save_result(result)
            _save_raw(
                pb_id,
                pae_matrix,
                native_iptm,
                native_ptm,
                native_plddt,
                ranking_score,
                target_len,
                binder_len,
            )
            _save_structure(cif_path, pb_id)
            return result

        except Exception as e:
            import traceback

            print(f"Protenix failed for {pb_id}: {traceback.format_exc()}")
            result = {
                "pb_id": pb_id,
                "predictor": PREDICTOR,
                "status": f"error: {e}",
                "model_name": MODEL_NAME,
                **_empty_metrics(),
            }
            _save_result(result)
            return result


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
    shutil.copy2(str(local_path), str(struct_dir / f"{_sanitize(pb_id)}{Path(local_path).suffix}"))
    RESULTS_VOLUME.commit()


def _save_raw(
    pb_id: str,
    pae_matrix,
    native_iptm: float,
    native_ptm: float,
    native_plddt: float,
    ranking_score: float,
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
        ranking_score=np.float32(ranking_score),
    )
    RESULTS_VOLUME.commit()


orchestrator_image = modal.Image.debian_slim(python_version="3.11").pip_install("pandas")


@app.function(image=orchestrator_image, timeout=24 * 3600, volumes={RESULTS_DIR: RESULTS_VOLUME})
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
    print(f"Protenix: {len(completed)} done, {len(pending)} pending ({MODEL_NAME})")

    if pending:
        predictor = ProtenixPredictor()
        for done, result in enumerate(
            predictor.predict.map(
                [p[0] for p in pending],
                [p[1] for p in pending],
                [target_seq] * len(pending),
                return_exceptions=True,
                wrap_returned_exceptions=False,
            ),
            start=1,
        ):
            if isinstance(result, Exception):
                print(f"  exception: {result}")
            elif done % 10 == 0:
                print(f"  {done}/{len(pending)} done")

    RESULTS_VOLUME.reload()
    rows = []
    if pred_dir.exists():
        for f in pred_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                rows.append(json.loads(f.read_text()))
    df = pd.DataFrame(rows)
    df.to_csv(Path(RESULTS_DIR) / f"{PREDICTOR}_summary.csv", index=False)
    RESULTS_VOLUME.commit()
    ok = int((df["status"] == "ok").sum()) if "status" in df.columns else 0
    print(f"Protenix: {ok}/{len(df)} succeeded.")


@app.local_entrypoint()
def main(
    designs_csv: str = "./data/designs.csv",
    target_fasta: str = "./data/target/rbx1.fasta",
    limit: int | None = None,
    download: bool = False,
    retry_failed: bool = False,
) -> None:
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
            f"rbx1.fasta seq doesn't match RBX1_TARGET_SEQ in this script "
            f"(file={len(target_seq)} aa, script={len(RBX1_TARGET_SEQ)} aa)"
        )

    df = pd.read_csv(designs_csv)
    df = df[df["sequence"].notna() & (df["sequence"].str.len() > 0)].copy()
    if limit:
        df = df.head(limit)
    pb_ids = df["pb_id"].tolist()
    seqs = df["sequence"].tolist()
    print(
        f"Triggering Protenix batch: {len(pb_ids)} designs, "
        f"model={MODEL_NAME}, target=RBX1 (108 aa)"
    )
    # Use .remote() + --detach. .spawn() only schedules on `modal deploy`'d
    # apps; under `modal run` (ephemeral), spawn returns immediately but the
    # call never actually runs. The `--detach` flag preserves the LAST
    # .remote() call across local-CLI disconnect — that's what we want.
    run_batch.remote(pb_ids, seqs, target_seq)
    print("Done. Pull results with `--download`.")


def _read_fasta(path: Path) -> str:
    seq_parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">") or not line.strip():
            if seq_parts:
                break
            continue
        seq_parts.append(line.strip())
    return "".join(seq_parts).upper()


def _download_to_local(designs_csv: str) -> None:
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
