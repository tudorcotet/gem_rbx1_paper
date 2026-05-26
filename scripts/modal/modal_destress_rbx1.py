# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "modal>=1.0",
#     "pandas",
# ]
# ///
"""DE-STRESS (Rosetta + EvoEF2 + BUDEFF) stability scoring for RBX1 designs.

DE-STRESS (https://pragmaticproteindesign.bio.ed.ac.uk/de-stress/) is the
Wood lab's stability scorer, wrapping PyRosetta total_score, EvoEF2,
BUDE force field, and biophysical descriptors. AGGRESCAN3D is Python 2.7
only (skipped); DFIRE2 only ships inside LightDock (too heavy).

For RBX1 we score every available structure for every design:

* ``data/structures/esmfold/{pb_id}.cif`` — binder monomer (ESMFold)
* ``data/structures/boltz2/{pb_id}.cif``  — RBX1+binder complex (Boltz-2)
* ``data/structures/chai/{pb_id}.cif``    — RBX1+binder complex (Chai-1)
* ``data/structures/protenix/{pb_id}.cif``— RBX1+binder complex (Protenix)
* ``data/structures/af2m/{pb_id}.cif``    — RBX1+binder complex (AF2-Multimer)

Missing structures are skipped silently. With 322 designs and 5 predictors
this is up to 1,610 scoring calls.

Outputs land on Modal Volume ``rbx1-rerun-results``:

* ``destress/{pb_id}.json`` — one entry per source predictor (merged)

Usage (Modal CLI auth via ``modal token set ...`` or env vars)::

    cd <repo_root>
    modal run --detach scripts/modal/modal_destress_rbx1.py

    # Download merged per-pb_id JSONs into data/metrics/destress/
    modal run scripts/modal/modal_destress_rbx1.py --download

    # Score a single predictor
    modal run --detach scripts/modal/modal_destress_rbx1.py --predictors boltz2

    # Tune concurrency / per-call timeout (minutes)
    CONCURRENCY=20 TIMEOUT=45 modal run --detach scripts/modal/modal_destress_rbx1.py
"""

from __future__ import annotations

import contextlib
import json
import os
from pathlib import Path

import modal

TIMEOUT_MIN = int(os.environ.get("TIMEOUT", 30))
CONCURRENCY = int(os.environ.get("CONCURRENCY", 10))

# Modal app / volume names -- overridable so external users can run with
# their own workspace naming without editing source.
APP_NAME = os.environ.get("MODAL_APP_NAME", "rbx1-destress")
RESULTS_VOLUME_NAME = os.environ.get("MODAL_RESULTS_VOLUME", "rbx1-rerun-results")

# The five structure sources we iterate over. ESMFold is monomer-only;
# the other four are complexes (chain A = RBX1 target, chain B = binder).
PREDICTORS: tuple[str, ...] = ("esmfold", "boltz2", "chai", "protenix", "af2m")

# ---------------------------------------------------------------------------
# Modal image -- EvoEF2 + BuDEff + (optional) PyRosetta
# ---------------------------------------------------------------------------
# PyRosetta requires a licence and pyrosettacolabsetup downloads a wheel at
# build time, which breaks in Modal's sandboxed image builder (SystemExit).
# We use a conda-based install via micromamba. If it fails the image
# still builds -- Rosetta metrics will be None at runtime.
# ---------------------------------------------------------------------------


def _build_evoef2() -> None:
    """Clone and compile EvoEF2 from source."""
    import subprocess

    subprocess.run(
        "git clone https://github.com/tommyhuangthu/EvoEF2.git /opt/evoef2"
        " && cd /opt/evoef2 && g++ -O3 -o EvoEF2 src/*.cpp -I src/",
        shell=True,
        check=True,
    )


destress_image = (
    modal.Image.micromamba(python_version="3.11")
    .apt_install("git", "wget", "gcc", "g++", "dssp")
    .micromamba_install(
        "pyrosetta",
        channels=["https://conda.graylab.jhu.edu", "conda-forge"],
    )
    .pip_install(
        "biopython>=1.80", "pandas",
    )
    # BuDEff/AMPAL require Cython + complex C build chain.
    # Install separately with Cython pre-built to avoid build failures.
    .pip_install("Cython")
    .run_commands("pip install ampal BUDEFF || echo 'BuDEff install failed, will skip BuDEff scoring'")
    .run_function(_build_evoef2)
)

app = modal.App(APP_NAME)

RESULTS_VOLUME = modal.Volume.from_name(RESULTS_VOLUME_NAME, create_if_missing=True)
RESULTS_DIR = f"/{RESULTS_VOLUME_NAME}"
PREDICTOR_OUT = "destress"


# ---------------------------------------------------------------------------
# Per-(pb_id, predictor) DE-STRESS scoring -- engine ported verbatim
# ---------------------------------------------------------------------------


@app.function(
    image=destress_image,
    cpu=4,
    timeout=TIMEOUT_MIN * 60,
    max_containers=CONCURRENCY,
    volumes={RESULTS_DIR: RESULTS_VOLUME},
)
def score_destress(pb_id: str, predictor: str) -> dict:
    """Compute Rosetta, EvoEF2, BuDEff, and biophysical metrics for one structure.

    Reads ``{RESULTS_DIR}/structures/{predictor}/{pb_id}.{cif,pdb}`` from
    the shared Modal Volume. Returns ``{"status": "no_structure", ...}`` if
    no structure exists. Results are NOT cached at this level -- the per-
    ``pb_id`` JSON aggregation happens in the orchestrator.
    """
    import shutil
    import subprocess
    import tempfile

    safe = _sanitize(pb_id)
    empty = _empty_destress()

    try:
        RESULTS_VOLUME.reload()
        struct_dir = Path(RESULTS_DIR) / "structures" / predictor

        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path: Path | None = None
            for ext in (".pdb", ".cif"):
                src = struct_dir / f"{safe}{ext}"
                if src.exists():
                    dest = Path(tmpdir) / f"{safe}{ext}"
                    shutil.copy2(str(src), str(dest))
                    # Rosetta needs PDB; convert CIF -> PDB.
                    if ext == ".cif":
                        from Bio.PDB import MMCIFParser, PDBIO

                        parser = MMCIFParser(QUIET=True)
                        structure = parser.get_structure("s", str(dest))
                        pdb_dest = Path(tmpdir) / f"{safe}.pdb"
                        io = PDBIO()
                        io.set_structure(structure)
                        io.save(str(pdb_dest))
                        pdb_path = pdb_dest
                    else:
                        pdb_path = dest
                    break

            if pdb_path is None or not pdb_path.exists():
                return {
                    "pb_id": pb_id,
                    "predictor": predictor,
                    "status": "no_structure",
                    **empty,
                }

            metrics: dict = {}

            # --- Rosetta scoring ---
            try:
                import pyrosetta as pr

                pr.init(
                    extra_options=(
                        "-ignore_unrecognized_res -ignore_zero_occupancy -mute all "
                        "-corrections::beta_nov16 true"
                    )
                )

                pose = pr.pose_from_pdb(str(pdb_path))
                sfxn = pr.get_fa_scorefxn()
                total_score = sfxn(pose)
                n_res = pose.total_residue()

                metrics["rosetta_total_per_aa"] = total_score / n_res if n_res > 0 else 0.0

                # Extract individual Rosetta energy terms
                energies = pose.energies()
                terms = {
                    "fa_atr": pr.rosetta.core.scoring.ScoreType.fa_atr,
                    "fa_rep": pr.rosetta.core.scoring.ScoreType.fa_rep,
                    "fa_sol": pr.rosetta.core.scoring.ScoreType.fa_sol,
                    "fa_elec": pr.rosetta.core.scoring.ScoreType.fa_elec,
                    "hbond_sr_bb": pr.rosetta.core.scoring.ScoreType.hbond_sr_bb,
                    "hbond_lr_bb": pr.rosetta.core.scoring.ScoreType.hbond_lr_bb,
                    "hbond_bb_sc": pr.rosetta.core.scoring.ScoreType.hbond_bb_sc,
                    "hbond_sc": pr.rosetta.core.scoring.ScoreType.hbond_sc,
                    "rama_prepro": pr.rosetta.core.scoring.ScoreType.rama_prepro,
                    "fa_dun": pr.rosetta.core.scoring.ScoreType.fa_dun,
                }
                for term_name, term_type in terms.items():
                    term_total = sum(
                        energies.residue_total_energies(i)[term_type]
                        for i in range(1, n_res + 1)
                    )
                    metrics[f"rosetta_{term_name}_per_aa"] = (
                        term_total / n_res if n_res > 0 else 0.0
                    )

            except Exception as e:
                print(f"  Rosetta failed for {pb_id}/{predictor}: {e}")
                metrics["rosetta_total_per_aa"] = None

            # --- EvoEF2 scoring ---
            try:
                evoef2_bin = "/opt/evoef2/EvoEF2"
                if Path(evoef2_bin).exists():
                    result_ef = subprocess.run(
                        f'{evoef2_bin} --command=ComputeStability --pdb="{pdb_path}"',
                        shell=True,
                        capture_output=True,
                        text=True,
                        timeout=120,
                        cwd=str(pdb_path.parent),
                    )
                    # Parse EvoEF2 output
                    for line in result_ef.stdout.split("\n"):
                        if "Total" in line and "=" in line:
                            try:
                                val = float(line.split("=")[-1].strip())
                                # Get residue count from PDB
                                from Bio.PDB import PDBParser

                                parser = PDBParser(QUIET=True)
                                structure = parser.get_structure("s", str(pdb_path))
                                n_res_bio = sum(1 for _ in structure.get_residues())
                                metrics["evoef2_total_per_aa"] = (
                                    val / n_res_bio if n_res_bio > 0 else 0.0
                                )
                            except (ValueError, Exception):
                                pass
                            break
                else:
                    print("  EvoEF2 binary not found")

            except Exception as e:
                print(f"  EvoEF2 failed for {pb_id}/{predictor}: {e}")

            # --- BuDEff scoring (BUDE force field) ---
            try:
                import ampal
                import budeff

                assembly = ampal.load_pdb(str(pdb_path))
                if hasattr(assembly, "__iter__") and not isinstance(assembly, ampal.Assembly):
                    # load_pdb can return list of assemblies; take first
                    assembly = list(assembly)[0] if assembly else None

                if assembly is not None:
                    # Ensure assembly is an Assembly (not AmpalContainer)
                    if isinstance(assembly, ampal.AmpalContainer):
                        assembly = assembly[0] if len(assembly) > 0 else None

                if assembly is not None and isinstance(assembly, ampal.Assembly):
                    buff_score = budeff.get_internal_energy(assembly)
                    n_res_buff = sum(1 for _ in assembly.get_monomers())
                    if n_res_buff > 0:
                        metrics["budeff_total_per_aa"] = buff_score.total_energy / n_res_buff
                        metrics["budeff_steric_per_aa"] = buff_score.steric / n_res_buff
                        metrics["budeff_desolvation_per_aa"] = (
                            buff_score.desolvation / n_res_buff
                        )
                        metrics["budeff_charge_per_aa"] = buff_score.charge / n_res_buff

                    # Also compute interaction energy between chains if 2+ chains
                    chains = list(assembly)
                    if len(chains) >= 2:
                        inter_score = budeff.get_interaction_energy(chains)
                        metrics["budeff_interaction_total"] = inter_score.total_energy
                        metrics["budeff_interaction_steric"] = inter_score.steric
                        metrics["budeff_interaction_desolvation"] = inter_score.desolvation
                        metrics["budeff_interaction_charge"] = inter_score.charge

            except ImportError:
                print("  BuDEff/AMPAL not available, skipping")
            except Exception as e:
                print(f"  BuDEff failed for {pb_id}/{predictor}: {e}")

            # --- Biophysical properties (same as DE-STRESS) ---
            try:
                from Bio.PDB import PDBParser
                from Bio.PDB.Polypeptide import protein_letters_3to1
                from Bio.SeqUtils.ProtParam import ProteinAnalysis

                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("s", str(pdb_path))

                # Extract sequence from structure
                seq = ""
                for residue in structure.get_residues():
                    resname = residue.get_resname().strip().upper()
                    if resname in protein_letters_3to1:
                        seq += protein_letters_3to1[resname]

                if seq:
                    pa = ProteinAnalysis(seq)
                    metrics["isoelectric_point"] = pa.isoelectric_point()
                    metrics["molecular_weight"] = pa.molecular_weight()
                    metrics["num_residues"] = len(seq)
                    with contextlib.suppress(Exception):
                        metrics["instability_index"] = pa.instability_index()
                    with contextlib.suppress(Exception):
                        metrics["gravy"] = pa.gravy()

            except Exception as e:
                print(f"  Biophysical props failed for {pb_id}/{predictor}: {e}")

            return {
                "pb_id": pb_id,
                "predictor": predictor,
                "status": "ok",
                **{**empty, **metrics},
            }

    except Exception as e:
        import traceback

        print(f"DE-STRESS failed for {pb_id}/{predictor}: {traceback.format_exc()}")
        return {
            "pb_id": pb_id,
            "predictor": predictor,
            "status": f"error: {e}",
            **empty,
        }


def _empty_destress() -> dict:
    return {
        "rosetta_total_per_aa": None,
        "rosetta_fa_atr_per_aa": None,
        "rosetta_fa_rep_per_aa": None,
        "rosetta_fa_sol_per_aa": None,
        "rosetta_fa_elec_per_aa": None,
        "rosetta_hbond_sr_bb_per_aa": None,
        "rosetta_hbond_lr_bb_per_aa": None,
        "rosetta_hbond_bb_sc_per_aa": None,
        "rosetta_hbond_sc_per_aa": None,
        "rosetta_rama_prepro_per_aa": None,
        "rosetta_fa_dun_per_aa": None,
        "evoef2_total_per_aa": None,
        "budeff_total_per_aa": None,
        "budeff_steric_per_aa": None,
        "budeff_desolvation_per_aa": None,
        "budeff_charge_per_aa": None,
        "budeff_interaction_total": None,
        "budeff_interaction_steric": None,
        "budeff_interaction_desolvation": None,
        "budeff_interaction_charge": None,
        "isoelectric_point": None,
        "molecular_weight": None,
        "num_residues": None,
        "instability_index": None,
        "gravy": None,
    }


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


def _save_pb_id(pb_id: str, by_predictor: dict[str, dict]) -> None:
    """Write merged per-pb_id JSON: {predictor: {metrics...}}."""
    out_dir = Path(RESULTS_DIR) / PREDICTOR_OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {"pb_id": pb_id, "results": by_predictor}
    (out_dir / f"{_sanitize(pb_id)}.json").write_text(json.dumps(payload))
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
def run_batch(jobs: list[tuple[str, str]]) -> None:
    """Fan out (pb_id, predictor) jobs and aggregate per-pb_id JSON.

    A job is ``(pb_id, predictor)``. Results are grouped by ``pb_id`` and
    written to ``destress/{pb_id}.json`` as ``{"results": {predictor: {...}}}``.
    """
    import pandas as pd

    RESULTS_VOLUME.reload()
    dest_dir = Path(RESULTS_DIR) / PREDICTOR_OUT
    dest_dir.mkdir(parents=True, exist_ok=True)

    # Load any pre-existing per-pb_id aggregates so we can extend rather
    # than overwrite (and skip predictors already scored successfully).
    aggregates: dict[str, dict[str, dict]] = {}
    if dest_dir.exists():
        for f in dest_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                d = json.loads(f.read_text())
                pid = d.get("pb_id") or f.stem
                aggregates[pid] = d.get("results", {}) or {}

    pending: list[tuple[str, str]] = []
    for pb_id, predictor in jobs:
        existing = aggregates.get(pb_id, {}).get(predictor, {})
        if existing.get("status") == "ok":
            continue
        pending.append((pb_id, predictor))

    print(
        f"DE-STRESS: {len(jobs) - len(pending)} cached, {len(pending)} pending "
        f"({len({p for p, _ in pending})} unique pb_ids)"
    )

    if pending:
        for done, result in enumerate(
            score_destress.map(
                [j[0] for j in pending],
                [j[1] for j in pending],
                return_exceptions=True,
                wrap_returned_exceptions=False,
            ),
            start=1,
        ):
            if isinstance(result, Exception):
                print(f"  exception: {result}")
                continue
            pid = result["pb_id"]
            pred = result["predictor"]
            aggregates.setdefault(pid, {})[pred] = result
            # Persist incrementally so partial progress survives.
            _save_pb_id(pid, aggregates[pid])
            if done % 25 == 0:
                print(f"  {done}/{len(pending)} done")

    # Write a flat summary CSV alongside the per-pb_id JSONs.
    RESULTS_VOLUME.reload()
    rows: list[dict] = []
    if dest_dir.exists():
        for f in dest_dir.glob("*.json"):
            with contextlib.suppress(Exception):
                payload = json.loads(f.read_text())
                pid = payload.get("pb_id") or f.stem
                for pred, entry in (payload.get("results") or {}).items():
                    rows.append({"pb_id": pid, **entry})
    if rows:
        df = pd.DataFrame(rows)
        csv_path = Path(RESULTS_DIR) / f"{PREDICTOR_OUT}_summary.csv"
        df.to_csv(csv_path, index=False)
        RESULTS_VOLUME.commit()
        ok = int((df["status"] == "ok").sum()) if "status" in df.columns else 0
        print(f"DE-STRESS: {ok}/{len(df)} (pb_id, predictor) entries OK.")


@app.local_entrypoint()
def main(
    designs_csv: str = "./data/designs.csv",
    structures_root: str = "./data/structures",
    predictors: str = ",".join(PREDICTORS),
    limit: int | None = None,
    download: bool = False,
    retry_failed: bool = False,
) -> None:
    """Trigger a DE-STRESS batch on Modal, or pull results to local disk.

    With ``--detach`` the batch survives client disconnects (Modal best
    practice for multi-hour jobs).

    Args:
        designs_csv: Path to ``data/designs.csv`` (322 rows).
        structures_root: Local root that contains ``{predictor}/{pb_id}.cif``.
        predictors: Comma-separated subset of ``esmfold,boltz2,chai,protenix,af2m``.
        limit: Optional cap on number of designs.
        download: If true, pull ``destress/{pb_id}.json`` to local disk.
        retry_failed: If true, delete failed/no_structure entries from the
            volume so they get re-attempted on the next run.
    """
    import pandas as pd

    selected = [p.strip() for p in predictors.split(",") if p.strip()]
    for p in selected:
        if p not in PREDICTORS:
            raise SystemExit(
                f"Unknown predictor {p!r}; valid: {', '.join(PREDICTORS)}"
            )

    if download:
        _download_to_local(designs_csv)
        return

    if retry_failed:
        _clear_failed()
        return

    df = pd.read_csv(designs_csv)
    if "pb_id" not in df.columns:
        raise SystemExit(f"designs csv {designs_csv} missing pb_id column")
    df = df[df["pb_id"].notna()].copy()
    if limit:
        df = df.head(limit)

    pb_ids: list[str] = df["pb_id"].astype(str).tolist()
    print(f"Loaded {len(pb_ids)} pb_ids from {designs_csv}")

    # Discover (pb_id, predictor) pairs that have a local structure file.
    # Missing structures are skipped silently per task spec.
    root = Path(structures_root)
    jobs: list[tuple[str, str]] = []
    uploads: dict[str, list[tuple[Path, str]]] = {}  # predictor -> [(local, remote_rel)]
    found: dict[str, int] = {p: 0 for p in selected}
    for pb_id in pb_ids:
        safe = _sanitize(pb_id)
        for pred in selected:
            local_cif = root / pred / f"{pb_id}.cif"
            local_pdb = root / pred / f"{pb_id}.pdb"
            local = local_cif if local_cif.exists() else (local_pdb if local_pdb.exists() else None)
            if local is None:
                continue
            jobs.append((pb_id, pred))
            uploads.setdefault(pred, []).append(
                (local, f"structures/{pred}/{safe}{local.suffix}")
            )
            found[pred] += 1

    print("Per-predictor coverage:")
    for pred in selected:
        print(f"  {pred:>10s}: {found[pred]:>4d}/{len(pb_ids)} structures")

    if not jobs:
        print("No structures found locally; nothing to do.")
        return

    # Push local structures up to the Modal Volume so the worker can read
    # them. ``batch_upload`` is the bulk-upload API.
    print(f"Syncing {sum(len(v) for v in uploads.values())} structure files to volume...")
    with RESULTS_VOLUME.batch_upload(force=True) as batch:
        for files in uploads.values():
            for local_path, remote_path in files:
                batch.put_file(str(local_path), remote_path)
    print("  upload complete")

    print(f"Triggering DE-STRESS batch: {len(jobs)} (pb_id, predictor) jobs")
    # Use .remote() + --detach. .spawn() only schedules on `modal deploy`'d
    # apps; under `modal run` (ephemeral), spawn returns immediately but the
    # call never actually runs. The `--detach` flag preserves the LAST
    # .remote() call across local-CLI disconnect.
    run_batch.remote(jobs)
    print("Done. Pull results with `--download`.")


def _download_to_local(designs_csv: str) -> None:
    """Pull all per-pb_id JSON into ``./data/metrics/destress/``."""
    import pandas as pd

    df = pd.read_csv(designs_csv)
    expected_ids = set(df["pb_id"].astype(str).tolist())

    json_out = Path("./data/metrics") / PREDICTOR_OUT
    json_out.mkdir(parents=True, exist_ok=True)

    n_json = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(PREDICTOR_OUT):
            if not entry.path.endswith(".json"):
                continue
            payload = b"".join(RESULTS_VOLUME.read_file(entry.path))
            with contextlib.suppress(Exception):
                data = json.loads(payload)
                pb_id = data.get("pb_id") or Path(entry.path).stem
                if pb_id not in expected_ids:
                    continue
                (json_out / f"{pb_id}.json").write_bytes(payload)
                n_json += 1
    except Exception as e:
        print(f"  no metric JSONs to download: {e}")

    print(f"Downloaded {n_json} JSON for {PREDICTOR_OUT}.")


def _clear_failed() -> None:
    """Delete per-(pb_id, predictor) entries whose status is not ``ok``.

    A pb_id with only some predictors failing is rewritten with only its
    OK entries; the failing predictors are dropped so they get re-tried.
    """
    cleared = 0
    rewritten = 0
    try:
        for entry in RESULTS_VOLUME.iterdir(PREDICTOR_OUT):
            if not entry.path.endswith(".json"):
                continue
            payload = b"".join(RESULTS_VOLUME.read_file(entry.path))
            data = json.loads(payload)
            results = data.get("results") or {}
            kept = {
                pred: r for pred, r in results.items() if r.get("status") == "ok"
            }
            if len(kept) == len(results):
                continue
            cleared += len(results) - len(kept)
            if kept:
                data["results"] = kept
                # Rewrite via the worker mount path -- iterdir returns paths
                # relative to the volume root, e.g. "destress/foo.json".
                pb_id = data.get("pb_id") or Path(entry.path).stem
                (Path(RESULTS_DIR) / PREDICTOR_OUT / f"{_sanitize(pb_id)}.json").write_text(
                    json.dumps(data)
                )
                rewritten += 1
            else:
                RESULTS_VOLUME.remove_file(entry.path)
        RESULTS_VOLUME.commit()
    except Exception as e:
        print(f"  retry-failed scan failed: {e}")
    print(f"  cleared {cleared} non-OK entries; rewrote {rewritten} pb_id files")
