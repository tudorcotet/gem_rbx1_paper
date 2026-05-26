"""Modal apps for RBX1 complex prediction.

Each module here is a self-contained ``modal run`` script. They are *not*
intended to be imported — the local driver (``scripts/data/run_rescoring.py``)
shells out to ``modal run`` for fan-out + Volume downloads.

Conventions match the upstream complex re-scoring stack so a single pipeline produces
parity outputs across Boltz-2, Protenix, and Chai-1:

* Target chain ``A`` (RBX1, 108 aa with target MSA).
* Binder chain ``B`` (single-sequence, no MSA / empty MSA).
* Outputs land in a shared Modal Volume ``rbx1-rerun-results`` under
  ``{predictor}/{pb_id}.{cif,json}``.

The IPSAE / pDockQ / LIS scoring helpers are inlined per-file to keep each
Modal image self-contained — copying them across files is intentional and
mirrors ``rescoring stack``.
"""
