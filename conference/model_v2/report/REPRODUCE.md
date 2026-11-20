# Reproduce J11

Run from `/home/dashii/smt-isef` with the runtime versions pinned in INPUT_LOCK.json.
The lock is an explicit report snapshot. Capture refuses to replace it; generate
and verify reject drift. Missing input files raise; hashes are never refreshed.
Initial capture (already performed for this package):

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python conference/model_v2/j11_report.py capture
```

Regenerate and verify:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python conference/model_v2/j11_report.py generate
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python conference/model_v2/j11_report.py verify
```

Inline reads are capped at 2 MiB/file, except the named J7 v3 and J8 v2 metadata
markers (16 MiB/file), and 64 MiB total. Run with a 6 GiB memory ceiling.
For larger inputs, use `sbatch conference/model_v2/report/generate.sbatch`; it is supplied and has
not been submitted by J11. It requeues on preemption and restarts from the
unchanged input lock; a final manifest is the publication checkpoint. No model
weights or genome shards are loaded. Byte-identical PNGs require the pinned
Python, NumPy, PyArrow and Matplotlib versions and the same font environment.
The manifest excludes its own hash; its files include all input citations,
figures, model cards, reproduction instructions and the Slurm template.
Wet-lab files are never opened. Audit documentation is hashed without parsing.
J9 uses read_dataset() and requires the live pointer to still select the locked
generation; pointer drift fails closed and requires an explicit new snapshot.
