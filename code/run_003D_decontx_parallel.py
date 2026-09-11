#!/usr/bin/env python3
import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
JOBS = ROOT / "runs/003D_decontx/jobs.tsv"
VAR = ROOT / "inputs/002B_asset_audit/TAURUS_actual_var_features_full.tsv"
SCRIPT = ROOT / "code/run_003D_decontx_sample.R"
LOGDIR = ROOT / "runs/003D_decontx/logs"
COORD = ROOT / "runs/003D_decontx/coordinator.log"
WORKERS = int(os.environ.get("DECONTX_WORKERS", "8"))
LOGDIR.mkdir(parents=True, exist_ok=True)

with JOBS.open(newline="", encoding="utf-8") as f:
    jobs = list(csv.DictReader(f, delimiter="\t"))

env = os.environ.copy()
env.update({
    "RETICULATE_PYTHON": "/usr/bin/python3",
    "RCPP_PARALLEL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
})

def stamp():
    return datetime.now(timezone.utc).isoformat()

def run_one(job):
    sid = job["sample_id"]
    cmd = [
        "Rscript", str(SCRIPT), sid, str(ROOT / job["matrix_path"]),
        str(ROOT / job["metadata_path"]), str(VAR), str(ROOT / job["outdir"]),
        job["seed"],
    ]
    with (LOGDIR / f"{sid}.stdout.log").open("w") as out, (LOGDIR / f"{sid}.stderr.log").open("w") as err:
        rc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=out, stderr=err).returncode
    return sid, rc

with COORD.open("a", encoding="utf-8") as log:
    log.write(f"{stamp()} START workers={WORKERS} jobs={len(jobs)}\n")
failures = []
done = 0
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    futures = {ex.submit(run_one, j): j["sample_id"] for j in jobs}
    for fut in as_completed(futures):
        sid, rc = fut.result()
        done += 1
        line = f"{stamp()} FINISH {done}/{len(jobs)} {sid} rc={rc}"
        print(line, flush=True)
        with COORD.open("a", encoding="utf-8") as log:
            log.write(line + "\n")
        if rc != 0:
            failures.append(sid)
with COORD.open("a", encoding="utf-8") as log:
    log.write(f"{stamp()} END failures={len(failures)} {','.join(failures)}\n")
if failures:
    print("FAILED_SAMPLES", ",".join(failures), file=sys.stderr)
    sys.exit(1)
