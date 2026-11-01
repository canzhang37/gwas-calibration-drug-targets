"""Smoke check and merge for the MSI binary run."""
import json
import sys
import zipfile
from datetime import datetime
from pathlib import Path

JOBS = [x for x in Path("msi/jobs.txt").read_text().splitlines() if x.strip()]
METHODS = ("M1_b", "M2_b", "M2f_b", "M3_b", "M3s_b")


def load(path):
    return {(r["method"], r["trait"]): r for r in map(json.loads, open(path))}


def smoke():
    mine = {k: v for k, v in load("output_smoke/smoke.jsonl").items() if k[0] in METHODS}
    ref = {k: v for k, v in load("msi/inputs/smoke_laptop.jsonl").items() if k[0] in METHODS}
    assert set(mine) == set(ref), f"methods/traits differ: {sorted(set(ref) ^ set(mine))}"
    worst = 0.0
    for k, r in ref.items():
        for b, v in r["bins"].items():
            w = mine[k]["bins"][b]
            # same inputs: counts agree up to float noise
            assert abs(w["rej"]["0.0001"] - v["rej"]["0.0001"]) <= 1, (k, b, w["rej"], v["rej"])
            assert w["n"] == v["n"], (k, b, w["n"], v["n"])
            worst = max(worst, abs((w["lambda"] or 0) - (v["lambda"] or 0)))
    assert worst < 0.01, f"lambda differs by {worst}"
    print(f"smoke ok: {len(ref)} method-trait results match the laptop (largest lambda gap {worst:.2e})")


def merge():
    status, recs = [], []
    for i, traits in enumerate(JOBS, 1):
        f, done = Path(f"output/task_{i}.jsonl"), Path(f"output/task_{i}.done")
        n_want = len(traits.split(",")) * len(METHODS)
        if not (f.exists() and done.exists()):
            status.append(f"task {i}: not finished")
            continue
        rows = [json.loads(x) for x in open(f)]
        status.append(f"task {i}: {'ok' if len(rows) == n_want else f'incomplete ({len(rows)} of {n_want})'}")
        recs += rows
    Path("output/binary_summary.jsonl").write_text("".join(json.dumps(r) + "\n" for r in recs))
    ok = sum(s.endswith(": ok") for s in status)
    text = (f"binary run collected {datetime.now():%Y-%m-%d %H:%M}\ntasks: {len(JOBS)}, finished: {ok}, "
            f"records: {len(recs)} of {150 * len(METHODS)}\n" + "\n".join(status) + "\n")
    if ok < len(JOBS):
        text += "To rerun: bash msi/run_all_msi.sh rerun\n"
    Path("STATUS.txt").write_text(text)
    print(text)
    name = f"gwas_binary_results_{datetime.now():%Y%m%d_%H%M}.zip"
    with zipfile.ZipFile(name, "w", zipfile.ZIP_DEFLATED) as z:
        for p in (Path("output/binary_summary.jsonl"), Path("STATUS.txt"), Path("output_smoke/smoke.txt")):
            if p.exists():
                z.write(p, p.as_posix())
    print(f"packed {name}")


if __name__ == "__main__":
    {"smoke": smoke, "merge": merge}[sys.argv[1]]()
