"""Part A tables and the P1-P5 verdicts."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
W = ROOT / "data" / "work"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
rec = [json.loads(x) for x in open(W / "summary.jsonl")]
# binary traits come back from MSI
binfile = W / "binary_summary.jsonl"
if binfile.exists():
    rec += [json.loads(x) for x in open(binfile)]
# amendment 1, post hoc: population indicators
if (W / "posthoc.jsonl").exists():
    rec += [json.loads(x) for x in open(W / "posthoc.jsonl")]

rows = []
for r in rec:
    kind = "quantitative" if r["trait"].startswith("Q") else r["trait"][:3]
    for b, v in r["bins"].items():
        rows.append({"method": r["method"], "kind": kind, "rep": r["rep"], "bin": b, "n": v["n"],
                     "na": v["na"], "lambda": v["lambda"], **{f"rej_{a}": c for a, c in v["rej"].items()}})
d = pd.DataFrame(rows)
d = d[d["method"] != "K3_q"]

rng = np.random.default_rng(2027)
summ = []
for (m, k, b), g in d.groupby(["method", "kind", "bin"]):
    g = g.set_index("rep")
    reps = g.index.to_numpy()
    ratio = g["rej_0.0001"].sum() / (g["n"].sum() * 1e-4)
    boot = []
    for _ in range(2000):
        s = g.loc[rng.choice(reps, len(reps))]
        boot.append(s["rej_0.0001"].sum() / (s["n"].sum() * 1e-4))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    summ.append({"method": m.split("_")[0], "kind": k, "maf_bin": b, "traits": len(g),
                 "tests": int(g["n"].sum()), "na_share": g["na"].sum() / g["n"].sum(),
                 "ratio_1e-4": ratio, "lo": lo, "hi": hi, "median_lambda": g["lambda"].median(),
                 "rej_1e-6": int(g["rej_1e-06"].sum()),
                 "exp_1e-6": g["n"].sum() * 1e-6, "rej_5e-8": int(g["rej_5e-08"].sum()),
                 "exp_5e-8": g["n"].sum() * 5e-8})
s = pd.DataFrame(summ)
s.to_csv(OUT / "type1_error.csv", index=False)

truth = pd.read_csv(W / "truth.csv")
pw = []
for r in rec:
    if not r["trait"].startswith("Q") or r["method"] == "K3_q":
        continue
    t = truth[truth["rep"] == r["rep"]]
    for _, x in t.iterrows():
        p = r["causal"].get(x["id"])
        chi = stats.chi2.isf(p, 1) if p not in (None, 0) else np.nan
        pw.append({"method": r["method"].split("_")[0], "rep": r["rep"], "id": x["id"], "q": x["q"],
                   "p": p, "chisq": chi})
pw = pd.DataFrame(pw)
power = pw.groupby(["method", "q"]).agg(power_5e8=("p", lambda v: np.mean([pi is not None and pi < 5e-8 for pi in v])),
                                        mean_chisq=("chisq", "mean"), snps=("p", "size")).reset_index()
power.to_csv(OUT / "power.csv", index=False)


def get(m, k, b, col="ratio_1e-4"):
    return s[(s.method == m) & (s.kind == k) & (s.maf_bin == b)][col].iloc[0]


chi_ratio = pw[pw.method == "M3"]["chisq"].mean() / pw[pw.method == "M2"]["chisq"].mean()
v = []
p1 = all(get("M1", "quantitative", b) > 2 for b in ("low", "common")) and \
    s[(s.method == "M1") & (s.kind == "quantitative")]["median_lambda"].min() > 1.3
v.append(("P1", "M1 quantitative ratio > 2 in both bins and median lambda > 1.3", p1))
v.append(("P2", "M2 quantitative ratio > 1.25 (MAF >= 0.05)", get("M2", "quantitative", "common") > 1.25))
v.append(("P3", "M3 quantitative ratio in [0.8, 1.25], both bins",
          all(0.8 <= get("M3", "quantitative", b) <= 1.25 for b in ("low", "common"))))
if binfile.exists():
    p4 = not (0.5 <= get("M2", "B02", "low") <= 2) and all(0.5 <= get(m, "B02", "low") <= 2 for m in ("M3", "M3s"))
else:
    p4 = "pending (MSI)"
v.append(("P4", "K = 0.02, MAF 0.01-0.05: M2 outside [0.5, 2], M3 and M3s inside", p4))
v.append(("P5", f"M3/M2 mean chi-square at large-effect SNPs in [0.95, 1.10] (observed {chi_ratio:.3f})",
          0.95 <= chi_ratio <= 1.10))
vd = pd.DataFrame(v, columns=["id", "prediction", "met"])
vd.to_csv(OUT / "verdicts_part_a.csv", index=False)

pd.set_option("display.width", 200)
print(s.round(3).to_string(index=False))
print(power.round(3).to_string(index=False))
print(vd.to_string(index=False))
