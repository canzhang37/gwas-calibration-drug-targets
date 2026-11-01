"""K2: PLINK 2 vs OLS. K3: REGENIE w/o step 1 vs PCs."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

W = Path(__file__).resolve().parents[1] / "data" / "work"
rec = [json.loads(x) for x in open(W / "summary.jsonl")]
kept = {(r["method"], r["trait"]): r["keep"] for r in rec}

raw = pd.read_csv(W / "keep200.raw", sep=r"\s+")
raw.columns = [c.rsplit("_", 1)[0] if i >= 6 else c for i, c in enumerate(raw.columns)]
iid = raw["IID"].astype(str)
ph = pd.read_csv(W / "pheno_q.txt", sep="\t").set_index("IID").loc[iid]
pcs = pd.read_csv(W / "pcs.eigenvec", sep="\t").rename(columns={"#IID": "IID"}).set_index("IID").loc[iid]
PC = pcs[[f"PC{i}" for i in range(1, 11)]].to_numpy()

rng = np.random.default_rng(2027)
snps = list(raw.columns[6:])
traits = [c for c in ph.columns if c.startswith("Q")]
pairs = [(snps[i], traits[j]) for i, j in zip(rng.integers(0, len(snps), 200), rng.integers(0, len(traits), 200))]


def ols_p(y, g, Z):
    X = np.column_stack([np.ones_like(g), g] + ([Z] if Z is not None else []))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    df = len(y) - X.shape[1]
    s2 = res @ res / df
    se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
    return 2 * stats.t.sf(abs(beta[1] / se), df)


gap = {"M1_q": 0.0, "M2_q": 0.0}
for s, t in pairs:
    y, g = ph[t].to_numpy(), raw[s].to_numpy(dtype=float)
    for m, Z in (("M1_q", None), ("M2_q", PC)):
        mine = -np.log10(ols_p(y, g, Z))
        theirs = -np.log10(float(kept[(m, t)][s]))
        gap[m] = max(gap[m], abs(mine - theirs))
print("K2 largest |log10 p| gap:", {k: round(v, 6) for k, v in gap.items()},
      "PASS" if max(gap.values()) < 0.01 else "FAIL")

a, b = [], []
for t in traits:
    for s, v in kept[("K3_q", t)].items():
        if v != "NA" and kept[("M2_q", t)].get(s, "NA") != "NA":
            a.append(float(v))
            b.append(-np.log10(float(kept[("M2_q", t)][s])))
r = np.corrcoef(a, b)[0, 1]
print(f"K3 correlation of -log10 p (REGENIE --ignore-pred vs M2), {len(a)} pairs: {r:.6f}",
      "PASS" if r > 0.999 else "FAIL")
