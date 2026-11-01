"""50 quantitative and 150 binary traits on 1000G.

In: data/work/pool.raw, all_phase3.psam
Out: pheno_q.txt, pheno_b.txt, truth.csv,
large_snps.txt, k1_variance.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
W = ROOT / "data" / "work"
R = 50
N_POLY, Q_LARGE = 500, [0.01] * 5 + [0.02] * 5
V_POLY, V_ENV, V_NOISE = 0.25, 0.10, 0.50
PREV = {"B20": 0.20, "B05": 0.05, "B02": 0.02}

raw = pd.read_csv(W / "pool.raw", sep=r"\s+", low_memory=False)
iid = raw["IID"].astype(str).to_numpy()
snp_cols = list(raw.columns[6:])
ids = np.array([c.rsplit("_", 1)[0] for c in snp_cols])
G = raw[snp_cols].to_numpy(dtype=np.float32)
G = (G - G.mean(0)) / G.std(0)
n, m = G.shape

psam = pd.read_csv(W / "all_phase3.psam", sep="\t").rename(columns={"#IID": "IID"})
pop = psam.set_index("IID").loc[iid, "Population"].to_numpy()
pop_codes, pop_idx = np.unique(pop, return_inverse=True)

rng = np.random.default_rng(2027)
Y = np.empty((n, R))
truth, k1 = [], []
for r in range(R):
    pick = rng.choice(m, N_POLY + len(Q_LARGE), replace=False)
    poly_i, large_i = pick[:N_POLY], pick[N_POLY:]
    poly = G[:, poly_i] @ rng.normal(size=N_POLY)
    poly *= np.sqrt(V_POLY / poly.var())
    b_large = np.sqrt(Q_LARGE) * rng.choice([-1, 1], len(Q_LARGE))
    large = G[:, large_i] @ b_large
    env = rng.normal(size=len(pop_codes))[pop_idx]
    env *= np.sqrt(V_ENV / env.var())
    noise = rng.normal(0, np.sqrt(V_NOISE), n)
    y = poly + large + env + noise
    Y[:, r] = y
    tot = y.var()
    k1.append({"rep": r + 1, "poly": poly.var() / tot, "large": large.var() / tot,
               "env": env.var() / tot, "noise": noise.var() / tot})
    for j, (s, q, b) in enumerate(zip(ids[large_i], Q_LARGE, b_large)):
        truth.append({"rep": r + 1, "id": s, "role": "large", "q": q, "beta_std": b})

qn = [f"Q{r + 1:02d}" for r in range(R)]
pq = pd.DataFrame(Y, columns=qn)
pq.insert(0, "IID", iid)
pq.insert(0, "FID", iid)
# LF endings, else PLINK 2 reads Q50\r
pq.to_csv(W / "pheno_q.txt", sep="\t", index=False, float_format="%.6f", lineterminator="\n")

cols = {"FID": iid, "IID": iid}
for lab, k in PREV.items():
    n_case = round(k * n)
    for r in range(R):
        cut = np.sort(Y[:, r])[::-1][n_case - 1]
        cols[f"{lab}_{r + 1:02d}"] = (Y[:, r] >= cut).astype(int)
pb = pd.DataFrame(cols)
pb.to_csv(W / "pheno_b.txt", sep="\t", index=False, lineterminator="\n")

t = pd.DataFrame(truth)
t.to_csv(W / "truth.csv", index=False, lineterminator="\n")
t["id"].drop_duplicates().to_csv(W / "large_snps.txt", index=False, header=False, lineterminator="\n")
kd = pd.DataFrame(k1)
kd.to_csv(W / "k1_variance.csv", index=False)
target = {"poly": V_POLY, "large": sum(Q_LARGE), "env": V_ENV, "noise": V_NOISE}
mean = kd[list(target)].mean()
print(mean.round(4).to_dict())
gap = max(abs(mean[k] - v) for k, v in target.items())
print("K1", "PASS" if gap <= 0.02 else "FAIL", f"(largest gap {gap:.4f})")
print("cases per binary trait:", {lab: int(pb[f'{lab}_01'].sum()) for lab in PREV})
