"""One JSON line per trait from one result batch.

Stdlib only: runs in WSL beside the files,
which can then be deleted.
usage: a4_summarize.py METHOD FORMAT OUT FILE...
"""
import json
import os
import re
import sys
from pathlib import Path
from statistics import NormalDist

W = Path(os.environ.get("O_WORK", Path.home() / "o_data" / "work"))
ALPHAS = (1e-4, 1e-6, 5e-8)
ND = NormalDist()


def load_sets():
    test = set(open(W / "testset.snplist").read().split())
    maf = {}
    with open(W / "base.afreq") as f:
        h = f.readline().lstrip("#").split()
        i_id, i_f = h.index("ID"), h.index("ALT_FREQS")
        for ln in f:
            v = ln.split()
            if v[i_id] in test:
                a = float(v[i_f])
                maf[v[i_id]] = min(a, 1 - a)
    large = {}
    with open(W / "truth.csv") as f:
        h = f.readline().strip().split(",")
        for ln in f:
            d = dict(zip(h, ln.strip().split(",")))
            large.setdefault(int(d["rep"]), {})[d["id"]] = float(d["q"])
    keep = set(open(W / "keep200.txt").read().split())
    return maf, large, keep


def chisq(p):
    if p is None or p <= 0:
        return None
    return ND.inv_cdf(1 - p / 2) ** 2 if p < 1 else 0.0


def trait_of(path):
    return re.findall(r"(B\d\d_\d\d|Q\d\d)", Path(path).name)[-1]


def main():
    method, fmt, out = sys.argv[1], sys.argv[2], sys.argv[3]
    maf, large, keep = load_sets()
    with open(out, "a") as fo:
        for path in sys.argv[4:]:
            trait = trait_of(path)
            rep = int(trait[-2:])
            bins = {"low": [], "common": []}
            na = {"low": 0, "common": 0}
            causal, kept = {}, {}
            with open(path) as f:
                h = f.readline().lstrip("#").split()
                i_id = h.index("ID")
                i_test = h.index("TEST")
                i_p = h.index("P") if fmt == "plink" else h.index("LOG10P")
                for ln in f:
                    v = ln.split()
                    if v[i_test] != "ADD":
                        continue
                    sid, raw = v[i_id], v[i_p]
                    if raw == "NA":
                        p = None
                    else:
                        p = float(raw) if fmt == "plink" else 10 ** (-float(raw))
                    if sid in large.get(rep, {}):
                        causal[sid] = p
                    if sid in keep:
                        kept[sid] = raw
                    if sid in maf:
                        b = "low" if maf[sid] < 0.05 else "common"
                        if p is None:
                            na[b] += 1
                        else:
                            bins[b].append(p)
            rec = {"method": method, "trait": trait, "rep": rep, "bins": {}, "causal": causal, "keep": kept}
            for b, ps in bins.items():
                ps.sort()
                n_ok = len(ps)
                med = ps[n_ok // 2] if n_ok else None
                rec["bins"][b] = {"n": n_ok + na[b], "na": na[b], "median_p": med,
                                  "lambda": chisq(med) / 0.4549364 if med else None,
                                  "rej": {str(a): sum(1 for p in ps if p < a) for a in ALPHAS}}
            fo.write(json.dumps(rec) + "\n")
            print(method, trait, {b: rec["bins"][b]["rej"]["0.0001"] for b in bins}, flush=True)


if __name__ == "__main__":
    main()
