"""Pan-UKB gene windows (+/- 100 kb): LDL, IHD, T2D."""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "py"))
import tabix_range as tr

BASE = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/"
FILES = {"ldl": "biomarkers-30780-both_sexes-irnt.tsv.bgz",
         "ihd": "phecode-411-both_sexes.tsv.bgz",
         "t2d": "phecode-250.2-both_sexes.tsv.bgz"}
FLANK = 100_000
genes = list(csv.DictReader(open(ROOT / "inputs" / "genes.csv")))
out = ROOT / "data" / "panukb"
out.mkdir(parents=True, exist_ok=True)

for trait, fname in FILES.items():
    url = BASE + fname
    idx = tr.read_tbi(url + ".tbi")
    hdr = tr.header(url)
    keep = [i for i, h in enumerate(hdr) if h in ("chr", "pos", "ref", "alt") or h.endswith("_EUR")]
    for g in genes:
        s, e = int(g["start"]) - FLANK, int(g["end"]) + FLANK
        rows = tr.fetch(url, idx, g["chrom"], s, e)
        with open(out / f"{g['gene']}_{trait}.tsv", "w", newline="") as f:
            f.write("\t".join(hdr[i] for i in keep) + "\n")
            for r in rows:
                v = r.split("\t")
                f.write("\t".join(v[i] for i in keep) + "\n")
        print(trait, g["gene"], len(rows))
