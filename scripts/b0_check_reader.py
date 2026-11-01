"""K5: range reader vs a sequential scan."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py"))
import tabix_range as tr

BASE = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/"
url = BASE + "biomarkers-30780-both_sexes-irnt.tsv.bgz"
idx = tr.read_tbi(url + ".tbi")
hdr = tr.header(url)
print("columns:", len(hdr), hdr[:12])
chrom, start, end = "1", 55_405_221, 55_630_525
rows = tr.fetch(url, idx, chrom, start, end)
pos = [int(r.split("\t")[1]) for r in rows]

rid = idx["names"].index(chrom)
ioff = idx["refs"][rid][1]
early = ioff[max(((start - 1) >> 14) - 3, 0)] >> 16
scan = tr.fetch(url, idx, chrom, start, end, coff_from=early)

ok_sorted = pos == sorted(pos)
ok_range = all(start <= p <= end for p in pos)
ok_same = rows == scan
print(f"rows {len(rows)}, sorted {ok_sorted}, in range {ok_range}, same as scan {ok_same}")
print("K5", "PASS" if ok_sorted and ok_range and ok_same and rows else "FAIL")
