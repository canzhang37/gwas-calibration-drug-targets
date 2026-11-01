"""Read one range from a remote bgzipped, tabix file.

Fetches only the index and the blocks in range:
a few MB per gene from a 2 GB file.
"""
import gzip
import io
import struct
import urllib.request
import zlib

UA = {"User-Agent": "Mozilla/5.0"}


def _get(url, start=None, end=None):
    h = dict(UA)
    if start is not None:
        h["Range"] = f"bytes={start}-{end}"
    with urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=120) as r:
        return r.read()


def read_tbi(url):
    raw = gzip.decompress(_get(url))
    b = io.BytesIO(raw)
    magic, n_ref = struct.unpack("<4si", b.read(8))
    assert magic == b"TBI\x01"
    fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<7i", b.read(28))
    names = b.read(l_nm).split(b"\x00")[:-1]
    refs = []
    for _ in range(n_ref):
        bins = {}
        (n_bin,) = struct.unpack("<i", b.read(4))
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack("<Ii", b.read(8))
            bins[bin_id] = [struct.unpack("<QQ", b.read(16)) for _ in range(n_chunk)]
        (n_intv,) = struct.unpack("<i", b.read(4))
        ioff = struct.unpack(f"<{n_intv}Q", b.read(8 * n_intv))
        refs.append((bins, ioff))
    return {"names": [n.decode() for n in names], "refs": refs, "skip": skip,
            "meta": chr(meta), "col_seq": col_seq, "col_beg": col_beg}


def reg2bins(beg, end):
    """UCSC binning used by tabix (0-based, half-open)."""
    end -= 1
    out = [0]
    for shift, off in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        out.extend(range(off + (beg >> shift), off + (end >> shift) + 1))
    return out


def _inflate_blocks(data):
    """Inflate the whole BGZF blocks at the start.

    Returns text and compressed bytes used."""
    out, pos = [], 0
    while pos + 18 <= len(data):
        if data[pos:pos + 2] != b"\x1f\x8b":
            break
        bsize = struct.unpack("<H", data[pos + 16:pos + 18])[0] + 1
        if pos + bsize > len(data):
            break
        out.append(zlib.decompress(data[pos + 18:pos + bsize - 8], -15))
        pos += bsize
    return b"".join(out), pos


def header(url):
    text, _ = _inflate_blocks(_get(url, 0, 200000))
    return text.split(b"\n", 1)[0].decode().split("\t")


def fetch(url, idx, chrom, start, end, coff_from=None):
    """Rows on chrom with start <= pos <= end (1-based).

    coff_from: start at this compressed offset,
    not the index chunks (for the check)."""
    rid = idx["names"].index(str(chrom))
    bins, ioff = idx["refs"][rid]
    beg0 = start - 1
    if coff_from is None:
        lin = ioff[min(beg0 >> 14, len(ioff) - 1)] if ioff else 0
        chunks = [c for b in reg2bins(beg0, end) for c in bins.get(b, []) if c[1] > lin]
        if not chunks:
            return []
        vbeg = max(min(c[0] for c in chunks), lin)
        coff, uoff = vbeg >> 16, vbeg & 0xFFFF
    else:
        coff, uoff = coff_from, 0
    rows, cur, buf, first = [], coff, b"", True
    meta = idx["meta"].encode()
    while True:
        want = 1 << 20
        data = _get(url, cur, cur + want - 1)
        text, used = _inflate_blocks(data)
        last = len(data) < want or used == 0
        if first:
            text, first = text[uoff:], False
            if coff_from is not None:
                # a raw block start can fall mid-line
                text = text.split(b"\n", 1)[1]
        buf += text
        lines = buf.split(b"\n")
        # the tail may continue in the next block
        buf = b"" if last else lines.pop()
        done = False
        for ln in lines:
            if not ln or ln.startswith(meta):
                continue
            f = ln.split(b"\t")
            if f[idx["col_seq"] - 1].decode() != str(chrom):
                if rows:
                    done = True
                    break
                continue
            pos = int(f[idx["col_beg"] - 1])
            if pos > end:
                done = True
                break
            if pos >= start:
                rows.append(ln.decode())
        if done or last:
            break
        cur += used
    return rows
