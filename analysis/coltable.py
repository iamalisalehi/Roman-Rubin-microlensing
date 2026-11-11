"""The .colxz format: a lossless, column-wise xz container for the event tables.

A production table (`test<tag>.dat`) is one `#` header line naming the columns, then one row per
event, values separated by single spaces.  Stored column by column, the values of one column sit
next to each other and xz compresses them far better than the row-major text; the decoded bytes
are identical to the source.

File layout: ONE xz stream (or several, concatenated: `cat a.colxz b.colxz` is a valid .colxz,
which is how the pipeline merges chunks).  Inside the stream, a sequence of blocks of up to
`rows_per_block` data rows each:

  column block (every data line ends in "\\n", has exactly m single-space-separated, non-empty
  tokens -- no leading/trailing space, no empty line):
      #@colblock rows=<n> cols=<m> hdr=<k>\\n
      <the k leading '#' header lines of the source, each with its "\\n">
      <m lines; line j = the n tokens of column j joined by single spaces>

  row block (fallback for any block that cannot be transposed, e.g. ragged rows, a double space,
  a missing final newline; also the carrier of a header-only table):
      #@rowblock rows=<n> bytes=<b> hdr=<k>\\n
      <the k header lines>
      <the b raw bytes of the rows, verbatim>

The source's header lines are repeated in every block; `decode` writes them once (from the first
block) and raises ValueError if a later block's header differs (chunks of different runs).  A
zero-byte source encodes to a stream with no blocks.  Memory is bounded by one row block.

CLI:  coltable.py encode PLAIN OUT | decode IN (to stdout) | verify PLAIN IN (exit 0/1) | rows IN
"""
import itertools
import lzma
import os
import sys
from typing import BinaryIO, Iterator, List, NamedTuple, Optional

# Chosen by measurement (Task 2 codec scan, 244 MB bulge table): 50k rows beat 100k by 4.5 % in size;
# with dict_size given explicitly presets 6 and 9 build the same encoder (identical bytes), and
# EXTREME gained nothing (+0.03 %), so the plain preset is used. 16 MiB dict: -5.4 % vs 8 MiB at ~220 MB peak.
ROWS_PER_BLOCK = 50_000
XZ_PRESET = 6
XZ_DICT = 16 * 1024 * 1024


class Block(NamedTuple):
    """header_lines without their "\\n"; exactly one of columns / raw_rows is not None."""
    header_lines: List[bytes]
    columns: Optional[List[bytes]]   # m lines, each the n space-joined tokens of one column
    raw_rows: Optional[bytes]        # the verbatim bytes of the rows (fallback block)
    nrows: int


SLICE_ROWS = 2048             # rows tokenised at a time: bounds the per-token bytes objects
GROUP_BYTES = 8 * 1024 * 1024  # transposed bytes held at once: the block is re-read once per column group


def _slice_ok(part: List[bytes], m1: int) -> bool:
    """True if every line of `part` has m1 spaces and only single-space-separated, non-empty tokens."""
    joined = b"".join(part)
    if (not joined.endswith(b"\n") or joined.startswith(b" ") or joined.startswith(b"\n")
            or b"  " in joined or b" \n" in joined or b"\n " in joined or b"\n\n" in joined):
        return False
    return all(ln.count(b" ") == m1 for ln in part)


def _scan_block(src, rows_per_block: int):
    """Read up to rows_per_block lines from `src` without keeping them.

    Returns (nrows, nbytes, transposable, ncols); the file position ends after the block.
    """
    n = nbytes = 0
    ok, m1 = True, None
    while n < rows_per_block:
        part = list(itertools.islice(iter(src.readline, b""), min(SLICE_ROWS, rows_per_block - n)))
        if not part:
            break
        if m1 is None:
            m1 = part[0].count(b" ")
        ok = ok and _slice_ok(part, m1)
        n += len(part)
        nbytes += sum(map(len, part))
    return n, nbytes, ok, (m1 or 0) + 1


def _read_lines(src, start: int, n: int) -> Iterator[List[bytes]]:
    """The n lines starting at byte `start`, SLICE_ROWS at a time."""
    src.seek(start)
    for done in range(0, n, SLICE_ROWS):
        yield list(itertools.islice(iter(src.readline, b""), min(SLICE_ROWS, n - done)))


def _write_block(src, start: int, n: int, nbytes: int, m: int, ok: bool, head: bytes, nhdr: int,
                 comp, dst) -> None:
    """Compress the block of n rows at byte `start` (its text is re-read, never held whole)."""
    if not ok:
        dst.write(comp.compress(b"#@rowblock rows=%d bytes=%d hdr=%d\n" % (n, nbytes, nhdr) + head))
        src.seek(start)
        left = nbytes
        while left:
            chunk = src.read(min(left, 1 << 20))
            dst.write(comp.compress(chunk))
            left -= len(chunk)
        return
    dst.write(comp.compress(b"#@colblock rows=%d cols=%d hdr=%d\n" % (n, m, nhdr) + head))
    width = max(1, min(m, m * GROUP_BYTES // max(nbytes, 1)))  # columns per group
    for g0 in range(0, m, width):
        g1 = min(m, g0 + width)
        parts = [[] for _ in range(g1 - g0)]
        for lines in _read_lines(src, start, n):
            for j, col in enumerate(zip(*(ln[:-1].split(b" ")[g0:g1] for ln in lines))):
                parts[j].append(b" ".join(col))
        for j in range(g1 - g0):
            dst.write(comp.compress(b" ".join(parts[j]) + b"\n"))
            parts[j] = None


def _compress_file(plain_path: str, tmp: str, rows_per_block: int, preset: int, dict_size: int):
    """Stream `plain_path` into the xz file `tmp`; returns (rows, blocks, row_blocks).

    Kept apart from encode() so the compressor's memory (~12x the dictionary) is freed before
    the verify pass starts.
    """
    comp = lzma.LZMACompressor(format=lzma.FORMAT_XZ, filters=[
        {"id": lzma.FILTER_LZMA2, "preset": preset, "dict_size": dict_size}])
    blocks = row_blocks = rows = 0
    with open(plain_path, "rb") as src, open(tmp, "wb") as dst:
        hdr = []
        start = 0
        line = src.readline()
        while line.startswith(b"#") and line.endswith(b"\n"):
            hdr.append(line)
            start += len(line)
            line = src.readline()
        head = b"".join(hdr)
        src.seek(start)
        while True:
            n, nbytes, ok, m = _scan_block(src, rows_per_block)
            if n == 0:
                break
            end = src.tell()
            _write_block(src, start, n, end - start, m, ok, head, len(hdr), comp, dst)
            src.seek(end)
            start = end
            blocks += 1
            rows += n
            row_blocks += not ok
        if blocks == 0 and hdr:  # header-only table
            dst.write(comp.compress(b"#@rowblock rows=0 bytes=0 hdr=%d\n" % len(hdr) + head))
            blocks = row_blocks = 1
        dst.write(comp.flush())
    return rows, blocks, row_blocks


def _release_memory() -> None:
    """Hand freed heap back to the OS (glibc only) so the verify pass starts from a small footprint."""
    try:
        import ctypes
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except (OSError, AttributeError):
        pass


def encode(plain_path: str, out_path: str, rows_per_block: int = ROWS_PER_BLOCK,
           preset: int = XZ_PRESET, dict_size: int = XZ_DICT) -> dict:
    """Write `out_path` from the table `plain_path`; it appears only after `verify` passes."""
    tmp = out_path + ".tmp"
    try:
        rows, blocks, row_blocks = _compress_file(plain_path, tmp, rows_per_block, preset, dict_size)
        _release_memory()
        if not verify(plain_path, tmp):
            raise RuntimeError("colxz verify failed for %s" % plain_path)
        os.replace(tmp, out_path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    return {"rows": rows, "blocks": blocks, "row_blocks": row_blocks,
            "bytes_in": os.path.getsize(plain_path), "bytes_out": os.path.getsize(out_path)}


def _marker(f) -> Optional[dict]:
    line = f.readline()
    if not line:
        return None
    parts = line.split()
    if not parts or parts[0] not in (b"#@colblock", b"#@rowblock"):
        raise ValueError("not a .colxz block marker: %r" % line[:60])
    kv = {k.decode(): int(v) for k, v in (p.split(b"=") for p in parts[1:])}
    kv["kind"] = parts[0][2:].decode()
    return kv


def iter_blocks(path: str) -> Iterator[Block]:
    """The blocks of a .colxz file in order (concatenated streams read as one)."""
    with lzma.open(path, "rb") as f:
        while True:
            mk = _marker(f)
            if mk is None:
                return
            hdr = [f.readline().rstrip(b"\n") for _ in range(mk["hdr"])]
            if mk["kind"] == "colblock":
                cols = [f.readline().rstrip(b"\n") for _ in range(mk["cols"])]
                yield Block(hdr, cols, None, mk["rows"])
            else:
                raw = f.read(mk["bytes"])
                if len(raw) != mk["bytes"]:
                    raise ValueError("truncated row block in %s" % path)
                yield Block(hdr, None, raw, mk["rows"])


def _write_rows(cols: List[bytes], n: int, out: BinaryIO) -> None:
    """Write the n rows of a column block, SLICE_ROWS at a time (bounds the token objects)."""
    rest = list(cols)
    for done in range(0, n, SLICE_ROWS):
        k = min(SLICE_ROWS, n - done)
        toks = []
        for j, c in enumerate(rest):
            pieces = c.split(b" ", k)
            toks.append(pieces[:k])
            rest[j] = pieces[k] if len(pieces) > k else b""
        out.write(b"".join(b" ".join(r) + b"\n" for r in zip(*toks)))


def decode(path: str, out: BinaryIO) -> int:
    """Write the row-major text of `path` to `out`; return the number of data rows."""
    header = None
    rows = 0
    for hdr, cols, raw, n in iter_blocks(path):
        if header is None:
            header = hdr
            out.write(b"".join(h + b"\n" for h in hdr))
        elif hdr != header:
            raise ValueError("%s: chunks have different headers" % path)
        if cols is not None:
            _write_rows(cols, n, out)
        else:
            out.write(raw)
        rows += n
    return rows


class _Compare:
    """A write-only sink that compares what it is fed against a file, byte for byte."""

    def __init__(self, f):
        self.f, self.ok = f, True

    def write(self, data: bytes) -> None:
        if self.ok and self.f.read(len(data)) != data:
            self.ok = False


def verify(plain_path: str, colxz_path: str) -> bool:
    """True if decoding `colxz_path` reproduces `plain_path` byte for byte."""
    with open(plain_path, "rb") as f:
        cmp = _Compare(f)
        decode(colxz_path, cmp)
        return cmp.ok and f.read(1) == b""


def main(argv: List[str]) -> int:
    cmd, args = (argv[0], argv[1:]) if argv else ("", [])
    if cmd == "encode" and len(args) == 2:
        print(encode(*args))
    elif cmd == "decode" and len(args) == 1:
        decode(args[0], sys.stdout.buffer)
    elif cmd == "verify" and len(args) == 2:
        return 0 if verify(*args) else 1
    elif cmd == "rows" and len(args) == 1:
        print(sum(b.nrows for b in iter_blocks(args[0])))
    else:
        print(__doc__.split("CLI:")[1].strip(), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
