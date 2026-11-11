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

ROWS_PER_BLOCK = 50_000
XZ_PRESET = 9 | lzma.PRESET_EXTREME
XZ_DICT = 16 * 1024 * 1024


class Block(NamedTuple):
    """header_lines without their "\\n"; exactly one of columns / raw_rows is not None."""
    header_lines: List[bytes]
    columns: Optional[List[bytes]]   # m lines, each the n space-joined tokens of one column
    raw_rows: Optional[bytes]        # the verbatim bytes of the rows (fallback block)
    nrows: int


def _transpose(lines: List[bytes]) -> Optional[List[bytes]]:
    """The columns of `lines` as space-joined byte strings, or None if not transposable."""
    joined = b"".join(lines)
    if (not joined.endswith(b"\n") or joined.startswith(b" ") or b"  " in joined
            or b" \n" in joined or b"\n " in joined or b"\n\n" in joined or joined.startswith(b"\n")):
        return None
    m = lines[0].count(b" ") + 1
    if any(ln.count(b" ") != m - 1 for ln in lines):
        return None
    return [b" ".join(col) for col in zip(*(ln[:-1].split(b" ") for ln in lines))]


def encode(plain_path: str, out_path: str, rows_per_block: int = ROWS_PER_BLOCK,
           preset: int = XZ_PRESET, dict_size: int = XZ_DICT) -> dict:
    """Write `out_path` from the table `plain_path`; it appears only after `verify` passes."""
    comp = lzma.LZMACompressor(format=lzma.FORMAT_XZ, filters=[
        {"id": lzma.FILTER_LZMA2, "preset": preset, "dict_size": dict_size}])
    tmp = out_path + ".tmp"
    blocks = row_blocks = rows = 0
    try:
        with open(plain_path, "rb") as src, open(tmp, "wb") as dst:
            hdr = []
            first = src.readline()
            while first.startswith(b"#") and first.endswith(b"\n"):
                hdr.append(first)
                first = src.readline()
            head = b"".join(hdr)
            pending = [first] if first else []
            while True:
                lines = pending + list(itertools.islice(src, rows_per_block - len(pending)))
                pending = []
                if not lines:
                    break
                cols = _transpose(lines)
                if cols is not None:
                    text = b"#@colblock rows=%d cols=%d hdr=%d\n" % (len(lines), len(cols), len(hdr))
                    text += head + b"\n".join(cols) + b"\n"
                else:
                    raw = b"".join(lines)
                    text = b"#@rowblock rows=%d bytes=%d hdr=%d\n" % (len(lines), len(raw), len(hdr))
                    text += head + raw
                    row_blocks += 1
                dst.write(comp.compress(text))
                blocks += 1
                rows += len(lines)
            if blocks == 0 and hdr:  # header-only table
                dst.write(comp.compress(b"#@rowblock rows=0 bytes=0 hdr=%d\n" % len(hdr) + head))
                blocks = row_blocks = 1
            dst.write(comp.flush())
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
            out.write(b"".join(b" ".join(r) + b"\n" for r in zip(*(c.split(b" ") for c in cols))))
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
