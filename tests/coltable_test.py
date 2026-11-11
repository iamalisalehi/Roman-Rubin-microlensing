"""Round-trip tests for analysis/coltable.py (the .colxz column-block xz format).

Standalone script, not pytest:  .roman/bin/python tests/coltable_test.py
Scratch files go to the system temp dir (set TMPDIR to redirect).
"""
import io
import os
import random
import resource
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
import coltable as ct  # noqa: E402

SAMPLE = os.path.expanduser("~/.claude/jobs/b0aec636/tmp/kin_sample/rv/bulge/389/test5.dat")
HDR = b"# a b c\n"


def decoded(path):
    buf = io.BytesIO()
    n = ct.decode(path, buf)
    return buf.getvalue(), n


def roundtrip(tmp, name, data, **kw):
    """Encode bytes `data`, check verify() and decode() give them back; return encode's stats."""
    src, out = os.path.join(tmp, name + ".dat"), os.path.join(tmp, name + ".colxz")
    with open(src, "wb") as f:
        f.write(data)
    st = ct.encode(src, out, **kw)
    assert ct.verify(src, out), name
    got, _ = decoded(out)
    assert got == data, name
    assert not os.path.exists(out + ".tmp"), name
    return st, out


RSS_LIMIT_MB = 300


def memory_case(tmp):
    """Encode and decode a wide synthetic table (154 cols x 60,000 rows, 2 blocks) in child
    processes at the default constants; their peak RSS must stay under RSS_LIMIT_MB."""
    rng = random.Random(1)
    pool = ["%.6g" % (rng.random() * 10 ** rng.randint(-3, 4)) for _ in range(30_000)]
    src, out = os.path.join(tmp, "wide.dat"), os.path.join(tmp, "wide.colxz")
    with open(src, "w") as f:
        f.write("# " + " ".join("c%d" % j for j in range(154)) + "\n")
        for _ in range(60_000):
            f.write(" ".join(rng.choices(pool, k=154)) + "\n")
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis", "coltable.py")
    subprocess.run([sys.executable, script, "encode", src, out], check=True, stdout=subprocess.DEVNULL)
    enc = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
    with open(os.devnull, "wb") as null:
        subprocess.run([sys.executable, script, "decode", out], check=True, stdout=null)
    dec = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024  # max over both children
    assert ct.verify(src, out)
    print("wide table: %d -> %d bytes, peak RSS %.0f MB (encode), %.0f MB (max with decode)"
          % (os.path.getsize(src), os.path.getsize(out), enc, dec))
    assert enc <= RSS_LIMIT_MB and dec <= RSS_LIMIT_MB, (enc, dec)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        # first, while this process is still small: a child's ru_maxrss starts from its parent's size at fork
        memory_case(tmp)
        if os.path.exists(SAMPLE):
            data = open(SAMPLE, "rb").read()
            st, _ = roundtrip(tmp, "sample", data, rows_per_block=50_000)
            assert st["blocks"] == 1 and st["row_blocks"] == 0 and st["rows"] == 50_000, st
            assert st["bytes_in"] == len(data) and st["bytes_out"] < len(data) // 4, st
            st, _ = roundtrip(tmp, "sample7", data, rows_per_block=7_001)
            assert st["blocks"] == 8 and st["row_blocks"] == 0, st
        else:
            print("SKIP real-table cases: %s absent" % SAMPLE)

        st, _ = roundtrip(tmp, "hdronly", HDR)
        assert st["rows"] == 0, st
        st, _ = roundtrip(tmp, "one", HDR + b"1 2.5 -1\n")
        assert st["rows"] == 1 and st["row_blocks"] == 0, st
        st, _ = roundtrip(tmp, "empty", b"")
        assert st["rows"] == 0, st
        roundtrip(tmp, "crlf", HDR + b"1 2 3\r\n4 5 6\r\n")  # \r stays inside the last token
        st, _ = roundtrip(tmp, "nohdr", b"1 2\n3 4\n")
        assert st["rows"] == 2 and st["row_blocks"] == 0, st

        base = HDR + b"1 2 3\n4 5 6\n7 8 9\n"
        odd = {
            "tokens": HDR + b"1 2 3\n4 5\n7 8 9\n",
            "dblspace": HDR + b"1 2 3\n4  5 6\n7 8 9\n",
            "nonl": base[:-1],
            "trailsp": HDR + b"1 2 3 \n4 5 6\n",
            "leadsp": HDR + b" 1 2 3\n4 5 6\n",
            "blankline": HDR + b"1 2 3\n\n4 5 6\n",
            "hdr_nonl": b"# only a header with no newline",
            "hashrow": HDR + b"1 2 3\n#x 5 6\n7 8 9\n",
        }
        for name, data in odd.items():
            st, _ = roundtrip(tmp, name, data)
            assert st["row_blocks"] >= 1, (name, st)
        # a bad row only spoils its own block
        many = HDR + b"".join(b"%d 2 3\n" % i for i in range(10)) + b"4 5\n"
        st, _ = roundtrip(tmp, "mixed", many, rows_per_block=4)
        assert st["blocks"] == 3 and st["row_blocks"] == 1 and st["rows"] == 11, st

        # concatenated chunks form one table, header once
        a = HDR + b"1 2 3\n4 5 6\n"
        b = HDR + b"7 8 9\n"
        _, pa = roundtrip(tmp, "ca", a)
        _, pb = roundtrip(tmp, "cb", b)
        cat = os.path.join(tmp, "cat.colxz")
        with open(cat, "wb") as f:
            for p in (pa, pb):
                f.write(open(p, "rb").read())
        got, n = decoded(cat)
        assert got == HDR + b"1 2 3\n4 5 6\n7 8 9\n" and n == 3, (got, n)

        # different headers: decode raises
        _, pc = roundtrip(tmp, "cc", b"# x y z\n1 2 3\n")
        bad = os.path.join(tmp, "bad.colxz")
        with open(bad, "wb") as f:
            for p in (pa, pc):
                f.write(open(p, "rb").read())
        try:
            decoded(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("different headers must raise")

        # verify() detects a mismatch
        src = os.path.join(tmp, "ca.dat")
        with open(src, "ab") as f:
            f.write(b"0 0 0\n")
        assert not ct.verify(src, pa)
    print("coltable_test: all passed")


if __name__ == "__main__":
    main()
