"""romanlib reads .colxz exactly as it reads the plain table (dtypes included).

Standalone script, not pytest:  .roman/bin/python tests/coltable_romanlib_test.py
Scratch files go to the system temp dir (set TMPDIR to redirect).
"""
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

from pandas.testing import assert_frame_equal

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.join(HERE, "..", "analysis")
sys.path.insert(0, ANALYSIS)
import coltable as ct  # noqa: E402
import romanlib as R  # noqa: E402

SAMPLE = os.path.expanduser("~/.claude/jobs/b0aec636/tmp/kin_sample/rv/bulge/389/test5.dat")

CHILD = ("import sys; sys.path.insert(0, sys.argv[1]); import romanlib as R; "
         "uc = sys.argv[3].split(',') if len(sys.argv) > 3 else None; "
         "R.load_events(sys.argv[2], usecols=uc)")


def main():
    tmp = tempfile.mkdtemp(prefix="cx_rl_")
    try:
        plain = os.path.join(tmp, "t.dat")
        shutil.copy(SAMPLE, plain)
        with open(plain, "rb") as f, gzip.open(plain + ".gz", "wb") as g:
            shutil.copyfileobj(f, g)
        ct.encode(plain, plain + ".colxz")

        cols = open(plain).readline().lstrip("#").split()
        sub = [cols[i] for i in (0, 3, 7, 20, len(cols) - 1)]
        ref_all = R.load_events(plain)
        ref_sub = R.load_events(plain, usecols=sub)

        # resolver: bare name, order plain, .colxz, .gz
        assert R.resolve_table(plain) == plain
        os.rename(plain, plain + ".hide")
        assert R.resolve_table(plain) == plain + ".colxz", R.resolve_table(plain)
        os.rename(plain + ".colxz", plain + ".colxz.hide")
        assert R.resolve_table(plain) == plain + ".gz"
        os.rename(plain + ".colxz.hide", plain + ".colxz")
        os.rename(plain + ".hide", plain)

        for name in (plain + ".gz", plain + ".colxz", plain):
            assert_frame_equal(R.load_events(name), ref_all)
            assert_frame_equal(R.load_events(name, usecols=sub), ref_sub)
            assert_frame_equal(R.load_events(name, usecols=sub[::-1]), ref_sub)
        # via the bare name with only the .colxz (and .gz) present
        os.rename(plain, plain + ".hide")
        assert_frame_equal(R.load_events(plain), ref_all)
        os.rename(plain + ".hide", plain)

        # keep / chunksize / narrow
        kw = dict(keep=lambda c: c["tE"] > 20, chunksize=7000)
        assert_frame_equal(R.load_events(plain + ".colxz", **kw), R.load_events(plain, **kw))
        assert_frame_equal(R.load_events(plain + ".colxz", narrow=True, usecols=sub),
                           R.load_events(plain, narrow=True, usecols=sub))
        try:
            R.load_events(plain + ".colxz", usecols=["nope"])
            raise SystemExit("missing column accepted")
        except ValueError:
            pass

        # many blocks
        small = os.path.join(tmp, "small.colxz")
        ct.encode(plain, small, rows_per_block=1000)
        assert sum(1 for _ in ct.iter_blocks(small)) > 3
        assert_frame_equal(R.load_events(small), ref_all)
        assert_frame_equal(R.load_events(small, usecols=sub), ref_sub)

        # a row-block fallback among column blocks (a double space forces it)
        mixed = os.path.join(tmp, "mixed.dat")
        lines = open(plain).read().split("\n")
        lines[5] = lines[5].replace(" ", "  ", 1)
        open(mixed, "w").write("\n".join(lines))
        ct.encode(mixed, mixed + ".colxz", rows_per_block=1000)
        kinds = {b.raw_rows is not None for b in ct.iter_blocks(mixed + ".colxz")}
        assert kinds == {True, False}, kinds
        assert_frame_equal(R.load_events(mixed + ".colxz"), R.load_events(mixed))
        assert_frame_equal(R.load_events(mixed + ".colxz", usecols=sub), R.load_events(mixed, usecols=sub))

        # int in some blocks, float in others; nan tokens
        odd = os.path.join(tmp, "odd.dat")
        with open(odd, "w") as f:
            f.write("# a b c\n")
            for i in range(30):
                f.write("%d %s %s\n" % (i, ("%d" % i) if i < 20 else "%d.5" % i,
                                        "nan" if i % 7 == 0 else "-1.5e-3"))
        ct.encode(odd, odd + ".colxz", rows_per_block=10)
        assert_frame_equal(R.load_events(odd + ".colxz"), R.load_events(odd))

        # header-only table
        hdr = os.path.join(tmp, "hdr.dat")
        open(hdr, "w").write("# a b c\n")
        ct.encode(hdr, hdr + ".colxz")
        e = R.load_events(hdr + ".colxz")
        assert list(e.columns) == ["a", "b", "c"] and len(e) == 0

        # memory: the 5-column .colxz read is lighter than the all-column read
        sub_only = peak_rss_mb_isolated(plain + ".colxz", sub)
        all_only = peak_rss_mb_isolated(plain + ".colxz", None)
        print("peak RSS isolated: all %.0f MB, 5 columns %.0f MB" % (all_only, sub_only))
        assert sub_only < all_only
        print("coltable_romanlib_test: OK")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def peak_rss_mb_isolated(path, usecols):
    """Each measurement in its own grandchild-free process: a small parent that reports the child's maxrss."""
    code = ("import resource, subprocess, sys; "
            "subprocess.run([sys.executable, '-c', sys.argv[1]] + sys.argv[2:], check=True); "
            "print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0)")
    args = [ANALYSIS, path] + ([",".join(usecols)] if usecols else [])
    out = subprocess.run([sys.executable, "-c", code, CHILD] + args, check=True,
                         capture_output=True, text=True).stdout
    return float(out.strip().splitlines()[-1])


if __name__ == "__main__":
    main()
