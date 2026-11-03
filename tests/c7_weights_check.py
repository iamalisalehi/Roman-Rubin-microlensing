"""Check the two event-counting conventions in analysis/romanlib.py (every star / per object).

    .roman/bin/python tests/c7_weights_check.py

Checks: on the legacy run samples_bulge_S3c, attach_weight(convention="per-object") equals bc5ca0e's
attach_weight bit for bit; on stub runs of the new code (every-star run, and a --blend-thinning run)
each convention is the documented product of event_weight and pKeep; pKeep agrees with the legacy
rebuild acceptance_probability (1e-5 relative, same zero pattern); a missing provenance file and a
map without Rubin depths raise errors that say what to pass instead. Exits non-zero on failure.
"""
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "analysis"))
import romanlib as R          # noqa: E402

OLD_RUN = "/home/ali/Documents/PhD/Offline_project/Roman/runs/samples_bulge_S3c"
SCRATCH = "/home/ali/.claude/jobs/d695a888/tmp/c7"
NEW_RUN, LEG_RUN = f"{SCRATCH}/t3_new_bulge", f"{SCRATCH}/t3_leg_bulge"
fails = []


def check(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails.append(msg)


def paths(run):
    return (os.path.join(run, "test5.dat"),
            os.path.join(run, "files/MONTLMC/files/MapLMC5.dat"))


def load(run):
    table, mp = paths(run)
    return R.load_events(table, keep=R.keep_weightable(mp)), mp


def raises(fn, exc, *words):
    try:
        fn()
    except exc as e:
        miss = [w for w in words if w not in str(e)]
        return not miss, f"message lacks {miss}: {e}" if miss else ""
    except Exception as e:                       # noqa: BLE001
        return False, f"raised {type(e).__name__}: {e}"
    return False, "did not raise"


def old_module(tmp):
    src = subprocess.run(["git", "-C", ROOT, "show", "bc5ca0e:analysis/romanlib.py"],
                         check=True, capture_output=True, text=True).stdout
    path = os.path.join(tmp, "romanlib_bc5ca0e.py")
    with open(path, "w") as fh:
        fh.write(src)
    spec = importlib.util.spec_from_file_location("romanlib_bc5ca0e", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    tmp = tempfile.mkdtemp(prefix="c7_weights_")
    try:
        # 1. legacy run: the old numbers, bit for bit
        old = old_module(tmp)
        df, mp = load(OLD_RUN)
        w_old, _ = old.attach_weight(df, mp)
        w_new, label = R.attach_weight(df, mp, convention="per-object")
        check(np.array_equal(w_old.to_numpy(), w_new.to_numpy()),
              f"legacy run, per-object == bc5ca0e attach_weight bit for bit ({len(df)} rows)")
        check("per-object" in label, f"label names the convention: {label!r}")
        check(R.selection_mode(os.path.join(os.path.dirname(mp), "run_provenance.txt")) == "per-object",
              "run with no source_selection line is a per-object run")

        # 2. new-code stubs
        for run, mode in ((NEW_RUN, "every-star"), (LEG_RUN, "per-object")):
            df, mp = load(run)
            sl = R.load_sightlines(mp)
            w = R.event_weight(df, sl).to_numpy()
            pk = df["pKeep"].to_numpy()
            check(R.selection_mode(os.path.join(os.path.dirname(mp), "run_provenance.txt")) == mode,
                  f"{os.path.basename(run)}: selection_mode == {mode}")
            wa, label_a = R.attach_weight(df, mp, convention="every-star")
            if mode == "every-star":
                wd, label_d = R.attach_weight(df, mp)
                check(np.array_equal(wd.to_numpy(), wa.to_numpy()) and "every-star" in label_d,
                      f"default convention is every-star: {label_d!r}")
            wo, _ = R.attach_weight(df, mp, convention="per-object")
            if mode == "every-star":
                check(np.array_equal(wa.to_numpy(), w), f"{mode} run: every-star == event_weight")
                check(np.array_equal(wo.to_numpy(), w * pk), f"{mode} run: per-object == event_weight * pKeep")
            else:
                check(np.array_equal(wo.to_numpy(), w), f"{mode} run: per-object == event_weight")
                pos = (w > 0) & (pk > 0)
                det = (df["detJ"] == 1).to_numpy()
                check(bool((pk[det] > 0).all()) and
                      np.array_equal(wa.to_numpy()[pos], w[pos] / pk[pos]) and
                      bool((wa.to_numpy()[~pos] == 0).all()) and bool((~pos & (w > 0)).any()),
                      f"{mode} run: every-star == event_weight / pKeep; 0 for never-kept draws "
                      f"(no detections among them: {int((~pos & (w > 0)).sum())} rows)")
            # pKeep against the legacy rebuild
            leg = R.acceptance_probability(df, sl)
            check(np.array_equal(pk > 0, leg > 0),
                  f"{os.path.basename(run)}: pKeep and acceptance_probability have the same zero pattern")
            nz = pk > 0
            rel = np.abs(leg[nz] - pk[nz]) / pk[nz]
            check(bool(rel.max() < 1e-5), f"{os.path.basename(run)}: pKeep == acceptance_probability, "
                                          f"max rel diff {rel.max():.2e} < 1e-5")
            check(np.array_equal(R.keep_probability(df, sl, mode), pk),
                  f"{os.path.basename(run)}: keep_probability returns the pKeep column")

        # 3. keep_probability: legacy rebuild only for per-object runs
        df, mp = load(NEW_RUN)
        sl = R.load_sightlines(mp)
        ok, why = raises(lambda: R.keep_probability(df.drop(columns="pKeep"), sl, "every-star"),
                         ValueError, "pKeep")
        check(ok, f"every-star table without pKeep -> ValueError {why}")
        got = R.keep_probability(df.drop(columns="pKeep"), sl, "per-object")
        check(np.array_equal(got, R.acceptance_probability(df, sl)),
              "per-object table without pKeep -> acceptance_probability")
        ok, why = raises(lambda: R.count_weight(np.array([1.0, 2.0]), np.array([0.5, 0.0]),
                                                "per-object", "every-star"), ValueError, "pKeep")
        check(ok, f"count_weight: w > 0 with p_keep <= 0 -> ValueError {why}")
        got = R.count_weight(np.array([1.0, 2.0]), np.array([0.5, 0.0]), "per-object", "every-star",
                             detected=np.array([True, False]))
        check(got.tolist() == [2.0, 0.0], "count_weight: undetected never-kept draw -> 0")
        ok, why = raises(lambda: R.count_weight(np.array([1.0, 2.0]), np.array([0.5, 0.0]),
                                                "per-object", "every-star",
                                                detected=np.array([True, True])), ValueError, "pKeep")
        check(ok, f"count_weight: detected never-kept draw -> ValueError {why}")
        ok, why = raises(lambda: R.count_weight(np.ones(2), np.ones(2), "both", "every-star"), ValueError)
        check(ok, f"count_weight: unknown run mode -> ValueError {why}")
        ok, why = raises(lambda: R.count_weight(np.ones(2), np.ones(2), "every-star", "neither"), ValueError)
        check(ok, f"count_weight: unknown convention -> ValueError {why}")
        pr = os.path.join(tmp, "prov.txt")
        with open(pr, "w") as fh:
            fh.write("# source_selection    sideways   # x\n")
        ok, why = raises(lambda: R.selection_mode(pr), ValueError, "sideways")
        check(ok, f"selection_mode: unknown value -> ValueError {why}")

        # 4. no provenance file: raise, say what to pass
        d = os.path.join(tmp, "noprov")
        os.makedirs(d)
        shutil.copy(paths(OLD_RUN)[1], d)
        df, _ = load(OLD_RUN)
        mp_np = os.path.join(d, "MapLMC5.dat")
        ok, why = raises(lambda: R.attach_weight(df, mp_np, convention="per-object"),
                         FileNotFoundError, "run_provenance.txt", "mode", "--convention")
        check(ok, f"map without provenance -> FileNotFoundError naming file, mode, --convention {why}")
        w_exp, _ = R.attach_weight(df, mp_np, convention="per-object", mode="per-object")
        check(np.array_equal(w_exp.to_numpy(), w_old.to_numpy()), "explicit mode= skips the provenance file")

        # 5. legacy map without Rubin depths (samples_bulge_S3c's map has none): every-star cannot
        #    be rebuilt, per-object still works
        mp_nd = paths(OLD_RUN)[1]
        check(bool(np.isnan(R.load_sightlines(mp_nd)["depth5_r"]).all()), "the legacy map has no depth5_*")
        ok, why = raises(lambda: R.attach_weight(df, mp_nd, convention="every-star"),
                         ValueError, "--convention per-object")
        check(ok, f"legacy map without depth5_* asked for every-star -> ValueError {why}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\nFAIL: " + "; ".join(fails) if fails else "\nPASS: all checks")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
