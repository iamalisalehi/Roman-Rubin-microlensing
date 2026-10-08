#!/usr/bin/env python3
"""Concatenate the chunks of one population's scan into the single run directory the analysis reads.

    python3 pipeline/merge_chunks.py RUN_ROOT/<pop>        # reads <pop>/plan.txt, <pop>/chunks/NNNN/

Why this is exact: every sightline re-seeds the random-number generator from
(seed, its scan index), and every output file gets rows appended once per sightline (or per event).
So running [0,a), [a,b), ... in separate directories and concatenating the files in index order gives
the bytes a single unchunked run writes. A later chunk adds only its own copy of the leading '#' header
lines of each file; those are dropped here if they equal the first chunk's.

What is written into <pop>/ (each file to a .tmp name first, renamed when complete):
  test<tag>.dat, h3_pair.dat, ...         every regular file in a chunk's top directory
  files/MONTLMC/files/*                   every output file of the run, except run_provenance.txt
  run.log                                 the chunks' logs, in chunk order
  files/MONTLMC/files/run_provenance.txt  chunk 0's, with start_index/end_index set to the merged range,
                                          `merged_chunks`/`chunk_size` added, and the end-of-run
                                          sightline-outcome counts and areas summed over the chunks
  MERGED                                  written last: what was merged (read by `pipeline.sh status`)

NOT exact: EfLMC<tag>.dat / EfLMC<tag>B.dat (see CUMULATIVE below). They are concatenated so the file has its
usual shape (one block per sightline), but each block is cumulative over its own chunk only. No script in
analysis/ reads them.

Refuses (exit 1, listing them) if a chunk has no DONE file or the chunks do not tile [0, N).
Stdlib only.
"""
import argparse
import os
import re
import sys

BLOCK = 8 << 20
SKIP_TOP = {"DONE", "RUNNING", "run.log"}          # chunk bookkeeping, not outputs
PROV = "files/MONTLMC/files/run_provenance.txt"
# EfLMC<tag>.dat and EfLMC<tag>B.dat hold, per sightline, the detection efficiency accumulated over all
# sightlines so far. The running counts restart in each chunk and only percentages are written, so the
# chunks cannot be added up and the concatenation differs from the unchunked file.
CUMULATIVE = re.compile(r"(^|/)EfLMC.*\.dat$")


def read_kv(path):
    """key=value lines (plan.txt, DONE) -> dict of strings."""
    d = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, v = line.split("=", 1)
                d[k.strip()] = v.strip().strip('"')
    return d


def output_files(chunk_dir):
    """Relative paths of every output file a chunk wrote (the ones to concatenate)."""
    out = []
    for name in sorted(os.listdir(chunk_dir)):
        p = os.path.join(chunk_dir, name)
        if os.path.isfile(p) and not os.path.islink(p) and name not in SKIP_TOP:
            out.append(name)
    mdir = os.path.join(chunk_dir, "files/MONTLMC/files")
    if os.path.isdir(mdir):
        for name in sorted(os.listdir(mdir)):
            if os.path.isfile(os.path.join(mdir, name)) and name != "run_provenance.txt":
                out.append("files/MONTLMC/files/" + name)
    return out


def leading_header(path):
    """(header bytes, offset of the first non-header byte): the file's leading '#' lines."""
    hdr = b""
    with open(path, "rb") as fh:
        while True:
            line = fh.readline()
            if line.startswith(b"#"):
                hdr += line
            else:
                return hdr, len(hdr)


def merge_file(rel, chunk_dirs, dest_dir):
    """Concatenate one output file over the chunks. Returns (rows, bytes, n_chunks_with_file)."""
    dest = os.path.join(dest_dir, rel)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".tmp"
    header = None            # first non-empty chunk's header; later chunks must repeat it
    rows = nbytes = present = 0
    src_bytes = dropped = 0
    with open(tmp, "wb") as out:
        for cd in chunk_dirs:
            p = os.path.join(cd, rel)
            if not os.path.isfile(p) or os.path.getsize(p) == 0:
                continue
            present += 1
            src_bytes += os.path.getsize(p)
            hdr, off = leading_header(p)
            if header is None:
                header, off = hdr, 0                     # keep the first header
            elif hdr != header:
                os.remove(tmp)
                sys.exit(f"{p}: its leading '#' lines differ from the first chunk's; "
                         "chunks of one run must share a header. Not merging.")
            else:
                dropped += off
            last = b"\n"
            with open(p, "rb") as fh:
                fh.seek(off)
                while True:
                    buf = fh.read(BLOCK)
                    if not buf:
                        break
                    out.write(buf)
                    rows += buf.count(b"\n")
                    last = buf[-1:]
            if last != b"\n":
                os.remove(tmp)
                sys.exit(f"{p}: does not end in a newline (a chunk was cut mid-row?). Not merging.")
    nbytes = os.path.getsize(tmp)
    if nbytes != src_bytes - dropped:
        os.remove(tmp)
        sys.exit(f"{rel}: merged size {nbytes} != sum of chunk sizes {src_bytes} - dropped headers {dropped}")
    os.replace(tmp, dest)
    hrows = header.count(b"\n") if header else 0
    return rows, nbytes, present, hrows


def merge_log(chunk_dirs, dest_dir):
    dest = os.path.join(dest_dir, "run.log")
    with open(dest + ".tmp", "wb") as out:
        for cd in chunk_dirs:
            p = os.path.join(cd, "run.log")
            if os.path.isfile(p):
                with open(p, "rb") as fh:
                    while True:
                        buf = fh.read(BLOCK)
                        if not buf:
                            break
                        out.write(buf)
    os.replace(dest + ".tmp", dest)


def merge_provenance(chunk_dirs, dest_dir, n_total, chunk_size):
    """Chunk 0's run_provenance.txt, made to describe the merged run."""
    lines = open(os.path.join(chunk_dirs[0], PROV)).read().splitlines()
    marker = "# ---- sightline outcome"
    # sum the end-of-run block over all chunks (counts and areas)
    sums = {}
    for cd in chunk_dirs:
        txt = open(os.path.join(cd, PROV)).read().splitlines()
        if not any(l.startswith(marker) for l in txt):
            sys.exit(f"{cd}/{PROV}: no sightline-outcome block -- that chunk did not finish")
        seen = False
        for l in txt:
            if l.startswith(marker):
                seen = True
                continue
            m = re.match(r"^#\s+(sightlines_\w+|area_\w+_deg2)\s+(\S+)", l)
            if seen and m:
                sums[m.group(1)] = sums.get(m.group(1), 0.0) + float(m.group(2))
    out, in_block = [], False
    for l in lines:
        if l.startswith(marker):
            out += [f"# merged_chunks       {len(chunk_dirs)}   # {len(chunk_dirs)} chunks of "
                    f"{chunk_size} sightlines concatenated by pipeline/merge_chunks.py; "
                    "the block below is summed over them (areas are sums of rounded values)",
                    l]
            in_block = True
            continue
        if in_block:
            m = re.match(r"^#\s+(sightlines_\w+|area_\w+_deg2)\s+", l)
            if m and m.group(1) in sums:
                v = sums[m.group(1)]
                val = f"{int(round(v))}" if m.group(1).startswith("sightlines_") else f"{v:.6g}"
                out.append(f"# {m.group(1):<23} {val}")
                continue
            out.append(l)
            continue
        m = re.match(r"^#\s+start_index\s+", l)
        if m:
            out.append("# start_index         0   # merged run: the whole scan, assembled from chunks")
            continue
        m = re.match(r"^#\s+end_index\s+", l)
        if m:
            out.append(f"# end_index           {n_total}   # merged range [0, {n_total}); "
                       "see merged_chunks below")
            continue
        out.append(l)
    p = os.path.join(dest_dir, PROV)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p + ".tmp", "w") as fh:
        fh.write("\n".join(out) + "\n")
    os.replace(p + ".tmp", p)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pop_dir", help="RUN_ROOT/<pop> (holds plan.txt and chunks/)")
    a = ap.parse_args()
    pop = os.path.abspath(a.pop_dir)
    plan = read_kv(os.path.join(pop, "plan.txt"))
    n_total, csize, nchunks = int(plan["n_sightlines"]), int(plan["chunk_size"]), int(plan["n_chunks"])
    cdirs = [os.path.join(pop, "chunks", f"{k:04d}") for k in range(nchunks)]

    # 1. every chunk finished, and the ranges tile [0, N)
    bad, edges = [], []
    for k, cd in enumerate(cdirs):
        done = os.path.join(cd, "DONE")
        if not os.path.isfile(done):
            bad.append(f"chunk {k}: no DONE file ({cd})")
            continue
        d = read_kv(done)
        s, e = int(d["start"]), int(d["end"])
        if (s, e) != (k * csize, min((k + 1) * csize, n_total)):
            bad.append(f"chunk {k}: covers [{s},{e}), expected "
                       f"[{k * csize},{min((k + 1) * csize, n_total)})")
        edges.append((s, e))
    if bad:
        print(f"merge: {os.path.basename(pop)}: cannot merge, {len(bad)} problem(s):", file=sys.stderr)
        for b in bad:
            print("   " + b, file=sys.stderr)
        sys.exit(1)
    pos = 0
    for s, e in edges:
        if s != pos:
            sys.exit(f"chunks do not tile: gap/overlap at sightline {pos} (next chunk starts at {s})")
        pos = e
    if pos != n_total:
        sys.exit(f"chunks end at sightline {pos}, the scan has {n_total}")

    # 2. the files: the union over chunks (a chunk with no events may lack some), concatenated in order
    names = []
    for cd in cdirs:
        for rel in output_files(cd):
            if rel not in names:
                names.append(rel)
    print(f"merge: {os.path.basename(pop)}: {nchunks} chunks tile [0,{n_total}); "
          f"{len(names)} output files")
    summary = []
    for rel in names:
        rows, nbytes, present, hrows = merge_file(rel, cdirs, pop)
        summary.append((rel, rows, nbytes, present, hrows))
        note = ""
        if CUMULATIVE.search(rel):
            note = ("  NOTE: running efficiency that restarts in every chunk -- NOT equal to an unchunked "
                    "file; do not read an efficiency curve from it")
        print(f"   {rel:<42} {rows:>12,} lines {nbytes / 1e6:>10.1f} MB   (in {present}/{nchunks} chunks"
              f"{', header %d lines kept once' % hrows if hrows else ''}){note}")
    merge_log(cdirs, pop)
    merge_provenance(cdirs, pop, n_total, csize)
    with open(os.path.join(pop, "MERGED"), "w") as fh:
        fh.write(f"chunks={nchunks}\nn_sightlines={n_total}\nchunk_size={csize}\n")
        for rel, rows, nbytes, present, hrows in summary:
            fh.write(f"file {rel} lines={rows} bytes={nbytes}\n")
    print("merge: done")


if __name__ == "__main__":
    main()
