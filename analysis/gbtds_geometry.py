"""The GBTDS footprint and the simulator's scan region, in one place, for every Python script.

WHY THIS EXISTS. Until Deviation 69 five scripts each carried their own copy of the six field
circles (FIELDS, FoVRoman = 0.3003) and of Bulge.h's scan box. The simulator now takes the
footprint from the vendored adopted layout (Baseline/gbtds_layout/: field centres per roll, 18
detector rectangles per field) and builds its scan region by a distance rule, so the Python side
reads the same files and reproduces the same arithmetic here, once.

What it mirrors (Bulge.h / helper.cpp / Bulge_LSST.cpp):
  readGbtdsLayout, inDetector      -> sca_rects(), in_detector()
  SCAN_RUBIN_REACH + rField        -> scan_reach()
  the stratified sightline grid    -> scan_sightlines()   (checked against a run's own log by
                                                            the scripts that draw it)

Runs made before Deviation 69 used the notional layout: LEGACY_FIELDS / LEGACY_FOV_ROMAN, a
circle per field and no rolls. `run_geometry(run_dir)` tells the two apart from the run's
provenance (new runs record `# roman_layout`), so old tables are still read with the geometry
that produced them.
"""

import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAYOUT_DIR = os.path.join(ROOT, "Baseline", "gbtds_layout")
CENTERS_FILES = ("gbtds_spring_2026.4.3.centers", "gbtds_autumn_2026.4.3.centers")
SCA_FILES = ("sca_layout_spring.txt", "sca_layout_fall.txt")
LAYOUT_NAMES = ("spring", "autumn")

FOV_RUBIN = 1.75                     # Bulge.h FoV: Rubin's matching radius [deg]
SCAN_RUBIN_REACH = 2.0 * FOV_RUBIN   # Bulge.h SCAN_RUBIN_REACH
DD = 0.02                            # Bulge.h dd: the native grid unit [deg]

# The superseded notional layout (mtpenny/gbtds_optimizer layout_40395), for pre-Deviation-69 runs.
LEGACY_FIELDS = [(-0.417948, -1.2), (-0.008974, -1.2), (0.4, -1.2), (0.808974, -1.2),
                 (1.217948, -1.2), (0.0, -0.125)]
LEGACY_FOV_ROMAN = 0.3003


def read_centers(layout):
    """[(l, b)] x 6 for one roll (0 spring, 1 autumn): fields 1-5 then the Galactic Centre."""
    out = []
    with open(os.path.join(LAYOUT_DIR, CENTERS_FILES[layout])) as f:
        next(f)
        for line in f:
            p = line.split()
            if p:
                out.append((float(p[1]), float(p[2])))
    assert len(out) == 6, CENTERS_FILES[layout]
    return out


def placements():
    """The 12 (l, b, layout) placements RomanBaseline.dat visits, spring then autumn."""
    return [(l, b, k) for k in (0, 1) for l, b in read_centers(k)]


def sca_rects(layout):
    """(18, 4) array of detector rectangles (l0, l1, b0, b1), offsets from the field centre."""
    vert = {}
    with open(os.path.join(LAYOUT_DIR, SCA_FILES[layout])) as f:
        for line in f:
            p = line.split()
            if len(p) == 3:
                vert.setdefault(int(p[0]), []).append((float(p[1]), float(p[2])))
    r = np.array([[min(x for x, _ in v), max(x for x, _ in v),
                   min(y for _, y in v), max(y for _, y in v)] for _, v in sorted(vert.items())])
    assert r.shape == (18, 4), SCA_FILES[layout]
    return r


_RECTS = None


def _rects():
    global _RECTS
    if _RECTS is None:
        _RECTS = [sca_rects(0), sca_rects(1)]
    return _RECTS


def r_field():
    """Largest centre-to-detector-corner distance over both rolls (GbtdsLayout::rField)."""
    return max(np.hypot(r[:, [0, 1]][:, :, None], r[:, [2, 3]][:, None, :]).max() for r in _rects())


def sca_side():
    return min(min((r[:, 1] - r[:, 0]).min(), (r[:, 3] - r[:, 2]).min()) for r in _rects())


def outline_bbox(layout):
    r = _rects()[layout]
    return r[:, 0].min(), r[:, 1].max(), r[:, 2].min(), r[:, 3].max()


def scan_reach():
    return SCAN_RUBIN_REACH + r_field()


def in_detector(dl, db, layout):
    """Vectorised inDetector(): is the offset (dl, db) from a field centre on a detector?"""
    dl, db = np.asarray(dl, float), np.asarray(db, float)
    r = _rects()[layout]
    hit = np.zeros(np.broadcast(dl, db).shape, bool)
    for l0, l1, b0, b1 in r:
        hit |= (dl >= l0) & (dl <= l1) & (db >= b0) & (db <= b1)
    return hit


def on_detector(lon, lat):
    """(n, 2) bool: is each sky point on a detector in the spring / autumn roll (any field)?"""
    lon, lat = np.atleast_1d(np.asarray(lon, float)), np.atleast_1d(np.asarray(lat, float))
    out = np.zeros((lon.size, 2), bool)
    for l, b, k in placements():
        out[:, k] |= in_detector(lon - l, lat - b, k)
    return out


def scan_sightlines(stride, stride_roman=0):
    """Rebuild Bulge_LSST.cpp's stratified scan (full region, not --stub).

    Returns a dict of arrays: lon, lat, area (deg^2 represented), fine (in the footprint
    stratum), and the cells (lon, lat, representative index) for drawing. Same arithmetic and
    order as the C++ (grid points lonMin + i*fineStep, iLon-major, coarse-block representatives).
    """
    scaside = sca_side()
    grid = stride * DD
    if stride_roman == 0:
        stride_roman = stride
        if grid > scaside:
            stride_roman = max(k for k in range(stride, 0, -1)
                               if stride % k == 0 and k * DD <= scaside)
    ksub = stride // stride_roman
    fine = stride_roman * DD
    pl = placements()
    reach = scan_reach()
    fl = [p[0] for p in pl]
    fb = [p[1] for p in pl]
    lon_min = grid * np.floor((min(fl) - reach) / grid)
    lon_max = max(fl) + reach
    lat_min = grid * np.floor((min(fb) - reach) / grid)
    lat_max = max(fb) + reach
    n_lon = int(np.floor((lon_max - lon_min) / grid + 1e-9)) + 1
    n_lat = int(np.floor((lat_max - lat_min) / grid + 1e-9)) + 1
    bbox = [outline_bbox(k) for k in (0, 1)]

    def in_scan(l, b):
        return any(np.hypot(l - x, b - y) <= reach for x, y, _ in pl)

    def in_fine(l, b):
        for x, y, k in pl:
            dl0, dl1, db0, db1 = bbox[k]
            if dl0 - fine <= l - x <= dl1 + fine and db0 - fine <= b - y <= db1 + fine:
                return True
        return False

    nlf, nbf = n_lon * ksub, n_lat * ksub
    cnt, rep = {}, {}
    for i in range(nlf):
        lon = lon_min + i * fine
        for j in range(nbf):
            lat = lat_min + j * fine
            if not in_scan(lon, lat) or in_fine(lon, lat):
                continue
            blk = (i // ksub, j // ksub)
            rep.setdefault(blk, (i, j))
            cnt[blk] = cnt.get(blk, 0) + 1
    cell = fine * fine
    L, B, A, F, cells = [], [], [], [], []
    for i in range(nlf):
        lon = lon_min + i * fine
        for j in range(nbf):
            lat = lat_min + j * fine
            if not in_scan(lon, lat):
                continue
            if in_fine(lon, lat):
                cells.append((lon, lat, len(L)))
                L.append(lon); B.append(lat); A.append(cell); F.append(True)
            else:
                blk = (i // ksub, j // ksub)
                ri, rj = rep[blk]
                if (ri, rj) == (i, j):
                    L.append(lon); B.append(lat); A.append(cnt[blk] * cell); F.append(False)
                cells.append((lon, lat, None))
    L, B = np.array(L), np.array(B)
    # Resolve each coarse cell's representative index.
    index = {(round(l, 6), round(b, 6)): n for n, (l, b) in enumerate(zip(L, B))}
    out_cells = []
    for lon, lat, n in cells:
        if n is None:
            i = int(round((lon - lon_min) / fine)); j = int(round((lat - lat_min) / fine))
            ri, rj = rep[(i // ksub, j // ksub)]
            n = index[(round(lon_min + ri * fine, 6), round(lat_min + rj * fine, 6))]
        out_cells.append((lon, lat, n))
    return dict(lon=L, lat=B, area=np.array(A), fine=np.array(F), cells=out_cells,
                stride_roman=stride_roman, fine_step=fine,
                bounds=(lon_min, lon_max, lat_min, lat_max))


def run_geometry(run_dir):
    """'gbtds2026' for runs on the adopted layout, 'legacy' for earlier ones (provenance)."""
    for name in ("run_provenance.txt", os.path.join("files", "MONTLMC", "files", "run_provenance.txt")):
        p = os.path.join(run_dir, name)
        if os.path.exists(p):
            return "gbtds2026" if "# roman_layout" in open(p).read() else "legacy"
    raise FileNotFoundError(f"no run_provenance.txt under {run_dir}")


def legacy_in_footprint(lon, lat):
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    return np.any([np.hypot(lon - l, lat - b) <= LEGACY_FOV_ROMAN for l, b in LEGACY_FIELDS], axis=0)


def read_roman_visits(path=os.path.join(ROOT, "Baseline", "RomanBaseline.dat")):
    """A Roman visit list as a DataFrame (l, b, time, field, layout).

    Recognises both formats: the 9-column list of the adopted layout, and the 7-column list of
    runs before Deviation 69 (Baseline/legacy_layout40395/), for which layout = -1 marks
    "circle of LEGACY_FOV_ROMAN about (l, b)".
    """
    import pandas as pd
    with open(path) as f:
        ncol = len(f.readline().lstrip("#").split())
    if ncol == 9:
        v = pd.read_csv(path, sep=r"\s+", comment="#", header=None, usecols=[3, 4, 5, 7, 8],
                        names=["ID", "RA", "Dec", "l", "b", "time", "sig5", "field", "layout"])
    elif ncol == 7:
        v = pd.read_csv(path, sep=r"\s+", comment="#", header=None, usecols=[3, 4, 5],
                        names=["ID", "RA", "Dec", "l", "b", "time", "sig5"])
        v["layout"] = -1
        centres = v[["l", "b"]].drop_duplicates().reset_index(drop=True)
        v["field"] = v.merge(centres.reset_index(), on=["l", "b"], how="left")["index"].to_numpy()
    else:
        raise ValueError(f"{path}: {ncol} header columns, expected 7 (legacy) or 9")
    return v


def visit_covers(visits, l0, b0):
    """Boolean per visit: does it image the sky point (l0, b0)? Same test as the simulator."""
    lay = visits["layout"].to_numpy()
    dl = l0 - visits["l"].to_numpy()
    db = b0 - visits["b"].to_numpy()
    out = np.zeros(len(visits), bool)
    leg = lay < 0
    out[leg] = np.hypot(dl[leg], db[leg]) <= LEGACY_FOV_ROMAN
    for k in (0, 1):
        m = lay == k
        out[m] = in_detector(dl[m], db[m], k)
    return out


def field_label(lon, lat, visits):
    """Per sky point, the GBTDS field(s) that image it: 'F<i>' (same field in every roll that
    sees it), 'F<i>/F<j>' (spring / autumn differ), or 'outside'. Legacy lists: the nearest
    centre within LEGACY_FOV_ROMAN, as the pre-Deviation-69 F1 did."""
    lon, lat = np.atleast_1d(np.asarray(lon, float)), np.atleast_1d(np.asarray(lat, float))
    pl = visits[["l", "b", "field", "layout"]].drop_duplicates().to_numpy()
    if (pl[:, 3] < 0).all():
        best = np.full(lon.size, np.inf)
        lab = np.array(["outside"] * lon.size, dtype=object)
        for l, b, f, _ in pl:
            d = np.hypot(lon - l, lat - b)
            hit = (d < LEGACY_FOV_ROMAN) & (d < best)
            best[hit] = d[hit]
            lab[hit] = f"F{int(f)}({l:+.3f},{b:+.3f})"
        return lab
    per = {0: np.full(lon.size, -1), 1: np.full(lon.size, -1)}
    for l, b, f, k in pl:
        per[int(k)][in_detector(lon - l, lat - b, int(k))] = int(f)
    lab = []
    for s, a in zip(per[0], per[1]):
        ids = sorted({x for x in (s, a) if x >= 0})
        lab.append("outside" if not ids else f"F{ids[0]}" if len(ids) == 1 else f"F{s}/F{a}")
    return np.array(lab, dtype=object)


if __name__ == "__main__":
    print(f"r_field {r_field():.5f} deg, detector side {sca_side():.5f} deg, "
          f"scan reach {scan_reach():.4f} deg")
    for k in (0, 1):
        r = sca_rects(k)
        print(f"{LAYOUT_NAMES[k]}: area/field {((r[:, 1] - r[:, 0]) * (r[:, 3] - r[:, 2])).sum():.5f}"
              f" deg^2, outline {outline_bbox(k)}")
