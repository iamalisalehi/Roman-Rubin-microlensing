"""Regenerate the overview reports' two footprint figures. Run from the repo root.

    .roman/bin/python Report/overview/make_footprints.py            # adopted layout (Dev. 69)
        -> figures/footprint_20261001/: the scan and Roman's detectors as the simulator now builds
           them (analysis/gbtds_geometry), and their overlap with the real tiles; checked against a
           run's log with --log <run.log> once a run exists.
    .roman/bin/python Report/overview/make_footprints.py --legacy   # the old report's figures
        -> figures/footprint_20260929/, as described below (pre-Deviation-69 runs).

LEGACY MODE (the original description):

  figures/footprint_20260929/footprint_simulated.{pdf,png}
      The sky the production runs scanned: every 0.1-deg cell of the tiling, coloured by the
      number of Rubin visits at the sightline that represents it, with Roman's footprint (the
      fine-grid sightlines) and the six modelled GBTDS field circles drawn over it.
  figures/footprint_20260929/footprint_vs_gbtds.{pdf,png}
      Roman's simulated footprint drawn over an Aladin view of the real GBTDS tiles
      (Whitepaper/roman_967_both_aladinX.png: spring and autumn tile positions together).

THE TILING IS REBUILT, NOT ASSUMED. The grid below mirrors Bulge_LSST.cpp's stratified scan
(l1/l2/b1/b2/wid/lx/bx from Bulge.h; fine 0.1 deg, coarse 0.2 deg, kSub = 2): every fine cell that
survives the corner cut is either a footprint sightline itself or is represented by its coarse
block's first non-footprint cell. (The code's area bookkeeping has a grid point stand for the cell
extending from it in +l and +b; for DRAWING, each fine cell is centred on its own grid point, so the
figure shows where the sightlines actually sample the sky rather than a half-cell offset.)
The script refuses to draw unless the rebuilt sightlines are exactly the 1,829 the production run
entered (runs/prod_bulge_20260924/run.log), with the same Rubin/Roman coverage.

THE SCREENSHOT HAS NO WCS, so it is calibrated from the Galactic grid Aladin drew on it: vertical
lines at l = 1.5, 1.0, 0.5, 0.0, -0.5, -1.0 deg and horizontal ones at b = 0, -0.5, ..., -2.0 are
found as the image columns/rows that are green over most of their length, and a linear fit gives
pixels -> (l, b). The grid is rectilinear over this 2.8 x 2.4 deg view, so the fit is exact to a
pixel (1/348 deg). The real tiles are then the filled green detector outlines.
"""
import os
import re
import sys

import numpy as np
from PIL import Image
from scipy import ndimage

sys.path.insert(0, "analysis")
import plotstyle as ps                                    # noqa: E402
import matplotlib.pyplot as plt                           # noqa: E402
from matplotlib.patches import Circle, Rectangle          # noqa: E402
from matplotlib.collections import PatchCollection        # noqa: E402
from matplotlib.colors import LogNorm                     # noqa: E402

OUT = "figures/footprint_20260929"
LOG = "runs/prod_bulge_20260924/run.log"
IMG = "Whitepaper/roman_967_both_aladinX.png"

# ---- mirrored from Bulge.h / Bulge_LSST.cpp / Baseline/generateRomanBaseline.py ----
WID = 3.5 / 2
L1, L2 = -0.219 - 0.2 - WID, 1.4134 + 0.2 + WID
B1, B2 = -1.64 - 0.2 - WID, -0.85 + 0.2 + WID
LX, BX = 1.0053 - 0.2 - WID, -1.64 + 0.2 + WID           # corner cut: l < LX and b > BX
LON_MIN, LON_MAX, LAT_MIN, LAT_MAX = L1 - WID, L2 + WID, B1 - WID, B2 + WID
GRID, FINE, KSUB = 0.2, 0.1, 2
HALF = 0.5 * 0.1          # draw cells centred on their grid point
FOV_ROMAN = 0.3003
FIELDS = [(-0.417948, -1.2), (-0.008974, -1.2), (0.4, -1.2), (0.808974, -1.2),
          (1.217948, -1.2), (0.0, -0.125)]


def tiling():
    """Every fine cell kept, and the sightline (lon, lat) that represents it."""
    n_lon = int(np.floor((LON_MAX - LON_MIN) / GRID + 1e-9)) + 1
    n_lat = int(np.floor((LAT_MAX - LAT_MIN) / GRID + 1e-9)) + 1
    inside = lambda l, b: any(np.hypot(l - fl, b - fb) <= FOV_ROMAN for fl, fb in FIELDS)
    cells, rep = [], {}
    for i in range(n_lon * KSUB):
        lon = LON_MIN + i * FINE
        for j in range(n_lat * KSUB):
            lat = LAT_MIN + j * FINE
            if lon < LX and lat > BX:
                continue
            blk = (i // KSUB, j // KSUB)
            if inside(lon, lat):
                cells.append((lon, lat, (lon, lat), True))
            else:
                rep.setdefault(blk, (lon, lat))
                cells.append((lon, lat, rep[blk], False))
    return cells


def logged_sightlines():
    """{(lon, lat): (Rubin epochs, Roman epochs)} for every sightline the run entered."""
    import subprocess
    lines = subprocess.run(["grep", "-E", r"^longtitude:|^ndd \(", LOG], capture_output=True,
                           text=True, errors="replace", check=True).stdout.splitlines()
    out, key, nl = {}, None, None
    if True:
        for line in lines:
            m = re.match(r"^longtitude:\s*(\S+)\s+latitude:\s*(\S+)", line)
            if m:
                key = (round(float(m.group(1)), 3), round(float(m.group(2)), 3))
            elif line.startswith("ndd (LSST)") and key:
                nl = int(line.split()[2])
            elif line.startswith("ndd (Roman)") and key:
                out[key] = (nl, int(line.split()[2]))
                key = None
    return out


def calibrate(img):
    """Pixel -> (l, b) from the grid lines: returns (l(x), b(y)) as linear functions."""
    R, G, B = (img[..., k].astype(int) for k in range(3))
    green = (G > R + 20) & (G > B + 20)
    H, W = green.shape
    cols = [x for x in range(W) if green[:, x].sum() > 0.9 * H]
    rows = [y for y in range(H) if green[y, :].sum() > 0.9 * W]
    # Merge adjacent pixels of the same line, keep the centre.
    def centres(v):
        v, groups = sorted(v), []
        for p in v:
            if groups and p - groups[-1][-1] <= 2:
                groups[-1].append(p)
            else:
                groups.append([p])
        return [float(np.mean(g)) for g in groups]
    xs, ys = centres(cols), centres(rows)
    step_x, step_y = np.diff(xs).mean(), np.diff(ys).mean()
    if max(np.ptp(np.diff(xs)), np.ptp(np.diff(ys))) > 3 or abs(step_x - step_y) > 3:
        sys.exit(f"grid lines not regular: x {xs}, y {ys}")
    # Anchor: the column nearest x = 615 is l = 0; the row nearest y = 124 is b = 0 (checked
    # against Aladin's own labels, '000' and '+00').
    x0 = min(xs, key=lambda x: abs(x - 615))
    y0 = min(ys, key=lambda y: abs(y - 124))
    ppd = 0.5 * (step_x + step_y) / 0.5          # pixels per degree (lines every 0.5 deg)
    lfun = lambda x: -(np.asarray(x) - x0) / ppd
    bfun = lambda y: -(np.asarray(y) - y0) / ppd
    return lfun, bfun, ppd, xs, ys


def real_tiles(img):
    """Boolean mask of pixels inside a real GBTDS detector outline (filled)."""
    R, G, B = (img[..., k].astype(int) for k in range(3))
    edge = (G > 170) & (R < 190) & (B < 150) & (G > R + 40)
    H, W = edge.shape
    edge[:45, :] = edge[790:, :] = False            # Aladin's labels and scale bar
    edge[:, :45] = edge[:, 935:] = False
    # Dilate to close the outlines, fill, then erode the dilation back off: without the erosion
    # every detector comes out ~9% too large (a 1-px rim on a ~45-px square).
    edge = ndimage.binary_dilation(edge, iterations=1)
    return ndimage.binary_erosion(ndimage.binary_fill_holes(edge), iterations=1)


def main():
    if "--legacy" not in sys.argv:
        return main_gbtds()
    os.makedirs(OUT, exist_ok=True)
    ps.use_paper_style()

    # ---- rebuild and verify the scan ----
    cells = tiling()
    reps = {(round(r[0], 3), round(r[1], 3)) for _, _, r, _ in cells}
    logged = logged_sightlines()
    if reps != set(logged):
        sys.exit(f"rebuilt tiling does not match the run: {len(reps)} rebuilt, "
                 f"{len(logged)} logged, {len(reps ^ set(logged))} differ")
    foot = {(round(c[0], 3), round(c[1], 3)) for c in cells if c[3]}
    roman_logged = {k for k, (nl, nr) in logged.items() if nr > 0}
    if foot != roman_logged:
        sys.exit("rebuilt footprint differs from the sightlines with Roman epochs")
    area = len(cells) * FINE ** 2
    print(f"tiling verified: {len(reps)} sightlines, {len(foot)} in the footprint, "
          f"{area:.2f} deg^2 scanned, {len(foot) * FINE**2:.2f} deg^2 footprint")
    n_none = sum(1 for v in logged.values() if v[0] == 0 and v[1] == 0)
    n_rub = sum(1 for v in logged.values() if v[0] > 0)
    area_rub = sum(FINE**2 for c in cells
                   if logged[(round(c[2][0], 3), round(c[2][1], 3))][0] > 0)
    print(f"Rubin coverage: {n_rub} sightlines, {area_rub:.2f} deg^2; no coverage: {n_none}")

    # ---- figure 1: the simulated footprints ----
    fig, ax = ps.figure(width="double", height=5.6)
    nvis = np.array([logged[(round(c[2][0], 3), round(c[2][1], 3))][0] for c in cells])
    covered = [Rectangle((c[0] - HALF, c[1] - HALF), FINE, FINE)
               for c, n in zip(cells, nvis) if n > 0]
    empty = [Rectangle((c[0] - HALF, c[1] - HALF), FINE, FINE) for c, n in zip(cells, nvis) if n == 0]
    pc = PatchCollection(covered, cmap="Blues", norm=LogNorm(vmin=max(nvis[nvis > 0].min(), 1),
                                                             vmax=nvis.max()),
                         edgecolor="face", linewidth=0.3)   # own-colour edges: no PDF seams
    pc.set_array(nvis[nvis > 0])
    ax.add_collection(pc)
    ax.add_collection(PatchCollection(empty, facecolor="#eeeeee", edgecolor="#eeeeee",
                                      linewidth=0.3))
    ax.add_collection(PatchCollection(
        [Rectangle((c[0] - HALF, c[1] - HALF), FINE, FINE) for c in cells if c[3]],
        facecolor="#f97316", edgecolor="white", linewidth=0.2))
    for fl, fb in FIELDS:
        ax.add_patch(Circle((fl, fb), FOV_ROMAN, fill=False, ls="--", lw=0.6, color=ps.INK))
    sl = np.array(sorted(logged))
    ax.plot(sl[:, 0], sl[:, 1], ".", ms=0.8, color=ps.INK, alpha=0.5)
    cb = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("Rubin visits at the sightline")
    ax.set_xlim(LON_MAX + 0.2, LON_MIN - 0.1)             # l increases to the left
    ax.set_ylim(LAT_MIN - 0.1, LAT_MAX + 0.2)
    ax.set_aspect("equal")
    ax.set_xlabel(r"Galactic longitude $l$ [deg]")
    ax.set_ylabel(r"Galactic latitude $b$ [deg]")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Patch(facecolor="#f97316", label=f"Roman, as simulated ({len(foot)})"),
               Line2D([], [], ls="--", color=ps.INK, lw=0.6, label="modelled GBTDS fields"),
               Patch(facecolor=plt.cm.Blues(0.6), label=f"Rubin coverage ({n_rub:,})"),
               Patch(facecolor="#eeeeee", edgecolor="#cccccc", label=f"no coverage ({n_none})"),
               Line2D([], [], ls="", marker=".", ms=3, color=ps.INK, label="sightline")]
    # Upper right is the corner the scan cuts away (l < -0.94, b > 0.31): empty sky, room for it.
    ax.legend(handles=handles, loc="upper right", fontsize=6, frameon=False, labelcolor=ps.INK,
              title="sightlines", title_fontsize=6, borderaxespad=0.3)
    ps.stamp(fig, f"rebuilt from Bulge.h's scan and checked against {LOG}: {len(reps)} sightlines, "
                  f"{area:.2f} deg^2 scanned, {len(foot) * FINE**2:.2f} deg^2 Roman footprint")
    for p in ps.save_figure(fig, f"{OUT}/footprint_simulated"):
        print("wrote", p)
    plt.close(fig)

    # ---- figure 2: simulated Roman footprint over the real GBTDS tiles ----
    img = np.asarray(Image.open(IMG).convert("RGB"))
    lfun, bfun, ppd, xs, ys = calibrate(img)
    H, W = img.shape[:2]
    print(f"calibration: {ppd:.1f} px/deg; grid columns {[round(x) for x in xs]}, rows "
          f"{[round(y) for y in ys]}")
    ext = [float(lfun(-0.5)), float(lfun(W - 0.5)), float(bfun(H - 0.5)), float(bfun(-0.5))]
    # Real detectors: the union of the SEPARATE spring and autumn tile images (same framing, checked
    # below). Filling outlines on the combined image instead also fills the small pockets where the
    # two roll angles' tiles cross, which inflates the area by ~3%.
    tiles = np.zeros(img.shape[:2], bool)
    for season in ("spring", "autumn"):
        simg = np.asarray(Image.open(IMG.replace("both", season)).convert("RGB"))
        l2, b2, ppd2, xs2, ys2 = calibrate(simg)
        if (xs2, ys2) != (xs, ys):
            sys.exit(f"{season} image is framed differently from the combined one")
        t = real_tiles(simg)
        print(f"  {season} tiles: {t.sum() / ppd2**2:.3f} deg^2 (design: 1.7)")
        tiles |= t
    area_real = tiles.sum() / ppd**2
    # Overlap, on the image's own pixel grid.
    yy, xx = np.mgrid[0:H, 0:W]
    L, B = lfun(xx), bfun(yy)
    sim = np.zeros_like(tiles)
    for c in cells:
        if c[3]:
            sim |= ((L >= c[0] - HALF) & (L < c[0] + HALF)
                    & (B >= c[1] - HALF) & (B < c[1] + HALF))
    area_sim = sim.sum() / ppd**2
    both = (sim & tiles).sum() / ppd**2
    print(f"real tiles (spring+autumn union, in view): {area_real:.2f} deg^2; simulated footprint "
          f"in view: {area_sim:.2f} deg^2; overlap {both:.2f} deg^2 = {100*both/area_sim:.0f}% of "
          f"the simulated, {100*both/area_real:.0f}% of the real")
    for name, m in (("GC", B > -0.8), ("five-field", B <= -0.8)):
        t = tiles & m
        print(f"  real {name} block: l {L[t].min():.2f}..{L[t].max():.2f}, "
              f"b {B[t].min():.2f}..{B[t].max():.2f}, centre ({L[t].mean():.2f}, {B[t].mean():.2f})")

    fig, ax = ps.figure(width="double", height=5.9)
    ax.imshow(img, extent=ext, origin="upper", interpolation="bilinear")
    ax.add_collection(PatchCollection(
        [Rectangle((c[0] - HALF, c[1] - HALF), FINE, FINE) for c in cells if c[3]],
        facecolor="#f97316", alpha=0.35, edgecolor="#f97316", linewidth=0.5))
    for fl, fb in FIELDS:
        ax.add_patch(Circle((fl, fb), FOV_ROMAN, fill=False, ls="--", lw=1.0, color="#ffb070"))
        ax.plot(fl, fb, "+", color="#ffb070", ms=5, mew=1.0)
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_aspect("equal")
    ax.set_xlabel(r"Galactic longitude $l$ [deg]")
    ax.set_ylabel(r"Galactic latitude $b$ [deg]")
    handles = [Patch(facecolor="#f97316", alpha=0.5, edgecolor="#f97316",
                     label="Roman footprint as simulated (0.1$^\\circ$ cells)"),
               Line2D([], [], ls="--", color="#ffb070", lw=1.0,
                      label="modelled field circles and centres"),
               Patch(facecolor="none", edgecolor="#7CFC00",
                     label="GBTDS detectors, spring and autumn (Aladin)")]
    # Above the axes: over the image no legend text stays legible.
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3,
              fontsize=6, frameon=False, labelcolor=ps.INK)
    ps.stamp(fig, f"background: {IMG}, calibrated from its Galactic grid ({ppd:.0f} px/deg). "
                  f"Overlap {both:.2f} deg^2: {100*both/area_sim:.0f}% of the simulated footprint, "
                  f"{100*both/area_real:.0f}% of the real tiles")
    for p in ps.save_figure(fig, f"{OUT}/footprint_vs_gbtds"):
        print("wrote", p)
    plt.close(fig)


def main_gbtds():
    """The adopted layout: rebuilt with analysis/gbtds_geometry, the simulator's mirror."""
    import gbtds_geometry as G
    out = "figures/footprint_20261001"
    os.makedirs(out, exist_ok=True)
    ps.use_paper_style()
    g = G.scan_sightlines(10, 5)                     # the production flags' grid
    lon, lat, fine = g["lon"], g["lat"], g["fine"]
    fstep = g["fine_step"]
    # Rubin visits per sightline, with the simulator's circle test.
    rub = np.loadtxt("Baseline/BulgeBaseline.dat", comments="#", usecols=(3, 4, 5))
    nvis = np.array([np.unique(rub[np.hypot(rub[:, 0] - a, rub[:, 1] - b) <= G.FOV_RUBIN, 2]).size
                     for a, b in zip(lon, lat)])
    if "--log" in sys.argv:
        global LOG
        LOG = sys.argv[sys.argv.index("--log") + 1]
        logged = logged_sightlines()
        mine = {(round(a, 3), round(b, 3)) for a, b in zip(lon, lat)}
        if mine != set(logged):
            sys.exit(f"rebuilt scan differs from {LOG}: {len(mine ^ set(logged))} sightlines")
        bad = sum(logged[(round(a, 3), round(b, 3))][0] != n for a, b, n in zip(lon, lat, nvis))
        print(f"scan verified against {LOG}: {len(mine)} sightlines; Rubin counts differ at {bad}")
    od = G.on_detector(lon, lat)
    print(f"scan: {len(lon)} sightlines, {g['area'].sum():.2f} deg^2; footprint stratum "
          f"{fine.sum()} at {fstep} deg; on a detector: spring {od[fine, 0].sum()}, "
          f"autumn {od[fine, 1].sum()}, both {(od[fine, 0] & od[fine, 1]).sum()}")

    # ---- figure 1: the scan, Rubin visits, Roman detectors ----
    fig, ax = ps.figure(width="double", height=5.6)
    cells = g["cells"]
    vis = np.array([nvis[n] for _, _, n in cells])
    covered = [Rectangle((c[0], c[1]), fstep, fstep) for c, v in zip(cells, vis) if v > 0]
    pc = PatchCollection(covered, cmap="Blues", norm=LogNorm(vmin=max(vis[vis > 0].min(), 1),
                                                             vmax=vis.max()),
                         edgecolor="face", linewidth=0.3)
    pc.set_array(vis[vis > 0])
    ax.add_collection(pc)
    empty = [Rectangle((c[0], c[1]), fstep, fstep) for c, v in zip(cells, vis) if v == 0]
    ax.add_collection(PatchCollection(empty, facecolor="#eeeeee", edgecolor="#eeeeee", linewidth=0.3))
    colours = {0: "#f97316", 1: "#7c3aed"}
    for l0, b0, k in G.placements():
        for dl0, dl1, db0, db1 in G.sca_rects(k):
            ax.add_patch(Rectangle((l0 + dl0, b0 + db0), dl1 - dl0, db1 - db0, fill=False,
                                   lw=0.4, ec=colours[k]))
    th = np.linspace(0, 2 * np.pi, 400)
    for l0, b0, _ in G.placements():
        ax.plot(l0 + G.scan_reach() * np.cos(th), b0 + G.scan_reach() * np.sin(th), lw=0.15,
                color=ps.INK, alpha=0.3)
    cb = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("Rubin visits at the sightline")
    lo, hi, blo, bhi = g["bounds"]
    ax.set_xlim(hi + 0.1, lo - 0.1)
    ax.set_ylim(blo - 0.1, bhi + 0.1)
    ax.set_aspect("equal")
    ax.set_xlabel(r"Galactic longitude $l$ [deg]")
    ax.set_ylabel(r"Galactic latitude $b$ [deg]")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fill=False, ec=colours[0], label="Roman detectors, spring roll"),
                       Patch(fill=False, ec=colours[1], label="Roman detectors, autumn roll"),
                       Patch(facecolor=plt.cm.Blues(0.6), label="Rubin coverage"),
                       Patch(facecolor="#eeeeee", label="no coverage"),
                       Line2D([], [], lw=0.4, color=ps.INK, alpha=0.5,
                              label=f"scan reach ({G.scan_reach():.2f} deg)")],
              loc="upper right", fontsize=6, frameon=False, labelcolor=ps.INK)
    ps.stamp(fig, f"analysis/gbtds_geometry scan_sightlines(10, 5): {len(lon)} sightlines, "
                  f"{g['area'].sum():.2f} deg^2; layout gbtds_2026.4.3 (Baseline/gbtds_layout)")
    for p in ps.save_figure(fig, f"{out}/footprint_simulated"):
        print("wrote", p)
    plt.close(fig)

    # ---- figure 2: simulated detectors over the real tiles, and their overlap ----
    img = np.asarray(Image.open(IMG).convert("RGB"))
    lfun, bfun, ppd, xs, ys = calibrate(img)
    H, W = img.shape[:2]
    ext = [float(lfun(-0.5)), float(lfun(W - 0.5)), float(bfun(H - 0.5)), float(bfun(-0.5))]
    yy, xx = np.mgrid[0:H, 0:W]
    L, B = lfun(xx), bfun(yy)
    print(f"calibration: {ppd:.1f} px/deg")
    rows = []
    union_real = np.zeros((H, W), bool)
    union_sim = np.zeros((H, W), bool)
    for k, season in ((0, "spring"), (1, "autumn")):
        simg = np.asarray(Image.open(IMG.replace("both", season)).convert("RGB"))
        if calibrate(simg)[3:] != (xs, ys):
            sys.exit(f"{season} image is framed differently from the combined one")
        real = real_tiles(simg)
        sim = G.on_detector(L.ravel(), B.ravel())[:, k].reshape(H, W)
        both = (real & sim).sum() / ppd ** 2
        rows.append((season, real.sum() / ppd ** 2, sim.sum() / ppd ** 2, both))
        print(f"  {season}: real tiles {real.sum() / ppd**2:.3f} deg^2, simulated detectors "
              f"{sim.sum() / ppd**2:.3f}, overlap {both:.3f} = {100*both/(sim.sum()/ppd**2):.1f}% "
              f"of simulated, {100*both/(real.sum()/ppd**2):.1f}% of real")
        union_real |= real
        union_sim |= sim
    both = (union_real & union_sim).sum() / ppd ** 2
    ar, asim = union_real.sum() / ppd ** 2, union_sim.sum() / ppd ** 2
    print(f"  union: real {ar:.3f}, simulated {asim:.3f}, overlap {both:.3f} = "
          f"{100*both/asim:.1f}% / {100*both/ar:.1f}%")
    fig, ax = ps.figure(width="double", height=5.9)
    ax.imshow(img, extent=ext, origin="upper", interpolation="bilinear")
    for l0, b0, k in G.placements():
        for dl0, dl1, db0, db1 in G.sca_rects(k):
            ax.add_patch(Rectangle((l0 + dl0, b0 + db0), dl1 - dl0, db1 - db0, fill=False,
                                   lw=0.6, ec=colours[k], ls="--" if k else "-"))
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    ax.set_aspect("equal")
    ax.set_xlabel(r"Galactic longitude $l$ [deg]")
    ax.set_ylabel(r"Galactic latitude $b$ [deg]")
    ax.legend(handles=[Patch(fill=False, ec=colours[0], label="simulated detectors, spring"),
                       Patch(fill=False, ec=colours[1], ls="--", label="simulated detectors, autumn"),
                       Patch(facecolor="none", edgecolor="#7CFC00",
                             label="GBTDS detectors, spring and autumn (Aladin)")],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=6, frameon=False,
              labelcolor=ps.INK)
    ps.stamp(fig, f"background: {IMG} ({ppd:.0f} px/deg). Overlap (union of rolls) {both:.2f} "
                  f"deg^2: {100*both/asim:.0f}% of the simulated detectors, {100*both/ar:.0f}% of "
                  f"the real tiles")
    for p in ps.save_figure(fig, f"{out}/footprint_vs_gbtds"):
        print("wrote", p)
    plt.close(fig)


if __name__ == "__main__":
    main()
