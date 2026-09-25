#!/usr/bin/env python3
"""The Besancon-style density model of `Disk_model()`, ported to NumPy (Step W1).

WHY THIS EXISTS. The pooled event weight (Deviation 41) needs two per-sightline quantities that
the simulator computes but does not write per event: the number of source stars the sightline
holds, and the normaliser of `func_lens`'s lens-distance sampler. The first is in the map file
as `log10 Nstart` (one decimal, and the v3 file is missing rows); the second is nowhere. Both are
deterministic functions of (l, b), so they are recomputed here rather than stored.

THIS IS A PORT, NOT A SECOND MODEL. Every constant mirrors `Bulge.h` / `Bulge_LSST.cpp`. If the
C++ changes, this must change with it -- `check_against_map()` is the guard: it compares the
recomputed column densities against the map file's own and is the first thing to run after any
change to `Disk_model`.

    from galaxy_model import density_profile, lens_distance_norm
    prof = density_profile(l, b)
    Z = lens_distance_norm(prof, Ds_array)
"""

import numpy as np

# ---- mirrored from Bulge.h ----
NUM = 9500
MAXD = 12.0                      # kpc
STEP = MAXD / NUM                # kpc
DSUN = 8.0                       # kpc
RAA = 180.0 / np.pi
BINARY_FRACTION = 2.0 / 3.0
RHO0 = np.array([4.0, 7.9, 6.2, 4.0, 5.8, 4.9, 6.6, 3.96])
D0 = np.array([0.073117, 0.0216524, 0.0217405, 0.0217901,
               0.0218061, 0.0218118, 0.0218121, 0.0218121])
EPCI = np.array([0.014, 0.0268, 0.0375, 0.0551, 0.0696, 0.0785, 0.0791, 0.0791])
CORR = np.array([1.0, 7.9 / 4.48419, 6.2 / 3.52112, 4.0 / 2.27237,
                 5.8 / 3.29525, 4.9 / 2.78402, 6.6 / 3.74991, 3.96 / 2.24994])
BAR_MASS_RESCALE = 0.24529       # calibrated to the Han & Gould (2003) Baade's Window benchmark

# Mean stellar mass per component, the divisors in Nstari. Thin disk, thick disk and halo share
# the disk value in the C++; the bulge is lighter.
MBAR_THIN, MBAR_THICK, MBAR_HALO, MBAR_BULGE = 0.403445, 0.4542, 0.4542, 0.308571


class Profile:
    """Densities along one sightline, on the same grid the simulator integrates."""

    __slots__ = ("lon", "lat", "x", "rho_tot", "Rostart", "Nstart")

    def __init__(self, lon, lat, x, rho_tot, Rostart, Nstart):
        self.lon, self.lat = lon, lat
        self.x = x                   # kpc, distance along the sightline
        self.rho_tot = rho_tot       # Msun/pc^3, all four components
        self.Rostart = Rostart       # Msun/deg^2, the column
        self.Nstart = Nstart         # stars/deg^2, the column


def density_profile(lon, lat):
    """Port of `Disk_model(s, numt)`: densities and column totals toward (lon, lat) in degrees."""
    TET, FI = (360.0 - lon) / RAA, lat / RAA
    x = np.arange(1, NUM) * STEP
    zb = np.sin(FI) * x
    yb = np.cos(FI) * np.sin(TET) * x
    xb = x * np.cos(FI) * np.cos(TET) - DSUN
    Rb = np.hypot(xb, yb)

    # --- thin disk: eight age bins, each with its own scale height ---
    Rdd, Rhh = 2.17, 1.33
    thin = np.zeros_like(x)
    for ii in range(8):
        rdi = Rb**2 + zb**2 / EPCI[ii]**2
        if ii == 0:
            rho = np.exp(-rdi / 25.0) - np.exp(-rdi / 9.0)
        else:
            rho = (np.exp(-np.sqrt(0.25 + rdi / Rdd**2))
                   - np.exp(-np.sqrt(0.25 + rdi / Rhh**2)))
        rho *= 1.2                                   # total mass 4.25e10 Msun
        thin = np.abs(thin + RHO0[ii] * CORR[ii] * 0.001 * rho / D0[ii])

    # --- thick disk: parabolic near the plane, exponential beyond 400 pc ---
    nnf = 0.4 / 0.8
    rho00 = 1.34e-3 + 3.04e-4
    near = np.abs((rho00 / 0.999719) * np.exp(-(Rb - DSUN) / 2.5)
                  * (1.0 - zb**2 / (0.4 * 0.8 * (2.0 + nnf))))
    far = np.abs((rho00 / 0.999719) * np.exp(-(Rb - DSUN) / 2.5)
                 * np.exp(nnf) * np.exp(-np.abs(zb) / 0.8) / (1.0 + 0.5 * nnf))
    thick = np.where(np.abs(zb) < 0.4, near, far) * 2.67   # total mass 0.8e10 Msun

    # --- stellar halo: flattened power law, cored inside 0.5 kpc ---
    rdi = np.sqrt(Rb**2 + zb**2 / 0.76**2)
    halo = np.abs((0.932e-5 / 867.067) * np.power(np.maximum(rdi, 0.5) / DSUN, -2.44)) * 5281.0

    # --- bulge: two triaxial bars, rotated by alfa from the Sun-centre line ---
    alfa = 12.89 / RAA
    xf = xb * np.cos(alfa) + yb * np.sin(alfa)
    yf = -xb * np.sin(alfa) + yb * np.cos(alfa)
    zf = zb

    def _bar_radii(Rx0, Ry0, Rz0, cp, cn):
        r4 = (np.power(np.abs(np.power(np.abs(xf / Rx0), cn)
                              + np.power(np.abs(yf / Ry0), cn)), cp / cn)
              + np.power(np.abs(zf / Rz0), cp))
        return np.power(np.abs(r4), 1.0 / cp), np.hypot(xf, yf)

    r4, r2 = _bar_radii(1.46, 0.49, 0.39, 3.007, 3.329)
    mS = 35.45 / 3.84723 * BAR_MASS_RESCALE
    rhoS = mS / np.cosh(-r4)**2 * np.where(r2 <= 3.43, 1.0, np.exp(-4.0 * (r2 - 3.43)**2))

    r4, r2 = _bar_radii(4.44, 1.31, 0.80, 2.786, 3.917)
    mE = 2.27 / 87.0 * BAR_MASS_RESCALE
    rhoE = mE * np.exp(-r4) * np.where(r2 <= 6.83, 1.0, np.exp(-4.0 * (r2 - 6.83)**2))

    bulge = (np.abs(rhoS) + np.abs(rhoE)) * 0.45           # total mass 1.7e10 Msun

    rho_tot = thin + thick + bulge + halo
    # Shell volume per square degree: x^2 dx in kpc^3 -> pc^3, times (pi/180)^2 sr per deg^2.
    shell = x * x * STEP * 1.0e9 * (np.pi / 180.0)**2
    Rostart = float(np.sum(rho_tot * shell))
    Nstart = float(np.sum(BINARY_FRACTION * (thin / MBAR_THIN + thick / MBAR_THICK
                                             + halo / MBAR_HALO + bulge / MBAR_BULGE) * shell))
    return Profile(lon, lat, x, rho_tot, Rostart, Nstart)


def lens_distance_norm(prof, Ds):
    """Normaliser Z(Ds) of `func_lens`'s lens-distance sampler, for each source distance in Ds.

    `func_lens` accepts a lens at Dl with probability proportional to
    `rho(Dl) * sqrt((Ds - Dl) Dl / Ds)` over the grid points k = 1 .. nums-2. Z is the sum of
    that, so dividing the physical rate by the sampling density leaves Z as a factor -- see
    Deviation 41. Units are Msun/pc^3 * kpc^(3/2), which cancels in any weighted fraction.
    """
    Ds = np.atleast_1d(np.asarray(Ds, dtype=float))
    k = np.arange(1, NUM)[None, :]
    Dl = k * STEP
    out = np.empty(len(Ds))
    # In row blocks: the (rows x grid) matrix and its temporaries are ~0.4 MB per row, so one
    # post-extinction-fix sightline (thousands of draws) in a single pass took GBs and got y1
    # killed for memory (2026-09-25). Each row's sum is unchanged by the blocking.
    for i in range(0, len(Ds), BLOCK):
        d = Ds[i:i + BLOCK, None]
        nums = np.round(d / STEP).astype(int)
        inside = k <= nums - 2
        tt = np.sqrt(np.clip((d - Dl) * Dl / d, 0.0, None)) * prof.rho_tot[None, :] * inside
        out[i:i + BLOCK] = tt.sum(axis=1) * STEP
    return out


BLOCK = 512


def check_against_map(sightlines, tol=0.05):
    """Recompute every sightline's column densities and compare with the map file's own.

    The map file writes `log10 Rostart` and `log10 Nstart` to ONE decimal, so agreement is only
    ever testable at the 0.05 dex level -- that is the guard's resolution, not its precision.
    Returns (worst_dex, n_checked). Raises if the port has drifted from the C++.
    """
    worst, n = 0.0, 0
    for _, row in sightlines.iterrows():
        if not np.isfinite(row.get("lon", np.nan)):
            continue
        prof = density_profile(row["lon"], row["lat"])
        worst = max(worst,
                    abs(np.log10(prof.Rostart) - row["log10_Rostart"]),
                    abs(np.log10(prof.Nstart) - row["log10_Nstart"]))
        n += 1
    if n and worst > tol:
        raise AssertionError(
            f"galaxy_model.py disagrees with the map file by {worst:.3f} dex over {n} sightlines "
            f"(tolerance {tol}). The port has drifted from Disk_model() in Bulge_LSST.cpp.")
    return worst, n
