#!/usr/bin/env python3
"""The Besancon-style density model of `Disk_model()`, ported to NumPy.

The pooled event weight needs two per-sightline quantities that the simulator computes but does
not write per event: the number of source stars the sightline holds, and the normaliser of
`func_lens`'s lens-distance sampler. The first is in the map file as `log10 Nstart` (one decimal);
the second is not stored. Both are deterministic functions of (l, b), so they are recomputed here.

This is a port, not a second model: every constant mirrors `config/parameters.h`,
`src/galaxy/density.cpp` and `src/galaxy/kinematics.cpp`, and must change with them.
`check_against_map()` compares the recomputed column densities with the map file's own and is the
first thing to run after any change to `Disk_model`.

Three galaxy models exist in runs on disk (the GALAXIES table below). Every profile therefore takes the
run's `galaxy` key, and `galaxy_for_map()` reads which one a run used from its own map file:
  "current"   bar once + Besancon 1612's thin disc, thick disc and halo, in shape and local density, the Sun 15 pc
              above the plane (run_provenance.txt has a disc_norm line)
  "bar-once"  bar once + Robin et al. 2003's discs and halo (runs of 2026-11-03 on; bulge_norm line, no thin_disc)
  "legacy"    bar x 0.45 + Robin 2003's discs and halo (runs before 2026-11-03; no bulge_norm line)
Runs of 2026-11-07 (9170d33: a thin_disc line but no disc_norm line) had Besancon 1612's thin-disc shape with Robin
2003's normalisations; none was kept, so that model is not in the table and such a map matches nothing.
The kinematics in kinematics.cpp are not ported here: the weights read each event's own Vt column.

    from galaxy_model import galaxy_for_map, density_profile, lens_distance_norm
    prof = density_profile(l, b, galaxy=galaxy_for_map(sightlines))
    Z = lens_distance_norm(prof, Ds_array)
"""

import numpy as np

from cparams import P            # the C++ headers' own numbers (analysis/cparams.py)

# ---- read from config/parameters.h ----
NUM = P.Num
MAXD = P.MaxD                    # kpc
STEP = MAXD / NUM                # kpc
DSUN = P.Dsun                    # kpc
RAA = 180.0 / np.pi
BINARY_FRACTION = P.binary_fraction
# Besancon 1612 (the "current" galaxy): thin disc per age bin, thick disc's two episodes, halo.
THIN_RHO_LOCAL, EPCI = P.THIN_RHO_LOCAL, P.epci
THICK_RHO_LOCAL, THICK_HZ, THICK_HR = P.THICK_RHO_LOCAL, P.THICK_HZ, P.THICK_HR
HALO_FLATTEN, HALO_A_CORE, HALO_EXPONENT = P.HALO_FLATTEN, P.HALO_A_CORE, P.HALO_EXPONENT
HALO_A_MIN, HALO_RHO_LOCAL = P.HALO_A_MIN, P.HALO_RHO_LOCAL
BAR_MASS_RESCALE = P.BAR_MASS_RESCALE   # calibrated to the Han & Gould (2003) Baade's Window benchmark
Z_SUN = P.Z_SUN                         # kpc, the Sun's height above the plane (Besancon 1612)
# The older runs put the Sun in the plane, and their bar scale was calibrated that way.
BAR_MASS_RESCALE_SUN0 = 0.24529

# The older runs' thin disc ("bar-once", "legacy"): Robin et al. (2003) as the C++ had it until 2026-11-07, eight bins
# (the 8th: disc white dwarfs), local densities RHO0_ROBIN03 [1e-3 Msun/pc^3] x CORR_ROBIN03 x 1.2 / D0_ROBIN03, and
# Robin 2003's flattenings. Their hole term also divided the vertical term by the hole's scale length squared, which
# Besancon 1612 does not. config/parameters.h no longer has these; the analysis keeps them to weight those runs.
RHO0_ROBIN03 = (4.0, 7.9, 6.2, 4.0, 5.8, 4.9, 6.6, 3.96)
D0_ROBIN03 = (0.073117, 0.0216524, 0.0217405, 0.0217901, 0.0218061, 0.0218118, 0.0218121, 0.0218121)
CORR_ROBIN03 = (1.0, 7.9 / 4.48419, 6.2 / 3.52112, 4.0 / 2.27237, 5.8 / 3.29525, 4.9 / 2.78402, 6.6 / 3.74991,
                3.96 / 2.24994)
EPCI_ROBIN03 = (0.014, 0.0268, 0.0375, 0.0551, 0.0696, 0.0785, 0.0791, 0.0791)

# The bar's second factor in "legacy" runs ("total mass 1.7e10 Msun", computed without BAR_MASS_RESCALE, so the
# bar was normalised twice: 937 instead of 2,086 Msun/pc^2 toward Baade's window). config/parameters.h no
# longer has it; the analysis keeps it to weight those runs.
# key: (bar mass rescale, bar factor, discs and halo, the Sun's height [kpc])
GALAXIES = {
    "current":  (BAR_MASS_RESCALE, 1.0, "besancon1612", Z_SUN),
    "bar-once": (BAR_MASS_RESCALE_SUN0, 1.0, "robin2003", 0.0),
    "legacy":   (BAR_MASS_RESCALE_SUN0, 0.45, "robin2003", 0.0),
}

# Mean stellar mass per component, the divisors in Nstari: config/data_products.h's MEANMASS_*, the
# mean mass of each complete CMD population (provenance.txt, mean_mass_population). A run made with
# different CMD lists carries its own Nstart in its map file, which is what the weights read.
MBAR_THIN, MBAR_THICK, MBAR_HALO, MBAR_BULGE = P.MEANMASS_THIN, P.MEANMASS_THICK, P.MEANMASS_HALO, P.MEANMASS_BULGE


class Profile:
    """Densities along one sightline, on the same grid the simulator integrates."""

    __slots__ = ("lon", "lat", "x", "rho_tot", "Rostart", "Nstart")

    def __init__(self, lon, lat, x, rho_tot, Rostart, Nstart):
        self.lon, self.lat = lon, lat
        self.x = x                   # kpc, distance along the sightline
        self.rho_tot = rho_tot       # Msun/pc^3, all four components
        self.Rostart = Rostart       # Msun/deg^2, the column
        self.Nstart = Nstart         # stars/deg^2, the column


def density_profile(lon, lat, galaxy):
    """Port of `Disk_model(s, numt)`: densities and column totals toward (lon, lat) in degrees.

    galaxy: a key of GALAXIES; galaxy_for_map() says which a run used. There is no default, so a run is
    never weighted with another run's galaxy.
    """
    return _profile(_components(lon, lat, GALAXIES[check_galaxy(galaxy)][3]), galaxy)


def check_galaxy(galaxy):
    if galaxy not in GALAXIES:
        raise ValueError(
            f"galaxy must be one of {tuple(GALAXIES)}, got {galaxy!r}: 'current' if the run's run_provenance.txt "
            f"has a disc_norm line, 'bar-once' if it has a bulge_norm line but no thin_disc, 'legacy' if it has "
            f"neither. galaxy_for_map(sightlines) reads it from the run's map file.")
    return galaxy


def _profile(comp, galaxy):
    """A Profile from _components()'s output and the run's galaxy key."""
    lon, lat, x, discs, bar = comp
    bar_rescale, bar_factor, disc, _ = GALAXIES[check_galaxy(galaxy)]
    thin, thick, halo = discs[disc]
    bulge = bar * bar_rescale * bar_factor
    rho_tot = thin + thick + bulge + halo
    # Shell volume per square degree: x^2 dx in kpc^3 -> pc^3, times (pi/180)^2 sr per deg^2.
    shell = x * x * STEP * 1.0e9 * (np.pi / 180.0)**2
    Rostart = float(np.sum(rho_tot * shell))
    Nstart = float(np.sum(BINARY_FRACTION * (thin / MBAR_THIN + thick / MBAR_THICK
                                             + halo / MBAR_HALO + bulge / MBAR_BULGE) * shell))
    return Profile(lon, lat, x, rho_tot, Rostart, Nstart)


def _thin_bin_shape(ii, Rb, zb):
    """Shape of Besancon-1612 thin-disc age bin ii at (R, z) [kpc]; the hole has the disc's vertical scale."""
    e2 = EPCI[ii]**2
    if ii == 0:
        return (np.exp(-(Rb**2 + zb**2 / e2) / 25.0)
                - np.exp(-(Rb**2 / 9.0 + zb**2 / (e2 * 25.0))))
    return (np.exp(-np.sqrt(0.25 + (Rb**2 + zb**2 / e2) / 2.17**2))
            - np.exp(-np.sqrt(0.25 + Rb**2 / 1.33**2 + zb**2 / (e2 * 2.17**2))))


def thin_bin_density(ii, Rb, zb):
    """Port of `thinBinDensity()`: Msun/pc^3 of Besancon-1612 thin-disc age bin ii at (R, z) [kpc], the bin's
    local density THIN_RHO_LOCAL at (DSUN, 0)."""
    return THIN_RHO_LOCAL[ii] * _thin_bin_shape(ii, Rb, zb) / _thin_bin_shape(ii, DSUN, 0.0)


def thick_density(Rb, zb):
    """Port of `thickDensity()`: Msun/pc^3 of Besancon 1612's two sech^2 thick-disc episodes."""
    return sum(THICK_RHO_LOCAL[k] * np.exp(-(Rb - DSUN) / THICK_HR[k]) / np.cosh(zb / (2.0 * THICK_HZ[k]))**2
               for k in range(len(THICK_RHO_LOCAL)))


def halo_density(Rb, zb):
    """Port of `haloDensity()`: Msun/pc^3 of Besancon 1612's spheroid, 1 / (a (a_c + a)^n), a >= HALO_A_MIN."""
    def law(a):
        return 1.0 / (a * (HALO_A_CORE + a)**HALO_EXPONENT)
    a = np.maximum(np.sqrt(Rb**2 + zb**2 / HALO_FLATTEN**2), HALO_A_MIN)
    return HALO_RHO_LOCAL * law(a) / law(DSUN)


def _components(lon, lat, zsun):
    """The densities [Msun/pc^3] along (lon, lat), with the Sun `zsun` kpc above the plane: {model: (thin, thick,
    halo)} for both models (Besancon 1612, Robin 2003), and the bar before its mass rescale and any bar factor."""
    TET, FI = (360.0 - lon) / RAA, lat / RAA
    x = np.arange(1, NUM) * STEP
    zb = np.sin(FI) * x + zsun
    yb = np.cos(FI) * np.sin(TET) * x
    xb = x * np.cos(FI) * np.cos(TET) - DSUN
    Rb = np.hypot(xb, yb)

    # --- Besancon 1612: thin disc (seven age bins), thick disc (two episodes), halo ---
    thin_b = sum(thin_bin_density(ii, Rb, zb) for ii in range(len(THIN_RHO_LOCAL)))
    b1612 = (thin_b, thick_density(Rb, zb), halo_density(Rb, zb))

    # --- Robin et al. 2003 (the older runs) ---
    # thin disc: eight bins; the hole's vertical term divided by the hole's scale length
    Rdd, Rhh = 2.17, 1.33
    thin_r = np.zeros_like(x)
    for ii in range(8):
        rdi = Rb**2 + zb**2 / EPCI_ROBIN03[ii]**2
        if ii == 0:
            rho = np.exp(-rdi / 25.0) - np.exp(-rdi / 9.0)
        else:
            rho = (np.exp(-np.sqrt(0.25 + rdi / Rdd**2))
                   - np.exp(-np.sqrt(0.25 + rdi / Rhh**2)))
        rho *= 1.2
        thin_r = np.abs(thin_r + RHO0_ROBIN03[ii] * CORR_ROBIN03[ii] * 0.001 * rho / D0_ROBIN03[ii])

    # thick disc: parabolic near the plane, exponential beyond 400 pc
    nnf = 0.4 / 0.8
    rho00 = 1.34e-3 + 3.04e-4
    near = np.abs((rho00 / 0.999719) * np.exp(-(Rb - DSUN) / 2.5)
                  * (1.0 - zb**2 / (0.4 * 0.8 * (2.0 + nnf))))
    far = np.abs((rho00 / 0.999719) * np.exp(-(Rb - DSUN) / 2.5)
                 * np.exp(nnf) * np.exp(-np.abs(zb) / 0.8) / (1.0 + 0.5 * nnf))
    thick_r = np.where(np.abs(zb) < 0.4, near, far) * 2.67   # total mass 0.8e10 Msun

    # stellar halo: flattened power law, cored inside 0.5 kpc
    rdi = np.sqrt(Rb**2 + zb**2 / 0.76**2)
    halo_r = np.abs((0.932e-5 / 867.067) * np.power(np.maximum(rdi, 0.5) / DSUN, -2.44)) * 5281.0
    r2003 = (thin_r, thick_r, halo_r)

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
    mS = 35.45 / 3.84723
    rhoS = mS / np.cosh(-r4)**2 * np.where(r2 <= 3.43, 1.0, np.exp(-4.0 * (r2 - 3.43)**2))

    r4, r2 = _bar_radii(4.44, 1.31, 0.80, 2.786, 3.917)
    mE = 2.27 / 87.0
    rhoE = mE * np.exp(-r4) * np.where(r2 <= 6.83, 1.0, np.exp(-4.0 * (r2 - 6.83)**2))

    bar = np.abs(rhoS) + np.abs(rhoE)
    return lon, lat, x, {"besancon1612": b1612, "robin2003": r2003}, bar


def lens_distance_norm(prof, Ds):
    """Normaliser Z(Ds) of `func_lens`'s lens-distance sampler, for each source distance in Ds.

    `func_lens` accepts a lens at Dl with probability proportional to
    `rho(Dl) * sqrt((Ds - Dl) Dl / Ds)` over the grid points k = 1 .. nums-2. Z is the sum of
    that, so dividing the physical rate by the sampling density leaves Z as a factor. Units are Msun/pc^3 * kpc^(3/2), which cancels in any weighted fraction.
    """
    Ds = np.atleast_1d(np.asarray(Ds, dtype=float))
    k = np.arange(1, NUM)[None, :]
    Dl = k * STEP
    out = np.empty(len(Ds))
    # In row blocks: the (rows x grid) matrix and its temporaries are ~0.4 MB per row, so a
    # sightline with thousands of draws in one pass needs GBs. Blocking leaves each row's sum unchanged.
    for i in range(0, len(Ds), BLOCK):
        d = Ds[i:i + BLOCK, None]
        nums = np.round(d / STEP).astype(int)
        inside = k <= nums - 2
        tt = np.sqrt(np.clip((d - Dl) * Dl / d, 0.0, None)) * prof.rho_tot[None, :] * inside
        out[i:i + BLOCK] = tt.sum(axis=1) * STEP
    return out


BLOCK = 512


def _map_worst(sightlines):
    """Worst |dex| of the recomputed log10 Rostart and log10 Nstart against the map file's, as two
    {galaxy: dex} dicts, and the number of map rows compared."""
    wr, wn, n = {g: 0.0 for g in GALAXIES}, {g: 0.0 for g in GALAXIES}, 0
    for lon, lat, lr, ln in zip(sightlines["lon"], sightlines["lat"],
                                sightlines["log10_Rostart"], sightlines["log10_Nstart"]):
        if not np.isfinite(lon):
            continue
        # the expensive part, once per Sun height and shared by the galaxies that use it
        comps = {z: _components(lon, lat, z) for z in {spec[3] for spec in GALAXIES.values()}}
        for g, spec in GALAXIES.items():
            prof = _profile(comps[spec[3]], g)
            wr[g] = max(wr[g], abs(np.log10(prof.Rostart) - lr))
            wn[g] = max(wn[g], abs(np.log10(prof.Nstart) - ln))
        n += 1
    return wr, wn, n


# Detection compares the MASS column only: Nstart also depends on the star lists' mean masses, which
# older runs had different (their Nstart is off by up to 0.11 dex for that reason alone). The map's one
# decimal puts the right model within 0.05 dex; a wrong one is further off over a run's sightlines.
DETECT_TOL = 0.06


def galaxy_for_map(sightlines, tol=DETECT_TOL):
    """The galaxy a run used, from its map file (load_sightlines): the one whose mass columns match.

    Raises if none matches (the port has drifted from the C++) or if several do (a map with too few
    sightlines to tell them apart; pass galaxy explicitly then).
    """
    wr, _, n = _map_worst(sightlines)
    return _pick_galaxy(wr, n, tol)


def _pick_galaxy(worst, n, tol):
    fits = [g for g in GALAXIES if worst[g] <= tol]
    shown = ", ".join(f"{g}: {worst[g]:.3f} dex" for g in GALAXIES)
    if not n:
        raise ValueError("the map file has no sightline to read the galaxy model from; pass galaxy")
    if not fits:
        raise AssertionError(
            f"galaxy_model.py's mass column disagrees with the map file over {n} sightlines ({shown}; "
            f"tolerance {tol}). The port has drifted from Disk_model() in src/galaxy/density.cpp.")
    if len(fits) > 1:
        raise ValueError(
            f"the map file's {n} sightline(s) fit several galaxy models ({shown}), so it cannot say which "
            f"the run used. Pass galaxy: 'current' if its run_provenance.txt has a disc_norm line, "
            f"'bar-once' if it has a bulge_norm line only, else 'legacy'.")
    return fits[0]


def check_against_map(sightlines, tol=0.05, galaxy=None):
    """Recompute every sightline's column densities and compare with the map file's own.

    The map file writes `log10 Rostart` and `log10 Nstart` to ONE decimal, so agreement is only
    ever testable at the 0.05 dex level -- that is the guard's resolution, not its precision.
    galaxy None: the model galaxy_for_map() finds. Returns (worst_dex, n_checked, galaxy).
    Raises if the port has drifted from the C++.
    """
    wr, wn, n = _map_worst(sightlines)
    galaxy = _pick_galaxy(wr, n, DETECT_TOL) if galaxy is None else check_galaxy(galaxy)
    w = max(wr[galaxy], wn[galaxy])
    if n and w > tol:
        raise AssertionError(
            f"galaxy_model.py (galaxy {galaxy!r}) disagrees with the map file by {w:.3f} dex over {n} "
            f"sightlines (tolerance {tol}). The port has drifted from Disk_model() in src/galaxy/density.cpp.")
    return w, n, galaxy
