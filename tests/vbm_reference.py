"""Reference values for tests/vbm_test.cpp, the unit test of the vendored VBMicrolensing library.

    .roman/bin/python tests/vbm_reference.py -o tests/vbm_reference.h

writes the C++ header the test includes (the public .gitignore keeps data files out of the repository, so
the reference is a header). Deterministic, about a minute. Nothing here calls VBMicrolensing. Re-run it only
when a reference itself changes (grid, events, sightlines): updating VBMicrolensing must not need it,
since the test exists to show that the library's answers have not moved.

(a) Finite-source point-lens magnification, by direct integration over the source disc. In polar
    coordinates (d, phi) centred on the lens the point-source magnification depends only on d, and
    A_ps(d) d = (d^2 + 2) / sqrt(d^2 + 4) stays finite at d -> 0, so for a disc of radius rho at u
        A = 1 / (pi rho^2 Ibar) int dd A_ps(d) d int_{phi in disc} I(r) dphi,
    r^2 = d^2 + u^2 - 2 d u cos(phi), with linear limb darkening I = 1 - a1 (1 - sqrt(1 - r^2/rho^2)) and
    Ibar = 1 - a1/3 its mean over the disc. Both integrals by adaptive quadrature at 1e-11 relative;
    loosening that to 1e-9 moves no value by more than 3e-10. Grid: rho x z = u/rho x a1, and separately
    the source centre z = 0 and 1e-3, which the test prints without asserting.

(b) needs no reference: the test holds PSPLMag to our own point-source formula.

(c) Annual-parallax source track, from astropy's built-in ephemeris (ERFA epv00: Earth's heliocentric
    position within 11 km of JPL DE405 over 1900-2100), in the convention VBMicrolensing implements
    (ComputeParallax, PSPLAstroLightCurve) for PSPLLightCurveParallax with t_in_HJD = false and t0_par = t0:
      n            unit vector to the target; S, W the celestial South and West unit vectors on the sky there;
      E(t)         Earth's heliocentric position [AU] projected on (S, W), time t = JD - 2450000 read as TDB;
      lt(t)        = (r_E(t) . n) AU/c, the light-travel time between the Earth and the plane through the
                   Sun normal to n; t' = t + lt(t) is heliocentric time;
      Et(t)        = E(t) - E(t0) - V0 (t' - t0'): the Earth's offset from its straight-line motion at t0,
                   V0 = dE/dt at t0 by central difference with h = 0.5 d (the geocentric frame);
      tau          = (t' - t0') / tE + piE_N Et_S + piE_E Et_W,   beta = u0 + piE_N Et_W - piE_E Et_S,
      (y1, y2)     = (-tau, -beta), the source position relative to the lens in Einstein radii.
    The light-travel term reaches ~1.8e-4 Einstein radii here, so the test's 1e-4 bound resolves it.
"""
import argparse
import sys

import numpy as np
from scipy.integrate import quad
import astropy.units as u
from astropy.coordinates import SkyCoord, get_body_barycentric, solar_system_ephemeris
from astropy.time import Time

# ---- (a) finite source ------------------------------------------------------------------------------
RHO = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
Z = [0.01, 0.1, 0.5, 0.9, 0.99, 1.01, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0]   # z = u / rho
Z_CENTRE = [0.0, 1e-3]
A1 = [0.0, 0.5]
EPS = 1e-11


def espl(u_, rho, a1):
    """Magnification of a limb-darkened disc of radius rho at distance u_ from a point lens."""
    Ibar = 1.0 - a1 / 3.0

    def inner(d):                      # int over the arc of the circle of radius d inside the disc of I
        if u_ == 0.0:
            return 2.0 * np.pi * (1.0 - a1 * (1.0 - np.sqrt(1.0 - d * d / rho ** 2))) if d < rho else 0.0
        c = (d * d + u_ * u_ - rho * rho) / (2.0 * d * u_) if d > 0 else -2.0
        if c >= 1.0:
            return 0.0
        pmax = np.pi if c <= -1.0 else np.arccos(c)
        if a1 == 0.0:
            return 2.0 * pmax
        I = lambda p: 1.0 - a1 * (1.0 - np.sqrt(max(0.0, 1.0 - (d * d + u_ * u_ - 2 * d * u_ * np.cos(p)) / rho ** 2)))
        return 2.0 * quad(I, 0.0, pmax, epsabs=EPS * 1e-2, epsrel=EPS, limit=200)[0]

    lo, hi = max(0.0, u_ - rho), u_ + rho
    kink = [x for x in (rho - u_,) if lo < x < hi]        # where the circle of radius d first leaves the disc
    f = lambda d: (d * d + 2.0) / np.sqrt(d * d + 4.0) * inner(d)
    val = quad(f, lo, hi, points=kink or None, epsabs=EPS * 1e-2, epsrel=EPS, limit=500)[0]
    return val / (np.pi * rho * rho * Ibar)


# ---- (c) parallax -----------------------------------------------------------------------------------
SIGHTLINES = [(1.0, -1.5), (10.0, -3.0)]      # Galactic (l, b) [deg]: Roman's GBTDS fields, and one 9 deg away
# Simulation day 0 = MJD 61141.312002288 (TIME0_MJD in Baseline/readbaselineBulge.py), as JD - 2450000.
SIM_DAY0 = 61141.312002288 - 49999.5
T0_DAYS = [400.0, 1500.0]                       # t0 on the simulation clock: 2027-05, 2030-05
U0, TE, PIEN, PIEE = 0.15, 60.0, 0.3, -0.2
OFFSETS = np.arange(-300.0, 300.0 + 1e-9, 5.0)  # t - t0 [d]
AU_C = 0.005775518331436995                     # AU / c [d]


def coords(l, b):
    """'hh:mm:ss.ss -dd:mm:ss.ss' (ICRS) -- the string VBMicrolensing's SetObjectCoordinates parses."""
    return SkyCoord(l=l * u.deg, b=b * u.deg, frame='galactic').icrs.to_string('hmsdms', sep=':', precision=2)


def track(coord, t0, t):
    n = SkyCoord(coord, unit=(u.hourangle, u.deg)).cartesian.xyz.value   # the rounded string both sides use
    zax = np.array([0.0, 0.0, 1.0])
    north = zax - (zax @ n) * n
    north /= np.linalg.norm(north)
    east = np.cross(zax, n)
    east /= np.linalg.norm(east)
    S, W = -north, -east

    def earth(tt):                    # projections on (S, W) [AU] and the light-travel time [d]
        tt = np.atleast_1d(tt)
        T = Time(np.full(tt.shape, 2450000.0), tt, format='jd', scale='tdb')
        r = (get_body_barycentric('earth', T) - get_body_barycentric('sun', T)).xyz.to(u.AU).value.T
        return np.stack([r @ S, r @ W], 1), (r @ n) * AU_C

    E, lt = earth(t)
    E0, lt0 = earth(t0)
    h = 0.5
    V0 = (earth(t0 + h)[0] - earth(t0 - h)[0]) / (2 * h)
    dtp = (t + lt) - (t0 + lt0)
    Et = E - E0 - V0 * dtp[:, None]
    par_tau = PIEN * Et[:, 0] + PIEE * Et[:, 1]
    par_beta = PIEN * Et[:, 1] - PIEE * Et[:, 0]
    y1 = -(dtp / TE + par_tau)
    y2 = -(U0 + par_beta)
    return y1, y2, np.hypot(par_tau, par_beta).max()


# ---- header -----------------------------------------------------------------------------------------
def arr(vals, per_line=4, indent='    '):
    """A C++ brace initialiser holding every double exactly (repr round-trips); nested per dimension."""
    vals = np.asarray(vals, dtype=float)
    if vals.ndim > 1:
        return '{\n' + ',\n'.join(indent + arr(v, per_line, indent + '    ') for v in vals) + '}'
    v = [repr(float(x)) for x in vals]
    lines = [', '.join(v[i:i + per_line]) for i in range(0, len(v), per_line)]
    return '{' + (',\n' + indent).join(lines) + '}'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('-o', '--output', default='-', help="header to write ('-' = stdout)")
    args = ap.parse_args()

    espl_grid = np.array([[[espl(z * rho, rho, a1) for z in Z] for rho in RHO] for a1 in A1])
    espl_centre = np.array([[[espl(z * rho, rho, a1) for z in Z_CENTRE] for rho in RHO] for a1 in A1])

    with solar_system_ephemeris.set('builtin'):
        cstr = [coords(l, b) for l, b in SIGHTLINES]
        t0s = [SIM_DAY0 + d for d in T0_DAYS]
        times = np.array([t0 + OFFSETS for t0 in t0s])
        Y1, Y2, PAR = [], [], []
        for c in cstr:
            r = [track(c, t0, times[k]) for k, t0 in enumerate(t0s)]
            Y1.append([x[0] for x in r]); Y2.append([x[1] for x in r]); PAR.append([x[2] for x in r])

    out = ['// GENERATED by tests/vbm_reference.py (.roman/bin/python tests/vbm_reference.py -o tests/vbm_reference.h).',
           '// Do not edit by hand. The references and their definitions are described in that script.',
           '#ifndef VBM_REFERENCE_H', '#define VBM_REFERENCE_H', '', 'namespace vbmref {', '',
           '// (a) Finite-source point-lens magnification by direct integration over the disc.',
           '//     ESPL[ia1][irho][iz] at u = Z[iz] * RHO[irho], limb-darkening coefficient A1[ia1];',
           '//     ESPL_CENTRE[ia1][irho][izc] the same at u = Z_CENTRE[izc] * RHO[irho], near the disc centre.',
           'constexpr int N_RHO = %d, N_Z = %d, N_A1 = %d, N_ZC = %d;' % (len(RHO), len(Z), len(A1), len(Z_CENTRE)),
           'constexpr double RHO[N_RHO] = %s;' % arr(RHO, 8),
           'constexpr double Z[N_Z] = %s;' % arr(Z, 8),
           'constexpr double A1[N_A1] = %s;' % arr(A1, 8),
           'constexpr double Z_CENTRE[N_ZC] = %s;' % arr(Z_CENTRE, 8),
           'constexpr double ESPL[N_A1][N_RHO][N_Z] = %s;' % arr(espl_grid),
           'constexpr double ESPL_CENTRE[N_A1][N_RHO][N_ZC] = %s;' % arr(espl_centre), '',
           '// (c) Annual-parallax source track (y1, y2) [Einstein radii] from astropy\'s ephemeris, geocentric',
           '//     frame at t0. Times are JD - 2450000. COORDS[is] is ICRS for Galactic (GAL_L, GAL_B)[is];',
           '//     T0[ie] is simulation day T0_DAY[ie]; PAR_MAX[is][ie] is the largest |parallax term|.',
           'constexpr int N_SIGHT = %d, N_EVENT = %d, N_T = %d;' % (len(cstr), len(t0s), len(OFFSETS)),
           'constexpr const char* COORDS[N_SIGHT] = {%s};' % ', '.join('"%s"' % c for c in cstr),
           'constexpr double GAL_L[N_SIGHT] = %s;' % arr([l for l, _ in SIGHTLINES], 8),
           'constexpr double GAL_B[N_SIGHT] = %s;' % arr([b for _, b in SIGHTLINES], 8),
           'constexpr double T0_DAY[N_EVENT] = %s;' % arr(T0_DAYS, 8),
           'constexpr double T0[N_EVENT] = %s;' % arr(t0s, 8),
           'constexpr double U0 = %r, TE = %r, PIEN = %r, PIEE = %r;' % (U0, TE, PIEN, PIEE),
           'constexpr double T[N_EVENT][N_T] = %s;' % arr(times),
           'constexpr double Y1[N_SIGHT][N_EVENT][N_T] = %s;' % arr(Y1),
           'constexpr double Y2[N_SIGHT][N_EVENT][N_T] = %s;' % arr(Y2),
           'constexpr double PAR_MAX[N_SIGHT][N_EVENT] = %s;' % arr(PAR), '',
           '}  // namespace vbmref', '', '#endif // VBM_REFERENCE_H', '']
    text = '\n'.join(out)
    if args.output == '-':
        sys.stdout.write(text)
    else:
        with open(args.output, 'w') as f:
            f.write(text)


if __name__ == '__main__':
    main()
