"""Reference values for tests/model_test.cpp, the unit test of the light-curve model's sky geometry
(src/events/vbm_model.cpp, the adapter to VBMicrolensing).

    .roman/bin/python tests/model_reference.py -o tests/model_reference.h

writes the C++ header the test includes (the public .gitignore keeps data files out of the repository, so
the reference is a header). Deterministic, a few seconds. Nothing here calls the C++ code or
VBMicrolensing. Re-run it only when a reference itself changes (sightlines, times); the constants it
reads from the C++ (the simulation clock, L2) come through analysis/cparams.py.

Five references, at three bulge sightlines (Galactic (l, b) in degrees):

(a) The ICRS direction of (l, b), twice. ASTROPY: astropy's Galactic frame, which defines the Galactic
    system through FK5 J2000 and the FK5-ICRS frame bias. HIPPARCOS: the rotation the Hipparcos catalogue
    defines directly in the ICRS (ESA 1997 vol. 1 sec. 1.5.3), rebuilt here from its three defining angles
    (north Galactic pole at (192.85948, +27.12825) deg, Galactic longitude of the north celestial pole
    122.93192 deg) and not from the nine matrix elements the C++ carries. The two differ by ~1e-7 rad
    (0.02 arcsec). The C++ uses the Hipparcos rotation, so it is held to HIPPARCOS at rounding level and to
    ASTROPY at the size of that convention difference.
(b) The 2x2 rotation from (e_l, e_b) to (North, East), from the same two sources. ASTROPY: the proper motion
    of a source moving 1 mas/yr toward increasing l (then b) at 1 kpc, transformed to the ICRS; its
    (pm_dec, pm_ra_cosdec) are the (North, East) components of e_l (then e_b).
(c) Where a Galactic Cartesian velocity (x toward the Galactic centre, y toward l = 90, z toward the north
    Galactic pole) goes on (e_l, e_b): SUN_PROJ[sightline][e_l or e_b][x, y, z], astropy's own conversion of
    a unit velocity at 1 kpc to (pm_l_cosb, pm_b). The test multiplies it by the Sun's velocity
    (VSunR, VSunT, VSunZ), so a change to those constants does not stale the header.
(d) The Earth's heliocentric velocity projected on (North, East), AU/day, from astropy's ephemeris (ERFA
    epv00), at three simulation days. The time is read as UTC (the Sun table is JPL Horizons' JDUT), so the
    comparison measures the library's daily-table interpolation and not a 69 s clock offset.
(e) How far a source track moves, in Einstein radii, when the observer is at L2 (L2_KM / AU_KM along the
    Sun-to-Earth direction from the Earth) instead of the Earth, for two parallax vectors: with the
    observer displaced by s from the Earth, the library's source-lens track changes by
        dy1 = piE_N s_N + piE_E s_E,      dy2 = piE_N s_E - piE_E s_N
    (y1 along the motion, y2 across, s the displacement's North and East components in AU).
"""
import argparse
import os
import sys
import warnings

import numpy as np
import astropy.units as u
from astropy.coordinates import (CartesianDifferential, Galactic, SkyCoord, SphericalCosLatDifferential,
                                 SphericalRepresentation, get_body_barycentric, get_body_barycentric_posvel,
                                 solar_system_ephemeris)
from astropy.time import Time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analysis'))
from cparams import P  # noqa: E402

SIGHTLINES = [(1.0, -1.5), (10.0, -3.0), (358.5, -4.2)]   # Galactic (l, b) [deg]: Roman's fields, and two more
T_DAYS = [400.0, 1500.0, 3000.37]                           # simulation days (the last is off the daily grid)
PIE = [(0.3, -0.2), (-0.1, 0.5)]                            # (piE_N, piE_E)
L2_AU = P.L2_KM / P.AU_KM
DAY0_JD = P.SIM_DAY0_JD2450000 + 2450000.0                  # JD of simulation day 0


def north_east(ra, dec):
    """Celestial North and East unit vectors (ICRS Cartesian) at (ra, dec) [rad]."""
    return (np.array([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)]),
            np.array([-np.sin(ra), np.cos(ra), 0.0]))


def hipparcos_axes():
    """Rows: the Galactic x, y, z axes in the ICRS, from the Hipparcos defining angles."""
    ag, dg, th = np.radians([192.85948, 27.12825, 122.93192])
    gz = np.array([np.cos(dg) * np.cos(ag), np.cos(dg) * np.sin(ag), np.sin(dg)])
    toward_ncp = (np.array([0.0, 0.0, 1.0]) - np.sin(dg) * gz) / np.cos(dg)   # in the plane, at longitude th
    w = np.cross(gz, toward_ncp)
    return np.stack([np.cos(th) * toward_ncp - np.sin(th) * w, np.sin(th) * toward_ncp + np.cos(th) * w, gz])


def gal_unit_vectors(l, b):
    """Galactic Cartesian unit vectors toward (l, b), and toward increasing l and increasing b there."""
    l, b = np.radians(l), np.radians(b)
    return (np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)]),
            np.array([-np.sin(l), np.cos(l), 0.0]),
            np.array([-np.sin(b) * np.cos(l), -np.sin(b) * np.sin(l), np.cos(b)]))


def hipparcos_frame(l, b):
    """(ra, dec) [rad] and the (e_l, e_b) -> (North, East) matrix from the Hipparcos rotation."""
    A = hipparcos_axes()
    n, el, eb = (A.T @ v for v in gal_unit_vectors(l, b))
    ra, dec = np.arctan2(n[1], n[0]) % (2 * np.pi), np.arctan2(n[2], np.hypot(n[0], n[1]))
    eN, eE = north_east(ra, dec)
    return ra, dec, np.array([[eN @ el, eN @ eb], [eE @ el, eE @ eb]])


def astropy_frame(l, b):
    """(ra, dec) [rad] and the (e_l, e_b) -> (North, East) matrix from astropy's Galactic frame."""
    c = SkyCoord(l=l * u.deg, b=b * u.deg, frame='galactic').icrs
    M = np.zeros((2, 2))
    for j, (pml, pmb) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        p = SkyCoord(l=l * u.deg, b=b * u.deg, distance=1 * u.kpc, pm_l_cosb=pml * u.mas / u.yr,
                     pm_b=pmb * u.mas / u.yr, radial_velocity=0 * u.km / u.s, frame='galactic').icrs
        M[0, j], M[1, j] = p.pm_dec.value, p.pm_ra_cosdec.value
    return c.ra.rad, c.dec.rad, M


def sun_projection(l, b):
    """[e_l, e_b][x, y, z]: the components on (e_l, e_b) of a unit Galactic Cartesian velocity at (l, b)."""
    M = np.zeros((2, 3))
    pos = SphericalRepresentation(l * u.deg, b * u.deg, 1 * u.kpc).to_cartesian()
    for j in range(3):
        vel = np.zeros(3)
        vel[j] = 1.0
        g = Galactic(pos.with_differentials(CartesianDifferential(vel * u.km / u.s)))
        s = g.represent_as(SphericalRepresentation, SphericalCosLatDifferential)
        d = s.differentials['s']
        M[0, j] = (d.d_lon_coslat * s.distance).to(u.km / u.s, u.dimensionless_angles()).value
        M[1, j] = (d.d_lat * s.distance).to(u.km / u.s, u.dimensionless_angles()).value
    return M


def earth_state(day):
    """Earth's heliocentric velocity [AU/day] and unit position vector (ICRS) at simulation day `day`."""
    t = Time(DAY0_JD + day, format='jd', scale='utc')
    pe, ve = get_body_barycentric_posvel('earth', t)
    ps, vs = get_body_barycentric_posvel('sun', t)
    r = (pe - ps).xyz.to(u.AU).value
    return (ve - vs).xyz.to(u.AU / u.day).value, r / np.linalg.norm(r)


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

    warnings.simplefilter('ignore')        # ERFA 'dubious year' for UTC dates beyond its leap-second table
    with solar_system_ephemeris.set('builtin'):
        ast = [astropy_frame(l, b) for l, b in SIGHTLINES]
        hip = [hipparcos_frame(l, b) for l, b in SIGHTLINES]
        sunp = [sun_projection(l, b) for l, b in SIGHTLINES]
        earth = [earth_state(d) for d in T_DAYS]
    # (d) and (e) at astropy's pointing; the C++ points at the Hipparcos one, 1e-7 rad away, which moves
    # the projected velocity by that fraction.
    vne, shift = [], []
    for ra, dec, _ in ast:
        eN, eE = north_east(ra, dec)
        vne.append([[v @ eN, v @ eE] for v, _ in earth])
        sN_sE = [[L2_AU * r @ eN, L2_AU * r @ eE] for _, r in earth]
        shift.append([[[pn * sn + pe * se, pn * se - pe * sn] for pn, pe in PIE] for sn, se in sN_sE])

    out = ['// GENERATED by tests/model_reference.py (.roman/bin/python tests/model_reference.py -o tests/model_reference.h).',
           '// Do not edit by hand. The references and their definitions are described in that script.',
           '#ifndef MODEL_REFERENCE_H', '#define MODEL_REFERENCE_H', '', 'namespace modelref {', '',
           'constexpr int N_SIGHT = %d, N_T = %d, N_PI = %d;' % (len(SIGHTLINES), len(T_DAYS), len(PIE)),
           'constexpr double GAL_L[N_SIGHT] = %s;' % arr([l for l, _ in SIGHTLINES], 8),
           'constexpr double GAL_B[N_SIGHT] = %s;' % arr([b for _, b in SIGHTLINES], 8), '',
           '// (a) ICRS direction [rad]: astropy\'s Galactic frame, and the Hipparcos rotation from its defining angles.',
           'constexpr double RA_ASTROPY[N_SIGHT] = %s;' % arr([a[0] for a in ast], 8),
           'constexpr double DEC_ASTROPY[N_SIGHT] = %s;' % arr([a[1] for a in ast], 8),
           'constexpr double RA_HIPPARCOS[N_SIGHT] = %s;' % arr([h[0] for h in hip], 8),
           'constexpr double DEC_HIPPARCOS[N_SIGHT] = %s;' % arr([h[1] for h in hip], 8), '',
           '// (b) (e_l, e_b) -> (North, East), TONE[sightline][N or E][l or b].',
           'constexpr double TONE_ASTROPY[N_SIGHT][2][2] = %s;' % arr([a[2] for a in ast]),
           'constexpr double TONE_HIPPARCOS[N_SIGHT][2][2] = %s;' % arr([h[2] for h in hip]), '',
           '// (c) A unit Galactic Cartesian velocity on (e_l, e_b): SUN_PROJ[sightline][e_l or e_b][x, y, z].',
           'constexpr double SUN_PROJ[N_SIGHT][2][3] = %s;' % arr(sunp), '',
           '// (d) Earth\'s heliocentric velocity on (North, East) [AU/day] at simulation days T_DAY.',
           'constexpr double T_DAY[N_T] = %s;' % arr(T_DAYS, 8),
           'constexpr double EARTH_V_NE[N_SIGHT][N_T][2] = %s;' % arr(vne), '',
           '// (e) The track shift [Einstein radii] of an observer at L2 relative to the Earth, for parallax vectors',
           '//     PIE[ip] = (piE_N, piE_E): L2_SHIFT[sightline][time][ip][y1, y2].',
           'constexpr double PIE[N_PI][2] = %s;' % arr(PIE),
           'constexpr double L2_SHIFT[N_SIGHT][N_T][N_PI][2] = %s;' % arr(shift), '',
           '}  // namespace modelref', '', '#endif // MODEL_REFERENCE_H', '']
    text = '\n'.join(out)
    if args.output == '-':
        sys.stdout.write(text)
    else:
        with open(args.output, 'w') as f:
            f.write(text)


if __name__ == '__main__':
    main()
