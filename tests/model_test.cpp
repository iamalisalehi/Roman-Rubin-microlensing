// Unit test of the light-curve model's sky geometry: the adapter to VBMicrolensing (src/events/vbm_model.cpp)
// and the Galactic axes of the kinematics (src/galaxy/kinematics.cpp), against references that use neither:
// tests/model_reference.h, written by tests/model_reference.py from astropy. Needs only committed files (the
// vendored tables and the Roman ephemeris); run from the repo root. Exit status = the number of failed
// checks (0 = all held).
//
// Pinned: (0) the setup fails with a message, not a crash, outside the repository root and then works;
// (a) the ICRS direction of Galactic (l, b); (b) the rotation from (e_l, e_b) to (North, East); (c) that the
// kinematics' projection axes (n1, n2) ARE (e_l, e_b), by projecting the Sun's velocity through the real
// vrel; (d) the Earth's velocity on (North, East) at t0; (e) that Roman's ephemeris is loaded and moves the
// track by the parallax of an observer at L2, in the right direction.
//
// Bounds on a measured difference are ~3x its value (the value is printed); the ones that compare two
// computations of the same rotation or projection are held to rounding.
//
// (a) and (b) are held to two references. The adapter uses the Hipparcos rotation, defined in the ICRS
// (ESA 1997); astropy's Galactic frame is defined through FK5 J2000 and the frame bias, and differs from it
// by ~1e-7 rad (0.02 arcsec). Against the Hipparcos rotation rebuilt independently from its defining angles
// the bound is rounding; against astropy it is that convention difference.
#include "common.h"
#include "types.h"
#include "events/vbm_model.h"
#include "galaxy/kinematics.h"
#include "model_reference.h"
#include <unistd.h>

static int nfail = 0;
// One check: the measured value, the bound it is held to, and ok/FAIL.
static void expect(bool ok, const char* what, double got, const char* rel, double want) {
    std::printf("%-76s got %10.3e  want %s %8.1e  %s\n", what, got, rel, want, ok ? "ok" : "FAIL");
    if (!ok) ++nfail;
}
static void bound(const char* what, double got, double limit) {
    expect(std::isfinite(got) && got < limit, what, got, "<", limit);
}

// Angle between two directions given as (RA, Dec) [rad], from the chord (accurate to rounding at small angles).
static double separation(double ra1, double dec1, double ra2, double dec2) {
    const double d[3] = {std::cos(dec1) * std::cos(ra1) - std::cos(dec2) * std::cos(ra2),
                         std::cos(dec1) * std::sin(ra1) - std::cos(dec2) * std::sin(ra2),
                         std::sin(dec1) - std::sin(dec2)};
    return 2.0 * std::asin(0.5 * std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
}

int main() {
    using namespace modelref;

    // (0) setup. Outside the repository root the data files are not there: vbmInit must throw, naming one,
    // and leave itself retryable.
    char here[4096];
    if (!getcwd(here, sizeof here)) { std::printf("cannot read the working directory\n"); return 1; }
    bool threw = false;
    if (chdir("/") == 0) {
        try { vbmInit(); } catch (const std::runtime_error& e) { threw = true; std::printf("  (from /: %s)\n", e.what()); }
        if (chdir(here) != 0) { std::printf("cannot return to %s\n", here); return 1; }
    }
    expect(threw, "vbmInit throws outside the repository root", threw, "=", 1.0);
    vbmInit();
    vbmInit();                                  // idempotent
    double v0[2];
    bool refused = false;
    try { earthVelocityNE(0.0, v0); } catch (const std::runtime_error&) { refused = true; }
    expect(refused, "earthVelocityNE refuses to run before a sightline is set", refused, "=", 1.0);

    SkyFrame fr[N_SIGHT];
    for (int i = 0; i < N_SIGHT; ++i) fr[i] = vbmSetSightline(GAL_L[i], GAL_B[i]);   // leaves the last set

    // (a) the ICRS direction.
    double wa = 0.0, wh = 0.0;
    for (int i = 0; i < N_SIGHT; ++i) {
        wa = std::fmax(wa, separation(fr[i].raDeg / RAa, fr[i].decDeg / RAa, RA_ASTROPY[i], DEC_ASTROPY[i]));
        wh = std::fmax(wh, separation(fr[i].raDeg / RAa, fr[i].decDeg / RAa, RA_HIPPARCOS[i], DEC_HIPPARCOS[i]));
    }
    bound("(a) RA, Dec of (l, b) vs the Hipparcos rotation (defining angles): max separation [rad]", wh, 1e-12);
    bound("(a) RA, Dec of (l, b) vs astropy's Galactic frame: max separation [rad]", wa, 3e-7);

    // (b) (e_l, e_b) -> (North, East).
    double ta = 0.0, th = 0.0, orth = 0.0, det = 0.0;
    for (int i = 0; i < N_SIGHT; ++i) {
        for (int r = 0; r < 2; ++r)
            for (int c = 0; c < 2; ++c) {
                ta = std::fmax(ta, std::fabs(fr[i].toNE[r][c] - TONE_ASTROPY[i][r][c]));
                th = std::fmax(th, std::fabs(fr[i].toNE[r][c] - TONE_HIPPARCOS[i][r][c]));
            }
        const double (*m)[2] = fr[i].toNE;
        orth = std::fmax(orth, std::fabs(m[0][0] * m[0][0] + m[0][1] * m[0][1] - 1.0));
        orth = std::fmax(orth, std::fabs(m[0][0] * m[1][0] + m[0][1] * m[1][1]));
        det = std::fmax(det, std::fabs(m[0][0] * m[1][1] - m[0][1] * m[1][0] + 1.0));
    }
    bound("(b) toNE vs the Hipparcos rotation (defining angles): max |entry difference|", th, 1e-12);
    bound("(b) toNE vs astropy's Galactic frame: max |entry difference|", ta, 2e-7);
    bound("(b) toNE orthogonal: max |M M^T - 1|", orth, 1e-14);
    bound("(b) toNE has determinant -1 (east-then-north vs north-then-east): max |det + 1|", det, 1e-14);

    // (c) The projection axes of the kinematics. vrel projects every velocity on the lens plane as
    // (n1, n2), and the Sun's own velocity with the formulas for the line of sight at TET = 360 - l, FI = b
    // (src/sim/sightline.cpp). If (n1, n2) are (e_l, e_b), the Sun's (n1, n2) are the components of its
    // Galactic velocity (U, V, W) = (VSunR, VSunT, VSunZ) on (e_l, e_b), with U toward the Galactic centre,
    // V toward Galactic rotation (l = 90) and W toward the north Galactic pole. vrel's draws are irrelevant
    // here: the Sun's projection does not depend on them.
    const double sun[3] = {VSunR, VSunT, VSunZ};
    double worstSun = 0.0;
    std::printf("  Sun's velocity (U, V, W) = (%.2f, %.2f, %.2f) km/s through vrel, vs astropy on (e_l, e_b):\n", sun[0], sun[1], sun[2]);
    std::printf("  %6s %6s %12s %12s %12s %12s\n", "l", "b", "n1", "e_l.v", "n2", "e_b.v");
    source s;
    lens l;
    for (int i = 0; i < N_SIGHT; ++i) {
        s.lon = GAL_L[i];
        s.lat = GAL_B[i];
        s.TET = (360.0 - s.lon) / RAa;          // as src/sim/sightline.cpp
        s.FI  = s.lat / RAa;
        s.struc = l.struc = GalacticComponent::BULGE;
        s.Ds = 8.0;
        l.Dl = 4.0;
        l.xls = l.Dl / s.Ds;
        vrel(s, l);
        double w1 = 0.0, w2 = 0.0;
        for (int j = 0; j < 3; ++j) { w1 += SUN_PROJ[i][0][j] * sun[j]; w2 += SUN_PROJ[i][1][j] * sun[j]; }
        std::printf("  %6.1f %6.1f %12.6f %12.6f %12.6f %12.6f\n", GAL_L[i], GAL_B[i], s.VSun_n1, w1, s.VSun_n2, w2);
        worstSun = std::fmax(worstSun, std::fmax(std::fabs(s.VSun_n1 - w1), std::fabs(s.VSun_n2 - w2)));
    }
    bound("(c) vrel's (n1, n2) of the Sun vs its (U, V, W) on (e_l, e_b): max |difference| [km/s]", worstSun, 1e-12);

    // (d) the Earth's velocity on (North, East) at t0, the v_E,perp of mu_geo = mu_hel - v_E,perp pi_rel. The
    // library interpolates a daily table (the second difference of a 1-day step on the orbit is ~5e-5 of the
    // velocity), so the difference is held to that, as a fraction of the Earth's speed.
    double worstV = 0.0, speed = 0.0;
    for (int i = 0; i < N_SIGHT; ++i) {
        vbmSetSightline(GAL_L[i], GAL_B[i]);
        for (int k = 0; k < N_T; ++k) {
            double v[2];
            earthVelocityNE(T_DAY[k], v);
            const double vs = std::hypot(EARTH_V_NE[i][k][0], EARTH_V_NE[i][k][1]);
            speed = std::fmax(speed, vs);
            worstV = std::fmax(worstV, std::hypot(v[0] - EARTH_V_NE[i][k][0], v[1] - EARTH_V_NE[i][k][1]) / vs);
        }
    }
    std::printf("  (largest projected speed in the reference: %.5f AU/day)\n", speed);
    bound("(d) Earth velocity on (N, E) at t0 vs astropy: max |dv| / |v|", worstV, 1.3e-4);

    // (e) Roman's ephemeris is loaded, selected by satellite = 1, and displaces the track as an observer at
    // L2 (L2_KM / AU_KM away from the Earth, away from the Sun) does.
    double worstL2 = 0.0, refMax = 0.0;
    int flags = 0;
    for (int i = 0; i < N_SIGHT; ++i) {
        vbmSetSightline(GAL_L[i], GAL_B[i]);
        for (int k = 0; k < N_T; ++k)
            for (int ip = 0; ip < N_PI; ++ip) {
                double dy[2];
                flags += vbmTestSatelliteShift(T_DAY[k], PIE[ip][0], PIE[ip][1], dy);
                const double* ref = L2_SHIFT[i][k][ip];
                const double rn = std::hypot(ref[0], ref[1]);
                refMax = std::fmax(refMax, rn);
                worstL2 = std::fmax(worstL2, std::hypot(dy[0] - ref[0], dy[1] - ref[1]) / rn);
            }
    }
    std::printf("  (largest shift in the reference: %.3e Einstein radii)\n", refMax);
    bound("(e) track shift of Roman vs the Earth against the L2 parallax: max |dy - ref| / |ref|", worstL2, 1e-4);
    expect(flags == 0, "(e) no extrapolation flag raised for any in-range time", flags, "=", 0.0);
    double out[2];
    const int flagOut = vbmTestSatelliteShift(-500.0, 0.3, -0.2, out);   // before Roman's table starts
    expect(flagOut == 2, "(e) a time before Roman's table raises the satellite extrapolation flag (2)", flagOut, "=", 2.0);

    std::printf("%s (%d failed)\n", nfail ? "FAILED" : "all held", nfail);
    return nfail;
}
