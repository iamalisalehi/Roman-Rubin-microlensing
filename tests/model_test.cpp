// Unit test of the light-curve model's sky geometry: the adapter to VBMicrolensing (src/events/vbm_model.cpp)
// and the Galactic axes of the kinematics (src/galaxy/kinematics.cpp), against references that use neither:
// tests/model_reference.h, written by tests/model_reference.py from astropy. Needs only committed files (the
// vendored tables and the Roman ephemeris); run from the repo root. Exit status = the number of failed
// checks (0 = all held).
//
// Pinned: (0) the setup fails with a message, not a crash, outside the repository root and then works;
// (a) the ICRS direction of Galactic (l, b); (b) the rotation from (e_l, e_b) to (North, East); (c) that the
// kinematics' projection axes (n1, n2) ARE (e_l, e_b), by projecting the Sun's velocity through the real
// vrel, and that the direction of motion xi vrel draws is the one geocentricMotion reproduces; (d) the Earth's
// velocity on (North, East) at t0; (e) that Roman's ephemeris is loaded and moves the track by the parallax
// of an observer at L2, in the right direction.
//
// The light-curve model proper, as the simulator calls it (vbmEvaluate, evaluateModel, evaluateTrack):
// (f) the side of the lens the source passes and the handedness of xi, against the straight line
// u0 (-sin xi, cos xi) + tau (cos xi, sin xi) on (e_l, e_b) when the parallax is switched off; (g) the
// magnification is magnifOf(u) bit for bit in every path; (h) the image centroid is thetaE u/(u^2 + 2) from
// the unlensed source, along the source-lens vector, and the lens is -thetaE u from it, also at thousands of
// epochs of one call (the library perturbs the parallax vector by 1e-10 at every epoch of an astrometric
// call; the adapter follows it); (i) the heliocentric-to-geocentric conversion of tE and xi against numpy
// with astropy's Earth velocity; (j) the Roman-Earth separation of the source-lens vectors at t0 (du_sat);
// (k) the annual parallax: the Earth's orbit, as the track traces it, in position and in shape (axis ratio
// and orientation, not the circle of the legacy model); (l) the source's unlensed track does not depend on
// the lens, and source and lens move by their heliocentric proper motions over a year, which ties the
// geocentric draw to the heliocentric velocities; (m) Roman at the centre of the Earth is the Earth, at L2 it
// is not, and nothing else is allowed; (n) the two seams agree with each other and with the straight line;
// (o) evaluateTrack with and without the absolute positions: the same u, A, ux, uy either way, the point-lens
// deflection and NaN positions without; (p) an epoch outside the ephemerides is refused.
//
// The finite source (the uniform disc, --source-model espl): (q) the adapter's magnification is the library's
// ESPLMag2 and matches the exact magnification of a uniform disc, computed here by quadrature; (r) it goes to
// the point source as the radius goes to 0; (s) at the peak of an event that crosses the disc the
// magnification is finite and below the point source's, the centroid shift smaller than the point lens's,
// and the two seams agree; (t) under --source-model pspl a source radius changes nothing, bit for bit; (u) the
// no-parallax magnification of the straight line is the same function; (v) the edges: a source of no size,
// one beyond the library's table, the source on the lens.
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
#include "events/lightcurve.h"
#include "galaxy/kinematics.h"
#include "model_reference.h"
#include "vbm_reference.h"
#include "VBMicrolensingLibrary.h"
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

// The exact magnification and flux-weighted image position of a UNIFORM DISC of radius rho [Einstein radii]
// whose centre is u from a point lens, by quadrature, with no code from the library. In polar coordinates
// (r, phi) centred on the lens, a point source at distance r has the magnification A(r) = (r^2 + 2) /
// (r sqrt(r^2 + 4)) and the flux-weighted position of its images A x_img = (r^2 + 3) / sqrt(r^2 + 4) along
// its own direction. The disc's values are the averages over it. The circle of radius r around the lens lies
// wholly inside the disc for r < rho - u (u < rho) and its arc inside the disc subtends 2 phi(r),
// cos phi = (r^2 + u^2 - rho^2) / (2 r u), for |u - rho| < r < u + rho. Then
//     A  = [ pi R1 sqrt(R1^2 + 4) + int 2 phi r A(r) dr ] / (pi rho^2),   R1 = max(rho - u, 0)
//          (r A(r) = (r^2 + 2)/sqrt(r^2 + 4) integrates to r sqrt(r^2 + 4)/2),
//     h  = int 2 sin(phi) r (r^2 + 3)/sqrt(r^2 + 4) dr / (pi rho^2)       (component along the axis lens -> disc centre),
// over the arcs. r = m - hw cos(theta) between the ends makes the integrands smooth (Simpson, 20000 steps).
// The centroid of the images lies h/A from the lens, and the deflection from the disc's centre is h/A - u.
struct DiscRef { double A, h; };
static DiscRef uniformDisc(double u, double rho) {
    const double a = std::fabs(u - rho), b = u + rho, m = 0.5 * (a + b), hw = 0.5 * (b - a);
    const double R1 = std::fmax(rho - u, 0.0);
    const int N = 20000;                                   // even
    double sumA = 0.0, sumH = 0.0;
    for (int k = 0; k <= N; ++k) {
        const double th = pi * k / N, ct = std::cos(th), st = std::sin(th);
        const double r = m - hw * ct;
        // cos phi, written without the cancellation of the plain expression
        double c;
        if (u < rho) c = (u * (1.0 + ct * ct) - 2.0 * rho * ct) / (2.0 * r);
        else         c = 1.0 - rho * rho * st * st / (2.0 * r * u);
        c = std::fmax(-1.0, std::fmin(1.0, c));
        const double phi = std::acos(c), sphi = std::sqrt(std::fmax(0.0, 1.0 - c * c));
        const double w = (k == 0 || k == N) ? 1.0 : (k % 2 ? 4.0 : 2.0);
        const double dr = hw * st;                         // dr/dtheta
        sumA += w * 2.0 * phi * (r * r + 2.0) / std::sqrt(r * r + 4.0) * dr;
        sumH += w * 2.0 * sphi * r * (r * r + 3.0) / std::sqrt(r * r + 4.0) * dr;
    }
    const double h3 = pi / N / 3.0, norm = pi * rho * rho;
    return DiscRef{(pi * R1 * std::sqrt(R1 * R1 + 4.0) + sumA * h3) / norm, sumH * h3 / norm};
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
    double worstSun = 0.0, worstXi = 0.0, worstMu = 0.0;
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

        // vrel's xi is the direction of the source-minus-lens relative velocity on (n1, n2); geocentricMotion,
        // fed the same draw's proper motions with the Earth's motion switched off (pi_rel = 0), must point
        // the same way and give the same mu_rel. Many draws, disk and bulge lenses and sources.
        vbmSetSightline(GAL_L[i], GAL_B[i]);
        for (int k = 0; k < 60; ++k) {
            s.struc = (k % 3 == 0) ? GalacticComponent::THIN_DISK : GalacticComponent::BULGE;
            l.struc = (k % 2 == 0) ? GalacticComponent::THICK_DISK : GalacticComponent::BULGE;
            s.Ds = 6.0 + 0.1 * k;
            l.Dl = 1.0 + 0.08 * k;
            l.xls = l.Dl / s.Ds;
            vrel(s, l);
            const GeocentricMotion g0 = geocentricMotion(1500.0, 1.0, 0.0, s.mus1, s.mus2, l.mul1, l.mul2);
            worstXi = std::fmax(worstXi, std::fabs(std::remainder(g0.xi - s.xi, 2.0 * pi)));
            worstMu = std::fmax(worstMu, std::fabs(g0.muGeo / l.murel - 1.0));
        }
    }
    bound("(c) vrel's (n1, n2) of the Sun vs its (U, V, W) on (e_l, e_b): max |difference| [km/s]", worstSun, 1e-12);
    bound("(c) xi of vrel vs geocentricMotion at pi_rel = 0 (the same direction of motion): max |dxi| [rad]", worstXi, 1e-12);
    bound("(c) |mu_rel| of vrel vs geocentricMotion at pi_rel = 0: max relative difference", worstMu, 1e-12);

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


    // ---------------------------------------------------------------------------------------------------
    // The light-curve model as the simulator calls it.
    // ---------------------------------------------------------------------------------------------------
    // A representative event: the ModelEvent the adapter takes, and the same in the simulator's structs.
    // thetaE 0.6 mas at Ds = 8 kpc, Dl = 4 kpc: pi_rel = 0.125 mas, piE = 0.208.
    struct Ev {
        source s;
        lens l;
        astromet as;
        SourceModel sm = SourceModel::PSPL;
        ModelEvent m() const { return ModelEvent{l.u0, l.tE, l.t0, l.piE, s.xi, l.tetE, s.mus1, s.mus2, 1.0 / s.Ds, s.rho, sm}; }
    };
    // Every event is a point source (s.rho = 0, the point-source model) unless a test says otherwise; the
    // seams read the run's source model, which the first tests leave at the point source.
    auto fill = [](Ev& e, double u0, double tE, double t0, double piE, double xi, double tetE) {
        e.l.u0 = u0;  e.l.tE = tE;  e.l.t0 = t0;  e.l.piE = piE;  e.l.tetE = tetE;  e.s.xi = xi;
        e.s.mus1 = 0.004;  e.s.mus2 = -0.003;  e.s.Ds = 8.0;
        e.s.fb = {0.6, 0.7};  e.s.fLens = {0.0, 0.0};  e.s.blendOff = {};
        e.s.rho = 0.0;  e.sm = SourceModel::PSPL;
        e.as.satScale = 1.0;
    };
    auto ev = std::make_unique<Ev>();
    vbmSetSourceModel(SourceModel::PSPL);

    // (f) The side of the lens the source passes, and the handedness of xi. With the parallax switched off
    // (piE = 1e-9) the library's source-lens vector must be the simulator's straight line on (e_l, e_b),
    //     u0 (-sin xi, cos xi) + tau (cos xi, sin xi),    tau = (t - t0) / tE,
    // for every xi and sightline. At t0 the two coincide to rounding (the Earth is at rest in the frame of
    // the line, and its light-travel time to the Sun is the one subtracted). Away from t0 the library reads
    // the clock as the time light reaches the Sun, t + (Earth's position along the line of sight)/c, so its
    // along-track coordinate carries the change of that delay, at most 2 AU/c = 0.0116 d, over tE: nothing in
    // the legacy straight line, which was in the Earth's own time. The coordinate across the track, u0, is
    // untouched. The mirror line (u0 of the other sign) is 2 u0 away, printed so that the test is seen to
    // discriminate.
    {
        const double LT_MAX = 2.0 * 0.005775518331436995;     // 2 AU / c [days]
        const double tE = 25.0, u0 = 0.37;
        double atT0 = 0.0, across = 0.0, along = 0.0, worstMirror = 1e300;
        for (int i = 0; i < N_SIGHT; ++i) {
            vbmSetSightline(GAL_L[i], GAL_B[i]);
            for (double xi : {0.3, 1.4, 2.9, 4.1, 5.6}) {
                fill(*ev, u0, tE, 1500.0, 1e-9, xi, 0.6);
                const double tt[7] = {1470.0, 1482.0, 1495.0, 1500.0, 1507.0, 1519.0, 1530.0};
                ModelEpoch out[7];
                vbmEvaluate(ev->m(), false, false, tt, 7, out);
                for (int k = 0; k < 7; ++k) {
                    const double tau = (tt[k] - 1500.0) / tE;
                    const double lx = -u0 * std::sin(xi) + tau * std::cos(xi), ly = u0 * std::cos(xi) + tau * std::sin(xi);
                    if (tt[k] == 1500.0) atT0 = std::fmax(atT0, std::hypot(out[k].ux - lx, out[k].uy - ly));
                    across = std::fmax(across, std::fabs(-std::sin(xi) * (out[k].ux - lx) + std::cos(xi) * (out[k].uy - ly)));
                    along  = std::fmax(along,  std::fabs(std::cos(xi) * (out[k].ux - lx) + std::sin(xi) * (out[k].uy - ly)) * tE / LT_MAX);
                    const double mx = +u0 * std::sin(xi) + tau * std::cos(xi), my = -u0 * std::cos(xi) + tau * std::sin(xi);
                    worstMirror = std::fmin(worstMirror, std::hypot(out[k].ux - mx, out[k].uy - my));
                }
            }
        }
        bound("(f) source-lens vector at t0 (piE = 1e-9) vs u0 (-sin xi, cos xi) on (e_l, e_b): max |difference|", atT0, 1e-12);
        bound("(f) ... across the track at |t - t0| <= 30 d vs the straight line: max |difference| [Einstein radii]", across, 1e-9);
        bound("(f) ... along the track vs the straight line, in units of the light-travel bound 2 AU/(c tE): max", along, 1.0);
        expect(worstMirror > 0.3, "(f) ... and it is far from the mirror line (u0 of the other sign): min distance", worstMirror, ">", 0.3);
    }

    // (g) A is magnifOf(u) bit for bit, in the library's own array, in both seams and for both observers.
    {
        int same = 0, total = 0;
        for (int i = 0; i < N_SIGHT; ++i) {
            vbmSetSightline(GAL_L[i], GAL_B[i]);
            fill(*ev, 0.21, 18.0, 1500.0, 0.21, 2.2, 0.6);
            std::vector<double> tt(300);
            for (int k = 0; k < 300; ++k) tt[k] = 1380.0 + 0.8 * k;
            std::vector<ModelEpoch> me(300);
            std::vector<TrackPoint> tp(300);
            std::vector<ModelPoint> mp(300);
            for (int tele = 0; tele < 2; ++tele) {
                vbmEvaluate(ev->m(), tele == 1, true, tt.data(), 300, me.data());
                evaluateTrack(ev->s, ev->l, ev->as, tele, tt.data(), 300, true, tp.data());
                for (int k = 0; k < 300; ++k) { total += 3; same += (me[k].A == magnifOf(me[k].u)) + (tp[k].A == magnifOf(tp[k].u)); }
                evaluateModel(ev->s, ev->l, ev->as, tele, tt.data(), 300, false, mp.data());
                for (int k = 0; k < 300; ++k) same += (mp[k].A == magnifOf(mp[k].u));
            }
        }
        expect(same == total, "(g) A == magnifOf(u) bitwise: vbmEvaluate, evaluateTrack, evaluateModel", same, "=", total);
    }

    // (h) The centroid of the two images lies thetaE u/(u^2 + 2) from the unlensed source along the
    // source-lens vector (u = source - lens), and the lens lies -thetaE u from it. Held at the first epoch of
    // a call and at the 3000th of one: the library perturbs piE_E by 1e-10 per epoch of an astrometric call
    // (~15 microradians of direction at piE = 0.02 after 3000), which the adapter follows.
    {
        double worst1 = 0.0, worstN = 0.0, worstLens = 0.0;
        for (int i = 0; i < N_SIGHT; ++i) {
            vbmSetSightline(GAL_L[i], GAL_B[i]);
            for (int tele = 0; tele < 2; ++tele) {
                fill(*ev, 0.4, 30.0, 1500.0, 0.3, 1.3, 0.6);
                for (double t : {1450.0, 1490.0, 1500.0, 1512.0, 1560.0}) {
                    ModelEpoch o;
                    vbmEvaluate(ev->m(), tele == 1, true, &t, 1, &o);
                    const double d = 0.6 / (o.u * o.u + 2.0);
                    worst1 = std::fmax(worst1, std::hypot(o.imgx - o.srcx - d * o.ux, o.imgy - o.srcy - d * o.uy) / (d * o.u));
                    worstLens = std::fmax(worstLens, std::hypot(o.lensx - o.srcx + 0.6 * o.ux, o.lensy - o.srcy + 0.6 * o.uy));
                }
                fill(*ev, 0.4, 30.0, 1500.0, 0.02, 1.3, 0.6);
                std::vector<double> tt(3000);
                for (int k = 0; k < 3000; ++k) tt[k] = 1400.0 + 0.07 * k;
                std::vector<ModelEpoch> o(3000);
                vbmEvaluate(ev->m(), tele == 1, true, tt.data(), 3000, o.data());
                for (int k = 0; k < 3000; ++k) {
                    const double d = 0.6 / (o[k].u * o[k].u + 2.0);
                    worstN = std::fmax(worstN, std::hypot(o[k].imgx - o[k].srcx - d * o[k].ux, o[k].imgy - o[k].srcy - d * o[k].uy) / (d * o[k].u));
                    worstLens = std::fmax(worstLens, std::hypot(o[k].lensx - o[k].srcx + 0.6 * o[k].ux, o[k].lensy - o[k].srcy + 0.6 * o[k].uy));
                }
            }
        }
        bound("(h) image centroid - unlensed source vs thetaE u/(u^2 + 2): max relative difference, 1 epoch per call", worst1, 1e-9);
        bound("(h) ... the same at 3000 epochs of one call (piE = 0.02)", worstN, 1e-9);
        bound("(h) lens - unlensed source = -thetaE u: max difference [mas]", worstLens, 1e-9);
    }

    // (i) The heliocentric-to-geocentric conversion against numpy with astropy's matrix and Earth velocity.
    // The library's velocity is its daily table (1.3e-4 of the Earth's speed off astropy's, test (d)), which
    // reaches tE through v_E pi_rel / |mu_geo|: up to ~6 for the third event, whose heliocentric motion is no
    // larger than the Earth's contribution (the geocentric tE there is ~6 times the heliocentric one).
    {
        double worstTe = 0.0, worstXiG = 0.0, ratioMax = 0.0;
        for (int i = 0; i < N_SIGHT; ++i) {
            vbmSetSightline(GAL_L[i], GAL_B[i]);
            for (int k = 0; k < N_T; ++k)
                for (int e = 0; e < N_GEO; ++e) {
                    const double* g = GEO_EVENT[e];
                    const GeocentricMotion m = geocentricMotion(T_DAY[k], g[5], g[4], g[0], g[1], g[2], g[3]);
                    worstTe = std::fmax(worstTe, std::fabs(m.tE / GEO_TE[i][k][e] - 1.0));
                    worstXiG = std::fmax(worstXiG, std::fabs(std::remainder(m.xi - GEO_XI[i][k][e], 2.0 * pi)));
                    ratioMax = std::fmax(ratioMax, m.tE / (g[5] / std::hypot(g[2] - g[0], g[3] - g[1])));
                }
        }
        std::printf("  (largest geocentric/heliocentric tE in the constructed events: %.1f)\n", ratioMax);
        bound("(i) geocentric tE vs numpy with astropy's Earth velocity: max relative difference", worstTe, 3e-4);
        bound("(i) direction of the geocentric source-lens motion xi vs numpy: max |dxi| [rad]", worstXiG, 8e-4);
        expect(ratioMax > 5.0, "(i) ... and the conversion is not a small effect for them (largest tE_geo / tE_hel)", ratioMax, ">", 5.0);
    }

    // (j) The Roman-Earth separation of the source-lens vectors at t0, du_sat: piE times the L2 offset
    // projected on the sky. Reference (e): |dy| = |piE| |s|; and as a vector, on (e_l, e_b), the track shift
    // dy (along and across piE_vec) rotated into the frame of the sightline.
    {
        double worstNorm = 0.0, worstVec = 0.0, zero = 0.0;
        for (int i = 0; i < N_SIGHT; ++i) {
            const SkyFrame fr = vbmSetSightline(GAL_L[i], GAL_B[i]);
            for (int k = 0; k < N_T; ++k)
                for (int ip = 0; ip < N_PI; ++ip) {
                    const double pn = PIE[ip][0], pe = PIE[ip][1], pm = std::hypot(pn, pe);
                    // xi with piE_vec = piE (-toNE (cos xi, sin xi)): (cos xi, sin xi) = -toNE^T piE_hat
                    const double xi = std::atan2(-(fr.toNE[0][1] * pn + fr.toNE[1][1] * pe),
                                                 -(fr.toNE[0][0] * pn + fr.toNE[1][0] * pe));
                    fill(*ev, 0.3, 40.0, T_DAY[k], pm, xi, 0.6);
                    TrackPoint e0, e1;
                    evaluateTrack(ev->s, ev->l, ev->as, 0, &ev->l.t0, 1, true, &e0);
                    evaluateTrack(ev->s, ev->l, ev->as, 1, &ev->l.t0, 1, true, &e1);
                    const double* ref = L2_SHIFT[i][k][ip];
                    const double rn = std::hypot(ref[0], ref[1]);
                    worstNorm = std::fmax(worstNorm, std::fabs(std::hypot(e1.ux - e0.ux, e1.uy - e0.uy) / rn - 1.0));
                    // dy = (y1, y2) along piE_hat and a quarter turn counterclockwise on (N, E) from it
                    const double hN = pn / pm, hE = pe / pm;
                    const double vN = ref[0] * hN - ref[1] * hE, vE = ref[0] * hE + ref[1] * hN;
                    const double vl = fr.toNE[0][0] * vN + fr.toNE[1][0] * vE, vb = fr.toNE[0][1] * vN + fr.toNE[1][1] * vE;
                    worstVec = std::fmax(worstVec, std::hypot(e1.ux - e0.ux - vl, e1.uy - e0.uy - vb) / rn);
                    // Roman at the centre of the Earth is the Earth
                    ev->as.satScale = 0.0;
                    TrackPoint e2;
                    evaluateTrack(ev->s, ev->l, ev->as, 1, &ev->l.t0, 1, true, &e2);
                    zero = std::fmax(zero, std::hypot(e2.ux - e0.ux, e2.uy - e0.uy));
                }
        }
        bound("(j) |u_Roman - u_Earth| at t0 (du_sat) vs piE x L2 offset projected: max relative difference", worstNorm, 1e-4);
        bound("(j) u_Roman - u_Earth as a vector on (e_l, e_b) vs the reference rotated into place: max, relative to |ref|", worstVec, 1e-4);
        bound("(j) satScale = 0: Roman's vector minus the Earth's, max |difference|", zero, 0.0 + 1e-300);
    }

    // (k) The annual parallax. The track's departure from the geocentric straight line is piE times the
    // observer's departure from its own straight line at t0, X(t) - X(t0) - V(t0) (t - t0), on (e_l, e_b); adding
    // back the Earth's velocity at t0 (earthVelocityNE, test (d)) leaves the orbit X(t) - X(t0), which is
    // compared with astropy's, point by point over a year (to the light-travel delay the library applies to the
    // clock, up to 2 AU/c x (v_E + 1/(piE tE)) = 3e-4 AU here, and its daily table). The orbit's shape then follows from its principal
    // axes: the sightline's ecliptic latitude makes it a thin ellipse (axis ratio ~sin|beta|, 0.094 at (1, -1.5)),
    // its major axis ~120 deg from +l. The legacy model's circle tilted by a fixed 60 deg had ratio ~0.88.
    {
        double worstOrbit = 0.0, worstRatio = 0.0, worstPA = 0.0, ratioMin = 1e300;
        for (int i = 0; i < N_SIGHT; ++i) {
            const SkyFrame fr = vbmSetSightline(GAL_L[i], GAL_B[i]);
            const double t0 = T_DAY[1], piE = 0.5, xi = 1.0, u0 = 0.5, tE = 400.0;
            fill(*ev, u0, tE, t0, piE, xi, 0.6);
            double vNE[2];
            earthVelocityNE(t0, vNE);
            const double vl = fr.toNE[0][0] * vNE[0] + fr.toNE[1][0] * vNE[1], vb = fr.toNE[0][1] * vNE[0] + fr.toNE[1][1] * vNE[1];
            auto orbit = [&](const double* tt, int n, double* X) {      // X(t) - X(t0), (l, b), AU
                std::vector<ModelEpoch> o(n);
                vbmEvaluate(ev->m(), false, false, tt, n, o.data());
                for (int k = 0; k < n; ++k) {
                    const double tau = (tt[k] - t0) / tE;
                    const double dx = o[k].ux - (-u0 * std::sin(xi) + tau * std::cos(xi));
                    const double dy = o[k].uy - (u0 * std::cos(xi) + tau * std::sin(xi));
                    X[2 * k]     = dx / piE + vl * (tt[k] - t0);
                    X[2 * k + 1] = dy / piE + vb * (tt[k] - t0);
                }
            };
            double tt[N_DEV], X[2 * N_DEV];
            for (int k = 0; k < N_DEV; ++k) tt[k] = t0 + DEV_OFFSET[k];
            orbit(tt, N_DEV, X);
            for (int k = 0; k < N_DEV; ++k)
                worstOrbit = std::fmax(worstOrbit, std::hypot(X[2 * k] - ORBIT_REF[i][k][0], X[2 * k + 1] - ORBIT_REF[i][k][1]));
            std::vector<double> t2(ELL_N), Y(2 * ELL_N);
            for (int k = 0; k < ELL_N; ++k) t2[k] = t0 + ELL_PERIOD * k / ELL_N;
            orbit(t2.data(), ELL_N, Y.data());
            double mx = 0.0, my = 0.0;
            for (int k = 0; k < ELL_N; ++k) { mx += Y[2 * k]; my += Y[2 * k + 1]; }
            mx /= ELL_N;  my /= ELL_N;
            double cxx = 0.0, cxy = 0.0, cyy = 0.0;
            for (int k = 0; k < ELL_N; ++k) {
                const double a = Y[2 * k] - mx, b = Y[2 * k + 1] - my;
                cxx += a * a;  cxy += a * b;  cyy += b * b;
            }
            cxx /= ELL_N;  cxy /= ELL_N;  cyy /= ELL_N;
            const double tr = cxx + cyy, df = std::hypot(cxx - cyy, 2.0 * cxy);
            const double wMax = 0.5 * (tr + df), wMin = 0.5 * (tr - df);
            const double ratio = std::sqrt(wMin / wMax);
            double pa = 0.5 * std::atan2(2.0 * cxy, cxx - cyy) * RAa;      // major axis from +l toward +b
            pa = std::fmod(pa + 180.0, 180.0);
            worstRatio = std::fmax(worstRatio, std::fabs(ratio - ELL_RATIO[i]));
            worstPA = std::fmax(worstPA, std::fabs(std::remainder(pa - ELL_PA[i], 180.0)));
            ratioMin = std::fmin(ratioMin, 0.3 - ratio);
            if (i == 0) std::printf("  (l, b) = (%.1f, %.1f): orbit axis ratio %.4f (astropy %.4f), major axis %.2f deg from +l (astropy %.2f)\n",
                                    GAL_L[i], GAL_B[i], ratio, ELL_RATIO[i], pa, ELL_PA[i]);
        }
        bound("(k) Earth's orbit X(t) - X(t0) from the track vs astropy, over a year: max |difference| [AU]", worstOrbit, 3e-4);
        bound("(k) orbit axis ratio vs astropy's: max |difference|", worstRatio, 1e-3);
        bound("(k) orbit major-axis angle vs astropy's: max |difference| [deg]", worstPA, 0.05);
        expect(ratioMin > 0.0, "(k) every orbit is a thin ellipse (the legacy circle had 0.88): 0.3 - max axis ratio", ratioMin, ">", 0.0);
    }

    // (l) What is the source and what is the lens. (1) The unlensed source's track is the source's own: the same,
    // relative to its position at the first epoch (the origin of the positions is the lens at t0, which moves
    // with u0 and thetaE), for two events that differ in every lens parameter (thetaE, tE, u0, xi, piE). (2) Over a sidereal year the
    // Earth returns to its place, so the parallax cancels and source and lens move by their HELIOCENTRIC proper
    // motions: the source's mus, the lens's mul, which func_lens drew before converting tE and xi to the
    // geocentric frame. This ties the conversion (the sign and size of v_E pi_rel, the axes, the unit of the
    // library's proper motions) to the heliocentric velocities it started from; without the v_E pi_rel term the
    // lens would be off by ~6 pi_rel mas, here ~0.8 mas.
    {
        double worstSrc = 0.0, worstMuS = 0.0, worstMuL = 0.0, wrongLens = 1e300;
        const double Y = 365.256363;                                  // sidereal year [days]
        for (int i = 0; i < N_SIGHT; ++i) {
            vbmSetSightline(GAL_L[i], GAL_B[i]);
            Ev a, b;
            fill(a, 0.3, 40.0, 1500.0, 0.25, 1.1, 0.6);
            fill(b, 0.9, 15.0, 1500.0, 0.05, 4.0, 2.0);
            std::vector<double> tt(100);
            for (int k = 0; k < 100; ++k) tt[k] = 1300.0 + 4.0 * k;
            std::vector<ModelEpoch> oa(100), ob(100);
            vbmEvaluate(a.m(), false, true, tt.data(), 100, oa.data());
            vbmEvaluate(b.m(), false, true, tt.data(), 100, ob.data());
            for (int k = 0; k < 100; ++k)
                worstSrc = std::fmax(worstSrc, std::hypot((oa[k].srcx - oa[0].srcx) - (ob[k].srcx - ob[0].srcx),
                                                          (oa[k].srcy - oa[0].srcy) - (ob[k].srcy - ob[0].srcy)));

            // a realistic event: heliocentric proper motions as drawn (mas/day), then converted
            fill(a, 0.3, 40.0, 1500.0, 0.0, 0.0, 0.6);
            a.s.mus1 = 0.05;  a.s.mus2 = -0.03;  a.l.mul1 = -0.03;  a.l.mul2 = 0.02;
            a.s.Ds = 8.0;  a.l.Dl = 4.0;
            a.l.pirel = 1.0 / a.l.Dl - 1.0 / a.s.Ds;
            a.l.piE = a.l.pirel / a.l.tetE;
            const GeocentricMotion g = geocentricMotion(a.l.t0, a.l.tetE, a.l.pirel, a.s.mus1, a.s.mus2, a.l.mul1, a.l.mul2);
            a.l.tE = g.tE;  a.s.xi = g.xi;
            const double t2[2] = {1500.0, 1500.0 + Y};
            ModelEpoch o[2];
            vbmEvaluate(a.m(), false, true, t2, 2, o);
            worstMuS = std::fmax(worstMuS, std::hypot(o[1].srcx - o[0].srcx - a.s.mus1 * Y, o[1].srcy - o[0].srcy - a.s.mus2 * Y));
            worstMuL = std::fmax(worstMuL, std::hypot(o[1].lensx - o[0].lensx - a.l.mul1 * Y, o[1].lensy - o[0].lensy - a.l.mul2 * Y));
            // the heliocentric conversion left out (tE and xi from mu_hel itself): the lens is then off
            const GeocentricMotion h = geocentricMotion(a.l.t0, a.l.tetE, 0.0, a.s.mus1, a.s.mus2, a.l.mul1, a.l.mul2);
            Ev w = a;  w.l.tE = h.tE;  w.s.xi = h.xi;
            ModelEpoch ow[2];
            vbmEvaluate(w.m(), false, true, t2, 2, ow);
            wrongLens = std::fmin(wrongLens, std::hypot(ow[1].lensx - ow[0].lensx - a.l.mul1 * Y, ow[1].lensy - ow[0].lensy - a.l.mul2 * Y));
        }
        bound("(l) unlensed source track relative to its first epoch, two events with different lens parameters: max |difference| [mas]", worstSrc, 1e-6);
        bound("(l) source displacement over a sidereal year vs mus x year: max |difference| [mas] (of ~21)", worstMuS, 1e-3);
        bound("(l) lens displacement over a sidereal year vs mul x year, after the conversion: max |difference| [mas] (of ~13)", worstMuL, 1e-3);
        expect(wrongLens > 0.3, "(l) ... and with tE, xi taken as heliocentric the lens would be off by: min [mas]", wrongLens, ">", 0.3);
    }

    // (m) Where Roman is. satScale is a switch: 0 puts Roman at the centre of the Earth, bit for bit Rubin's
    // geometry; 1 puts it at L2; nothing between is allowed.
    {
        vbmSetSightline(GAL_L[0], GAL_B[0]);
        fill(*ev, 0.3, 40.0, 1500.0, 0.25, 1.1, 0.6);
        std::vector<double> tt(50);
        for (int k = 0; k < 50; ++k) tt[k] = 1450.0 + k;
        std::vector<TrackPoint> r(50), q1(50), q0(50);
        evaluateTrack(ev->s, ev->l, ev->as, 0, tt.data(), 50, true, r.data());
        ev->as.satScale = 1.0;
        evaluateTrack(ev->s, ev->l, ev->as, 1, tt.data(), 50, true, q1.data());
        ev->as.satScale = 0.0;
        evaluateTrack(ev->s, ev->l, ev->as, 1, tt.data(), 50, true, q0.data());
        int same = 0, differ = 0;
        for (int k = 0; k < 50; ++k) {
            same += (q0[k].u == r[k].u) + (q0[k].pos1b == r[k].pos1b) + (q0[k].pos2b == r[k].pos2b) + (q0[k].lens1 == r[k].lens1);
            differ += (q1[k].u != r[k].u);
        }
        expect(same == 200, "(m) satScale = 0: Roman's u, source and lens positions are Rubin's, bitwise (of 200)", same, "=", 200.0);
        expect(differ == 50, "(m) satScale = 1: Roman's u differs from Rubin's at every epoch (of 50)", differ, "=", 50.0);
        bool threw = false;
        ev->as.satScale = 0.5;
        try { evaluateTrack(ev->s, ev->l, ev->as, 1, tt.data(), 50, true, q1.data()); } catch (const std::runtime_error&) { threw = true; }
        expect(threw, "(m) satScale = 0.5 is refused", threw, "=", 1.0);
        ev->as.satScale = 1.0;
    }

    // (n) The seams agree with each other and with the straight line. evaluateModel without centroids is the
    // photometric track alone (no per-epoch drift of piE_E), which is also where evaluateTrack's u comes from:
    // equal, bit for bit; and with centroids it gives evaluateTrack's measured centroid, bit for bit. "No parallax" is the geocentric line the
    // event is parametrised by: at t0 it is the closest approach, u0, and with the parallax off the track is on it
    // up to the light-travel delay the library applies to the clock (see (f)): |u - u_noplx| <= 2 AU/(c tE).
    {
        double worstU = 0.0, worstLine = 0.0, atT0 = 0.0, worstDef = 0.0;
        int posSame = 0, nPos = 0;
        for (int i = 0; i < N_SIGHT; ++i) {
            vbmSetSightline(GAL_L[i], GAL_B[i]);
            fill(*ev, 0.3, 40.0, 1500.0, 0.25, 1.1, 0.6);
            ev->s.fb = {0.5, 0.8};  ev->s.fLens = {0.1, 0.05};  ev->s.blendOff = {{{0.3, -0.2}, {0.1, 0.4}}};
            std::vector<double> tt(200);
            for (int k = 0; k < 200; ++k) tt[k] = 1400.0 + k;
            std::vector<TrackPoint> tp(200);
            std::vector<ModelPoint> ph(200), as(200);
            for (int tele = 0; tele < 2; ++tele) {
                evaluateTrack(ev->s, ev->l, ev->as, tele, tt.data(), 200, true, tp.data());
                evaluateModel(ev->s, ev->l, ev->as, tele, tt.data(), 200, false, ph.data());
                evaluateModel(ev->s, ev->l, ev->as, tele, tt.data(), 200, true, as.data());
                for (int k = 0; k < 200; ++k) {
                    worstU = std::fmax(worstU, std::fabs(ph[k].u / tp[k].u - 1.0));
                    nPos += 2;
                    posSame += (as[k].pos1c == tp[k].pos1c) + (as[k].pos2c == tp[k].pos2c);
                }
            }
            // the line: u0 at t0, the deflection thetaE u/(u^2 + 2) there (u0 along (-sin xi, cos xi))
            const double t0 = 1500.0;
            TrackPoint p;
            evaluateTrack(ev->s, ev->l, ev->as, 0, &t0, 1, true, &p);
            atT0 = std::fmax(atT0, std::fabs(p.uNoPlx - 0.3));
            worstDef = std::fmax(worstDef, std::hypot(p.def1a + 0.6 * 0.3 * std::sin(1.1) / (0.3 * 0.3 + 2.0),
                                                       p.def2a - 0.6 * 0.3 * std::cos(1.1) / (0.3 * 0.3 + 2.0)));
            // parallax off: the track is on the line
            fill(*ev, 0.3, 40.0, 1500.0, 1e-9, 1.1, 0.6);
            ModelPoint ml[200];
            std::vector<TrackPoint> tl(200);
            evaluateTrack(ev->s, ev->l, ev->as, 1, tt.data(), 200, true, tl.data());
            for (int k = 0; k < 200; ++k) worstLine = std::fmax(worstLine, std::fabs(tl[k].u - tl[k].uNoPlx));
            (void)ml;
        }
        bound("(n) evaluateModel (photometry only) vs evaluateTrack: u, max relative difference", worstU, 1e-300);
        expect(posSame == nPos, "(n) evaluateModel with centroids == evaluateTrack's measured centroid, bitwise", posSame, "=", nPos);
        bound("(n) uNoPlx at t0 vs u0: max |difference|", atT0, 1e-15);
        bound("(n) def_noplx at t0 vs thetaE u0/(u0^2 + 2) along (-sin xi, cos xi): max |difference| [mas]", worstDef, 1e-14);
        bound("(n) at piE = 1e-9 the track is the no-parallax line, to the light-travel bound 2 AU/(c tE): max |u - uNoPlx|", worstLine, 2.0 * 0.005775518331436995 / 40.0);
    }

    // (o) evaluateTrack with and without the absolute positions (the time loop asks for them only when a sample
    // dump is written, and a run must not depend on that). u, A, ux, uy, uNoPlx, ANoPlx and the no-parallax
    // deflection are the photometric call's either way, so bitwise equal, and equal evaluateModel's. The
    // positions come from a second, astrometric call, whose track differs from the photometric one by the
    // library's perturbation of piE_E along the array (1e-10 an epoch, 1e-7 at the 1024th epoch, the time
    // loop's block): printed here for the call itself, and held for the deflection, which also turns by up to
    // 1.01e-10 x 1024 / piE = 4e-7 rad with the astrometric call's frame. Without positions the deflection is
    // thetaE (ux, uy)/(u^2 + 2), the point lens's, and the six positions are NaN.
    {
        double worstCall = 0.0, worstDef = 0.0, worstDefOwn = 0.0;
        int same = 0, nSame = 0, photSame = 0, nPhot = 0, isNan = 0, nNan = 0;
        for (int i = 0; i < N_SIGHT; ++i) {
            vbmSetSightline(GAL_L[i], GAL_B[i]);
            fill(*ev, 0.3, 40.0, 1500.0, 0.25, 1.1, 0.6);
            const int n = 1024;
            std::vector<double> tt(n);
            for (int k = 0; k < n; ++k) tt[k] = 1100.0 + 0.9 * k;
            std::vector<TrackPoint> on(n), off(n);
            std::vector<ModelPoint> ph(n);
            std::vector<ModelEpoch> pc(n), ac(n);
            for (int tele = 0; tele < 2; ++tele) {
                evaluateTrack(ev->s, ev->l, ev->as, tele, tt.data(), n, true, on.data());
                evaluateTrack(ev->s, ev->l, ev->as, tele, tt.data(), n, false, off.data());
                evaluateModel(ev->s, ev->l, ev->as, tele, tt.data(), n, false, ph.data());
                vbmEvaluate(ev->m(), tele == 1, false, tt.data(), n, pc.data());
                vbmEvaluate(ev->m(), tele == 1, true, tt.data(), n, ac.data());
                for (int k = 0; k < n; ++k) {
                    nSame += 8;  same += (off[k].u == on[k].u) + (off[k].A == on[k].A) + (off[k].ux == on[k].ux)
                                       + (off[k].uy == on[k].uy) + (off[k].uNoPlx == on[k].uNoPlx)
                                       + (off[k].ANoPlx == on[k].ANoPlx) + (off[k].def1a == on[k].def1a)
                                       + (off[k].def2a == on[k].def2a);
                    nPhot += 2;  photSame += (on[k].u == ph[k].u) + (on[k].A == ph[k].A);
                    worstCall = std::fmax(worstCall, std::fmax(std::fabs(ac[k].u / pc[k].u - 1.0), std::fabs(ac[k].A / pc[k].A - 1.0)));
                    const double d = 0.6 / (off[k].u * off[k].u + 2.0);
                    worstDef = std::fmax(worstDef, std::hypot(off[k].def1c - on[k].def1c, off[k].def2c - on[k].def2c) / (d * off[k].u));
                    worstDefOwn = std::fmax(worstDefOwn, std::hypot(off[k].def1c - d * off[k].ux, off[k].def2c - d * off[k].uy) / (d * off[k].u));
                    const double pos[6] = {off[k].pos1b, off[k].pos2b, off[k].pos1c, off[k].pos2c, off[k].lens1, off[k].lens2};
                    for (double x : pos) { ++nNan; isNan += std::isnan(x); }
                }
            }
        }
        expect(same == nSame, "(o) positions on vs off: u, A, ux, uy, uNoPlx, ANoPlx, def1a, def2a are bitwise equal", same, "=", nSame);
        expect(photSame == nPhot, "(o) ... and u, A are evaluateModel's (the photometric call), bitwise", photSame, "=", nPhot);
        bound("(o) (the library's astrometric call itself: u and A vs its photometric call, max relative difference)", worstCall, 1.5e-7);
        bound("(o) deflection with positions (library) vs without (point lens): max |difference| / (thetaE u/(u^2 + 2))", worstDef, 1e-6);
        bound("(o) positions off: deflection vs thetaE (ux, uy)/(u^2 + 2), max relative difference", worstDefOwn, 1e-15);
        expect(isNan == nNan, "(o) positions off: source, lens and centroid positions are NaN", isNan, "=", nNan);
    }

    // (p) Epochs outside the ephemerides are refused rather than extrapolated: the Sun table well before
    // the survey, and Roman's table before its first row (the Earth's alone does not stop an Earth observer).
    {
        vbmSetSightline(GAL_L[0], GAL_B[0]);
        fill(*ev, 0.3, 40.0, 1500.0, 0.25, 1.1, 0.6);
        double first, last;
        vbmEarthSpan(&first, &last);
        const double bad[2] = {last + 500.0, first - 500.0};
        int refused = 0;
        for (double t : bad) {
            try { ModelEpoch o; vbmEvaluate(ev->m(), false, false, &t, 1, &o); } catch (const std::runtime_error&) { ++refused; }
        }
        const double tRoman = -400.0;
        try { ModelEpoch o; vbmEvaluate(ev->m(), true, false, &tRoman, 1, &o); } catch (const std::runtime_error&) { ++refused; }
        expect(refused == 3, "(p) Earth past either end of the Sun table, Roman before its table: refused (of 3)", refused, "=", 3.0);
        std::printf("  (the Earth's ephemeris covers simulation days %.0f to %.0f)\n", first, last);
    }

    // ---------------------------------------------------------------------------------------------------
    // The finite source: a uniform disc (--source-model espl). The seams read the run's source model and the
    // source's rho; the direct references are the library's own ESPLMag2 (a second object of it, set up as the
    // adapter's is: uniform disc, Tol = VBM_TOL), the exact values of tests/vbm_reference.h, and the
    // quadrature uniformDisc above, which needs nothing of the library.
    // ---------------------------------------------------------------------------------------------------
    auto lib = std::make_unique<VBMicrolensing>();
    lib->LoadESPLTable(PATH_VBM_ESPL_TABLE);
    lib->Tol = VBM_TOL;
    lib->a1 = 0.0;
    vbmSetSightline(GAL_L[0], GAL_B[0]);

    // (q) The adapter's magnification is the library's ESPLMag2 at the track's u, and the exact uniform-disc
    // magnification to the accuracy vbmtest holds the library to (1e-4). The grid is vbmtest's: rho 1e-3..0.1,
    // z = u/rho 0.01..30, at t0 of an event with the parallax off, where the track is at u0 = z rho. The
    // quadrature is checked against the same references first, which is what lets it stand in for them below.
    {
        int libSame = 0, seamSame = 0, nGrid = 0;
        double worstRef = 0.0, worstInt = 0.0, worstQuad = 0.0;
        vbmSetSourceModel(SourceModel::ESPL);
        for (int ir = 0; ir < vbmref::N_RHO; ++ir)
            for (int iz = 0; iz < vbmref::N_Z; ++iz) {
                const double rho = vbmref::RHO[ir], ref = vbmref::ESPL[0][ir][iz];
                fill(*ev, vbmref::Z[iz] * rho, 25.0, 1500.0, 1e-9, 1.1, 0.6);
                ev->s.rho = rho;  ev->sm = SourceModel::ESPL;
                const double t = ev->l.t0;
                ModelEpoch o;
                vbmEvaluate(ev->m(), false, false, &t, 1, &o);
                ++nGrid;
                libSame += (o.A == lib->ESPLMag2(o.u, std::exp(std::log(rho))));
                worstRef = std::fmax(worstRef, std::fabs(o.A / ref - 1.0));
                const DiscRef q = uniformDisc(o.u, rho);
                worstInt = std::fmax(worstInt, std::fabs(o.A / q.A - 1.0));
                worstQuad = std::fmax(worstQuad, std::fabs(q.A / ref - 1.0));
                TrackPoint tp;  ModelPoint mp;
                evaluateTrack(ev->s, ev->l, ev->as, 0, &t, 1, true, &tp);
                evaluateModel(ev->s, ev->l, ev->as, 0, &t, 1, false, &mp);
                seamSame += (tp.A == o.A) + (mp.A == o.A);
            }
        bound("(q) quadrature vs tests/vbm_reference.h, uniform disc, 65 (rho, z): max relative difference", worstQuad, 1e-7);
        expect(libSame == nGrid, "(q) adapter A == the library's ESPLMag2(u, rho), bitwise, of 65 (rho, z)", libSame, "=", nGrid);
        bound("(q) adapter A vs the exact uniform-disc magnification (references): max relative difference", worstRef, 1e-4);
        bound("(q) adapter A vs the exact uniform-disc magnification (quadrature): max relative difference", worstInt, 1e-4);
        expect(seamSame == 2 * nGrid, "(q) evaluateTrack and evaluateModel give the same A, bitwise, of 130", seamSame, "=", 2 * nGrid);
    }

    // (r) The finite source tends to the point source with the radius: at rho = 1e-6 (and 1e-4, printed: the
    // difference is the library's quadrupole term, ~rho^2 at u ~ 0.01) along a whole light curve with
    // parallax, in both observers.
    {
        double worst6 = 0.0, worst4 = 0.0, worstOwn = 0.0;
        for (double u0 : {0.01, 0.1, 0.3, 1.0, 3.0}) {
            fill(*ev, u0, 30.0, 1500.0, 0.25, 1.1, 0.6);
            std::vector<double> tt(300);
            for (int k = 0; k < 300; ++k) tt[k] = 1380.0 + 0.8 * k;
            std::vector<ModelEpoch> pt(300), fs(300);
            for (int tele = 0; tele < 2; ++tele)
                for (double rho : {1e-6, 1e-4}) {
                    ev->s.rho = rho;
                    ev->sm = SourceModel::PSPL;  vbmEvaluate(ev->m(), tele == 1, true, tt.data(), 300, pt.data());
                    ev->sm = SourceModel::ESPL;  vbmEvaluate(ev->m(), tele == 1, true, tt.data(), 300, fs.data());
                    for (int k = 0; k < 300; ++k) {
                        const double d = std::fabs(fs[k].A / pt[k].A - 1.0);
                        (rho < 1e-5 ? worst6 : worst4) = std::fmax(rho < 1e-5 ? worst6 : worst4, d);
                        worstOwn = std::fmax(worstOwn, std::fabs(fs[k].u / pt[k].u - 1.0) + std::fabs(fs[k].ux - pt[k].ux));
                    }
                }
        }
        bound("(r) rho = 1e-6: A of the disc vs the point source, max relative difference", worst6, 5e-9);
        std::printf("  (rho = 1e-4: the same, %.2e)\n", worst4);
        bound("(r) the track does not depend on the source: u, ux between the two models, max |difference|", worstOwn, 1e-15);
    }

    // (s) An event that crosses the disc: rho = 0.01, u0 = 0.3 rho. At the peak the magnification is finite and
    // below the point source's (333), close to the disc-centre value sqrt(1 + 4/rho^2) = 200; the centroid shift
    // is the finite-source one, smaller than the point lens's thetaE u/(u^2 + 2) and of the opposite sign (the
    // images' flux-weighted centroid lies on the lens side of the source); and the library's centroid matches
    // the quadrature. The two seams agree: the measured centroid of evaluateModel is evaluateTrack's, bit for
    // bit. (rho = 0.01 and z = 0.3 are both nodes of the library's table, which makes the peak values agree
    // with the quadrature to rounding; the light curve and a second event at other values do not sit on one.)
    // The whole light curve is held to the quadrature.
    {
        const double rho = 0.01, u0 = 0.3 * rho, tetE = 0.6;
        fill(*ev, u0, 20.0, 1500.0, 0.25, 1.1, tetE);
        ev->s.rho = rho;  ev->sm = SourceModel::ESPL;
        ev->s.fb = {0.5, 0.8};  ev->s.fLens = {0.1, 0.05};  ev->s.blendOff = {{{0.3, -0.2}, {0.1, 0.4}}};
        const double t0 = ev->l.t0;
        TrackPoint fsT, psT;
        ModelPoint fsM;
        vbmSetSourceModel(SourceModel::ESPL);
        evaluateTrack(ev->s, ev->l, ev->as, 0, &t0, 1, true, &fsT);
        evaluateModel(ev->s, ev->l, ev->as, 0, &t0, 1, true, &fsM);
        vbmSetSourceModel(SourceModel::PSPL);
        evaluateTrack(ev->s, ev->l, ev->as, 0, &t0, 1, true, &psT);
        const double centre = std::sqrt(1.0 + 4.0 / (rho * rho));
        std::printf("  (peak: A = %.3f for the disc, %.3f for the point source, disc centre %.3f)\n", fsT.A, psT.A, centre);
        expect(std::isfinite(fsT.A) && fsT.A < psT.A, "(s) peak A of the disc is finite and below the point source's", fsT.A, "<", psT.A);
        expect(fsT.A > 0.9 * centre && fsT.A < centre, "(s) ... and between 0.9 and 1 times the disc-centre value", fsT.A, "<", centre);
        const DiscRef q = uniformDisc(fsT.u, rho);
        bound("(s) peak A vs the quadrature: relative difference", std::fabs(fsT.A / q.A - 1.0), 1e-4);

        // the shift of the images' centroid from the source, along the lens -> source axis, in units of thetaE
        const double defFs = std::hypot(fsT.def1c, fsT.def2c) / tetE, defPs = std::hypot(psT.def1c, psT.def2c) / tetE;
        const double defRef = q.h / q.A - fsT.u;
        const double along = (fsT.def1c * fsT.ux + fsT.def2c * fsT.uy) / (tetE * fsT.u);
        std::printf("  (peak centroid shift / thetaE along the source-lens axis: disc %.5e, quadrature %.5e, point lens %.5e)\n",
                    along, defRef, defPs);
        expect(defFs < 0.9 * defPs, "(s) peak centroid shift of the disc is smaller than the point lens's (ratio)", defFs / defPs, "<", 0.9);
        expect(along < 0.0, "(s) ... and reversed: the centroid lies on the lens side of the source (shift along lens -> source)", along, "<", 0.0);
        bound("(s) ... along the source-lens vector: ||shift| - |its component along it|| / |shift|", std::fabs(defFs - std::fabs(along)) / defFs, 1e-9);
        bound("(s) ... and matches the quadrature's: |difference| / (point lens's shift)", std::fabs(along - defRef) / defPs, 2e-3);
        expect(fsM.A == fsT.A && fsM.pos1c == fsT.pos1c && fsM.pos2c == fsT.pos2c,
               "(s) the seams agree: A and the measured centroid of evaluateModel and evaluateTrack, bitwise",
               fsM.A == fsT.A && fsM.pos1c == fsT.pos1c && fsM.pos2c == fsT.pos2c, "=", 1.0);

        // the whole light curve, through the disc and out the other side, against the quadrature
        const int n = 400;
        std::vector<double> tt(n);
        for (int k = 0; k < n; ++k) tt[k] = t0 - 3.0 * ev->l.tE + 6.0 * ev->l.tE * k / (n - 1);
        std::vector<TrackPoint> lc(n);
        vbmSetSourceModel(SourceModel::ESPL);
        evaluateTrack(ev->s, ev->l, ev->as, 0, tt.data(), n, true, lc.data());
        double worstA = 0.0, worstD = 0.0;
        for (int k = 0; k < n; ++k) {
            const DiscRef r = uniformDisc(lc[k].u, rho);
            worstA = std::fmax(worstA, std::fabs(lc[k].A / r.A - 1.0));
            const double dRef = r.h / r.A - lc[k].u, d = (lc[k].def1c * lc[k].ux + lc[k].def2c * lc[k].uy) / (tetE * lc[k].u);
            worstD = std::fmax(worstD, std::fabs(d - dRef) / (lc[k].u / (lc[k].u * lc[k].u + 2.0)));   // of the point lens's shift
        }
        bound("(s) light curve, 400 epochs over t0 +- 3 tE: A vs the quadrature, max relative difference", worstA, 1e-5);
        bound("(s) light curve: centroid shift vs the quadrature, max |difference| / (point lens's shift)", worstD, 6e-5);

        // a second event, with rho and u0/rho off the table's nodes: the interpolation error shows
        const double rho2 = 0.0137;
        fill(*ev, 0.3725 * rho2, 20.0, 1500.0, 0.25, 1.1, tetE);
        ev->s.rho = rho2;  ev->sm = SourceModel::ESPL;
        TrackPoint off;
        evaluateTrack(ev->s, ev->l, ev->as, 0, &t0, 1, true, &off);
        const DiscRef q2 = uniformDisc(off.u, rho2);
        const double along2 = (off.def1c * off.ux + off.def2c * off.uy) / (tetE * off.u);
        bound("(s) rho = 0.0137, u0 = 0.3725 rho (off the table's nodes): peak A vs the quadrature", std::fabs(off.A / q2.A - 1.0), 2e-5);
        bound("(s) ... and its centroid shift: |difference| / (point lens's shift)", std::fabs(along2 - (q2.h / q2.A - off.u)) / (off.u / (off.u * off.u + 2.0)), 2e-5);
    }

    // (t) Under the point-source model a source radius changes nothing, bit for bit: every field of the
    // library's adapter and of both seams is the same with s.rho = 0.01 as with s.rho = 0, A is magnifOf(u)
    // and the no-parallax A is magnifOf(uNoPlx). Under ESPL the same event differs.
    {
        int same = 0, total = 0, magnif = 0, nMagnif = 0, differ = 0;
        vbmSetSourceModel(SourceModel::PSPL);
        for (double u0 : {0.003, 0.05, 0.4}) {
            Ev a, b;
            fill(a, u0, 20.0, 1500.0, 0.25, 1.1, 0.6);  fill(b, u0, 20.0, 1500.0, 0.25, 1.1, 0.6);
            b.s.rho = 0.01;
            std::vector<double> tt(500);
            for (int k = 0; k < 500; ++k) tt[k] = 1440.0 + 0.2 * k;
            std::vector<ModelEpoch> ma(500), mb(500);
            std::vector<TrackPoint> ta(500), tb(500);
            std::vector<ModelPoint> pa(500), pb(500);
            for (int tele = 0; tele < 2; ++tele)
                for (bool astro : {false, true}) {
                    vbmEvaluate(a.m(), tele == 1, astro, tt.data(), 500, ma.data());
                    vbmEvaluate(b.m(), tele == 1, astro, tt.data(), 500, mb.data());
                    for (int k = 0; k < 500; ++k) {
                        total += 3;  same += (ma[k].A == mb[k].A) + (ma[k].u == mb[k].u) + (ma[k].ux == mb[k].ux);
                        if (astro) { total += 2; same += (ma[k].imgx == mb[k].imgx) + (ma[k].lensy == mb[k].lensy); }
                    }
                }
            for (int tele = 0; tele < 2; ++tele) {
                evaluateTrack(a.s, a.l, a.as, tele, tt.data(), 500, true, ta.data());
                evaluateTrack(b.s, b.l, b.as, tele, tt.data(), 500, true, tb.data());
                evaluateModel(a.s, a.l, a.as, tele, tt.data(), 500, true, pa.data());
                evaluateModel(b.s, b.l, b.as, tele, tt.data(), 500, true, pb.data());
                for (int k = 0; k < 500; ++k) {
                    total += 5;
                    same += (ta[k].A == tb[k].A) + (ta[k].ANoPlx == tb[k].ANoPlx) + (ta[k].pos1c == tb[k].pos1c)
                          + (pa[k].A == pb[k].A) + (pa[k].pos2c == pb[k].pos2c);
                    nMagnif += 3;
                    magnif += (tb[k].A == magnifOf(tb[k].u)) + (tb[k].ANoPlx == magnifOf(tb[k].uNoPlx)) + (pb[k].A == magnifOf(pb[k].u));
                }
            }
            // and the disc is different at the peak
            vbmSetSourceModel(SourceModel::ESPL);
            const double t0 = 1500.0;
            TrackPoint tfs, tps;
            evaluateTrack(b.s, b.l, b.as, 0, &t0, 1, true, &tfs);
            vbmSetSourceModel(SourceModel::PSPL);
            evaluateTrack(b.s, b.l, b.as, 0, &t0, 1, true, &tps);
            differ += (tfs.A != tps.A);
        }
        expect(same == total, "(t) point-source model: s.rho = 0.01 vs 0 changes no field of the adapter or the seams (bitwise)", same, "=", total);
        expect(magnif == nMagnif, "(t) ... A == magnifOf(u) and the no-parallax A == magnifOf(uNoPlx), bitwise", magnif, "=", nMagnif);
        expect(differ == 3, "(t) ... and under the disc the peak A is different, for the three events (of 3)", differ, "=", 3.0);
    }

    // (u) The magnification of the straight line without parallax (TrackPoint::ANoPlx) is the same function of
    // its own u: ESPLMag2(uNoPlx, rho) for the disc, bitwise, magnifOf(uNoPlx) for the point.
    {
        fill(*ev, 0.02, 25.0, 1500.0, 0.3, 1.1, 0.6);
        ev->s.rho = 0.01;
        std::vector<double> tt(300);
        for (int k = 0; k < 300; ++k) tt[k] = 1440.0 + 0.4 * k;
        std::vector<TrackPoint> tp(300);
        int fsSame = 0, psSame = 0;
        vbmSetSourceModel(SourceModel::ESPL);
        evaluateTrack(ev->s, ev->l, ev->as, 1, tt.data(), 300, false, tp.data());
        for (int k = 0; k < 300; ++k) fsSame += (tp[k].ANoPlx == lib->ESPLMag2(tp[k].uNoPlx, 0.01));
        vbmSetSourceModel(SourceModel::PSPL);
        evaluateTrack(ev->s, ev->l, ev->as, 1, tt.data(), 300, false, tp.data());
        for (int k = 0; k < 300; ++k) psSame += (tp[k].ANoPlx == magnifOf(tp[k].uNoPlx));
        expect(fsSame == 300, "(u) disc: ANoPlx == ESPLMag2(uNoPlx, rho), bitwise, of 300 epochs through the peak", fsSame, "=", 300.0);
        expect(psSame == 300, "(u) point: ANoPlx == magnifOf(uNoPlx), bitwise, of 300", psSame, "=", 300.0);
    }

    // (v) The edges. A disc of no size is the point source (the library's point-source function). A radius beyond
    // the library's table (100) is held at 99, not read out of bounds. A source at the centre of the lens
    // (u0 = 1e-9 rho) has the finite magnification sqrt(1 + 4/rho^2).
    {
        fill(*ev, 0.3, 25.0, 1500.0, 0.25, 1.1, 0.6);
        const double tt[3] = {1480.0, 1500.0, 1520.0};
        ModelEpoch pt[3], z0[3], big[3], cap[3];
        ev->sm = SourceModel::PSPL;  vbmEvaluate(ev->m(), false, true, tt, 3, pt);
        ev->sm = SourceModel::ESPL;  ev->s.rho = 0.0;   vbmEvaluate(ev->m(), false, true, tt, 3, z0);
        ev->s.rho = 1e3;   vbmEvaluate(ev->m(), false, true, tt, 3, big);
        ev->s.rho = 99.0;  vbmEvaluate(ev->m(), false, true, tt, 3, cap);
        int zero = 0, held = 0, fin = 0;
        for (int k = 0; k < 3; ++k) {
            zero += (z0[k].A == pt[k].A) + (z0[k].imgx == pt[k].imgx);
            held += (big[k].A == cap[k].A) + (big[k].imgx == cap[k].imgx);
            fin += std::isfinite(big[k].A) && std::fabs(big[k].A - 1.0) < 1e-3 && std::isfinite(big[k].imgx);
        }
        expect(zero == 6, "(v) ESPL with rho = 0 is the point source, bitwise (of 6)", zero, "=", 6.0);
        expect(held == 6 && fin == 3, "(v) rho = 1e3 is held at 99: same as 99, finite, A within 1e-3 of 1", held + fin, "=", 9.0);
        const double rho = 0.01;
        fill(*ev, 1e-9 * rho, 25.0, 1500.0, 1e-9, 1.1, 0.6);
        ev->s.rho = rho;  ev->sm = SourceModel::ESPL;
        const double t = 1500.0;
        ModelEpoch o;
        vbmEvaluate(ev->m(), false, false, &t, 1, &o);
        bound("(v) u0 = 1e-9 rho: A vs the centre of the disc sqrt(1 + 4/rho^2), relative difference", std::fabs(o.A / std::sqrt(1.0 + 4.0 / (rho * rho)) - 1.0), 1e-4);
    }
    vbmSetSourceModel(SourceModel::ESPL);

    std::printf("%s (%d failed)\n", nfail ? "FAILED" : "all held", nfail);
    return nfail;
}
