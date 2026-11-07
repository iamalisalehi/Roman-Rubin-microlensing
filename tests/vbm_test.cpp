// Unit test for the vendored VBMicrolensing library (external/VBMicrolensing/), against references that do
// not use it: tests/vbm_reference.h, written by tests/vbm_reference.py. Needs only the vendored files; run
// from the repo root. Exit status = the number of failed checks (0 = all held).
//
// Pinned: (a) finite-source point-lens magnification against direct integrals over the source disc, on a
// grid rho = 1e-3..0.1 by z = u/rho = 0.01..30, as the worst relative error of ESPLMag (the tabulated
// uniform disc) and of ESPLMag2 (Tol = 1e-3) for a uniform and a linearly limb-darkened (a1 = 0.5) disc.
// The source centre, z = 0 and 1e-3, is printed, not asserted: ESPLMag2 builds a limb-darkened disc from
// annuli starting at the point-source magnification of the centre, which diverges as u -> 0, and its
// error grows there (NaN at u = 0 exactly), whatever Tol.
// (b) PSPLMag is our point-source formula (src/events/lightcurve.cpp, src/fisher/fisher.cpp) bit for bit,
// so the two are interchangeable, and its image centroid lies at u + u/(u^2 + 2) from the lens.
// (c) the annual-parallax source track of PSPLLightCurveParallax (geocentric frame at t0, times in
// JD - 2450000 with the light-travel correction to heliocentric time) against astropy's ephemeris, at two
// bulge sightlines and two epochs over t0 +- 300 d. The 1e-4 Einstein-radius bound, against parallax terms
// of ~1.9, resolves the light-travel correction (~1.7e-4). What remains (~5e-5) is mostly linear in
// t - t0: the library takes the Earth's velocity at t0 from its daily table, and a linear drift is
// absorbed by tE and the direction of motion.
#include "parameters.h"
#include "VBMicrolensingLibrary.h"
#include "vbm_reference.h"
#include <cstdio>
#include <cmath>
#include <memory>

static int nfail = 0;
// One check: the measured value, the exact value or bound it is held to, and ok/FAIL.
static void expect(bool ok, const char* what, double got, const char* rel, double want) {
    std::printf("%-68s got %11.4e  want %s %9.2e  %s\n", what, got, rel, want, ok ? "ok" : "FAIL");
    if (!ok) ++nfail;
}

static bool readable(const char* path) {
    std::FILE* f = std::fopen(path, "rb");
    if (!f) return false;
    std::fclose(f);
    return true;
}

// Relative error, infinite when the value is not finite, so a NaN can never pass as a small error.
static double relerr(double got, double want) {
    return std::isfinite(got) ? std::fabs(got / want - 1.0) : INFINITY;
}

// The relative error as text for the printed-only rows; a non-finite value is shown as itself.
static void errText(char* buf, std::size_t n, double got, double want) {
    if (std::isfinite(got)) std::snprintf(buf, n, "%.1e", relerr(got, want));
    else std::snprintf(buf, n, "%g", got);
}

int main() {
    using namespace vbmref;
    // The library only prints a warning for a missing table and then computes with an empty one.
    for (const char* p : {PATH_VBM_ESPL_TABLE, PATH_VBM_SUN_TABLE})
        if (!readable(p)) { std::printf("cannot read %s (run from the repository root)\n", p); return 1; }

    // Half a megabyte of tables inside the object: keep it off the stack.
    auto V = std::make_unique<VBMicrolensing>();
    V->LoadESPLTable(PATH_VBM_ESPL_TABLE);
    char sunPath[1024];
    std::snprintf(sunPath, sizeof sunPath, "%s", PATH_VBM_SUN_TABLE);   // LoadSunTable takes a char*
    V->LoadSunTable(sunPath);
    expect(!V->ESPLoff && V->suntable, "ESPL and Sun tables loaded", !V->ESPLoff && V->suntable, "=", 1.0);

    // (a) finite source. ESPLMag is the uniform disc by construction; ESPLMag2 integrates the profile
    // set by a1 over annuli of ESPLMag, to an absolute tolerance Tol on the magnification.
    V->Tol = 1e-3;
    const char* name[3] = {"ESPLMag  uniform disc", "ESPLMag2 uniform disc", "ESPLMag2 limb-darkened a1=0.5"};
    const double bound[3] = {1e-4, 1e-4, 1e-3};
    for (int c = 0; c < 3; ++c) {
        const int ia = (c == 2) ? 1 : 0;
        V->a1 = A1[ia];
        double worst = 0.0, wrho = 0.0, wz = 0.0;
        for (int ir = 0; ir < N_RHO; ++ir)
            for (int iz = 0; iz < N_Z; ++iz) {
                const double u = Z[iz] * RHO[ir];
                const double got = (c == 0) ? V->ESPLMag(u, RHO[ir]) : V->ESPLMag2(u, RHO[ir]);
                const double e = relerr(got, ESPL[ia][ir][iz]);
                if (e > worst) { worst = e; wrho = RHO[ir]; wz = Z[iz]; }
            }
        char what[96];
        std::snprintf(what, sizeof what, "%s: worst rel. err. (rho %.0e, z %g)", name[c], wrho, wz);
        expect(worst < bound[c], what, worst, "<", bound[c]);
    }
    // The source centre: printed so the behaviour is on record, not asserted.
    std::printf("source centre, not asserted:   a1    rho      z     reference  rel.err: ESPLMag  ESPLMag2\n");
    for (int ia = 0; ia < N_A1; ++ia) {
        V->a1 = A1[ia];
        for (int ir = 0; ir < N_RHO; ++ir)
            for (int iz = 0; iz < N_ZC; ++iz) {
                const double u = Z_CENTRE[iz] * RHO[ir], ref = ESPL_CENTRE[ia][ir][iz];
                char e1[16] = "-", e2[16];
                if (ia == 0) errText(e1, sizeof e1, V->ESPLMag(u, RHO[ir]), ref);
                errText(e2, sizeof e2, V->ESPLMag2(u, RHO[ir]), ref);
                std::printf("%33.1f %6.0e %6g %13.6f %16s %9s\n", A1[ia], RHO[ir], Z_CENTRE[iz], ref, e1, e2);
            }
    }
    V->a1 = 0.0;

    // (b) point source, against the expression our light curves and Fisher matrices use.
    const int NU = 50;
    int same = 0, centroid = 0;
    V->astrometry = true;                       // PSPLMag then also sets astrox1
    for (int i = 0; i < NU; ++i) {
        const double u = std::pow(10.0, -4.0 + 5.0 * i / (NU - 1));     // 1e-4 .. 10
        const double ours = (u * u + 2.0) / std::sqrt(u * u * (u * u + 4.0));
        if (V->PSPLMag(u) == ours) ++same;
        if (std::fabs(V->astrox1 / (u + u / (u * u + 2.0)) - 1.0) < 1e-14) ++centroid;
    }
    V->astrometry = false;
    expect(same == NU, "PSPLMag == (u^2+2)/sqrt(u^2(u^2+4)) bit for bit, of 50 u in 1e-4..10", same, "=", NU);
    expect(centroid == NU, "PSPLMag centroid astrox1 = u + u/(u^2+2) to 1e-14, of 50 u", centroid, "=", NU);

    // (c) annual parallax: u0, ln tE, t0, piE_N, piE_E; times in JD - 2450000, not heliocentric.
    V->t_in_HJD = false;
    for (int is = 0; is < N_SIGHT; ++is) {
        char coord[64];
        std::snprintf(coord, sizeof coord, "%s", COORDS[is]);
        V->SetObjectCoordinates(coord);
        if (!V->AreCoordinatesSet()) { std::printf("SetObjectCoordinates refused \"%s\"\n", COORDS[is]); return 1; }
        for (int ie = 0; ie < N_EVENT; ++ie) {
            double pr[5] = {U0, std::log(TE), T0[ie], PIEN, PIEE};
            double t[N_T], mag[N_T], y1[N_T], y2[N_T];
            for (int k = 0; k < N_T; ++k) t[k] = T[ie][k];
            V->PSPLLightCurveParallax(pr, t, mag, y1, y2, N_T);
            double worst = 0.0;
            for (int k = 0; k < N_T; ++k) {
                const double d = std::hypot(y1[k] - Y1[is][ie][k], y2[k] - Y2[is][ie][k]);
                worst = std::isfinite(d) ? std::fmax(worst, d) : INFINITY;
            }
            char what[96];
            std::snprintf(what, sizeof what, "parallax (l,b)=(%g,%g) t0=day %g: max|y-y_ref| [term to %.2f]",
                          GAL_L[is], GAL_B[is], T0_DAY[ie], PAR_MAX[is][ie]);
            expect(worst < 1e-4 && V->parallaxextrapolation == 0, what, worst, "<", 1e-4);
        }
    }

    std::printf("%s (%d failed)\n", nfail ? "FAILED" : "all held", nfail);
    return nfail;
}
