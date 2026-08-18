// tests/fisher_fixture.cpp
//
// A deliberately simple, self-contained Fisher-matrix regression fixture.
//
// WHY THIS EXISTS
// ---------------
// Every numerical acceptance test in Phase C of JOINT_FIT_REFACTOR_PLAN.md has been throttled
// by the live Monte Carlo's detection efficiency: the production binary needs thousands of star
// draws (measured: 1 detected event per 886, 2600, and 10735 draws on three separate runs) and
// tens of minutes to yield a single Fisher-scored event. That is far too slow to check that a
// change to FisherM moved sigma in the expected direction, and far too few events to see a
// distribution.
//
// This fixture bypasses the Monte Carlo entirely. It hand-builds a small set of microlensing
// events with known parameters, synthesises their light curves on representative Roman and Rubin
// cadences, runs the real FisherM and ErrorCal on them, and prints the resulting sigmas in a
// stable, diffable format. It runs in seconds.
//
// DELIBERATELY SELF-CONTAINED: this fixture reads NO data files. The Baseline/, CMD/, and
// extinction inputs the production binary needs are gitignored and absent after a fresh clone,
// so depending on them would make the test unrunnable for anyone but the original author.
// Source magnitudes, blend fractions, cadences and photometric errors are therefore hardcoded
// representative values, not draws from the real population.
//
// WHAT IT IS AND IS NOT
// ---------------------
// IS:     a regression harness. Run it before and after a change to FisherM; diff the output.
//         Sigmas should move in the direction the change predicts, and nowhere else.
// IS NOT: a validation of absolute precision forecasts. The events are hand-picked, not
//         population-weighted, and the error model is a flat stand-in rather than errlsstM /
//         errRomanM. Absolute sigma values here are indicative, not publishable.
//
// USAGE
//   make fishertest && ./fishertest
//
// Note for whoever reads the git history: an earlier throwaway version of this idea (during
// Step C1) produced wildly degenerate sigmas and was misdiagnosed at the time as a badly
// conditioned fixture. It was not -- it was correctly reporting the Fisher accumulation bug
// later fixed in commit d5c8867 (FisherM overwrote instead of summing over epochs, leaving a
// rank-1 matrix). The fixture idea was sound; the code under test was broken.

#include "Bulge.h"

namespace {

// ---------------------------------------------------------------------------------------------
// Sight line. Fixed to the same test field the production runs have been using, so that the
// parallax geometry (deltao, FI) matches what the real code would compute there.
// ---------------------------------------------------------------------------------------------
constexpr double kLon = 0.5;   // Galactic longitude [deg]
constexpr double kLat = -1.0;  // Galactic latitude  [deg]

// ---------------------------------------------------------------------------------------------
// Cadence model. Simplified stand-ins for the two surveys' observing patterns -- enough to
// reproduce the structural feature that matters (Roman dense-but-seasonal, Rubin sparse-but-
// continuous), without needing the real visit lists.
// ---------------------------------------------------------------------------------------------

// Roman GBTDS: 6 high-cadence seasons of ~72 days, front- and back-loaded around a multi-year
// mid-mission gap. Real cadence in F146 is ~12-15 min; sampled here at 4 h to keep the fixture
// fast. That reduces the absolute information content but preserves the seasonal structure,
// which is what the gap-filling science actually turns on.
constexpr double kRomanSeasonStarts[] = {0.0, 180.0, 1600.0, 1780.0, 3300.0, 3480.0};
constexpr double kRomanSeasonLength   = 72.0;      // days
constexpr double kRomanCadence        = 4.0 / 24.0; // days (4 hours)

// Rubin LSST: ~3-day cadence, but the bulge is only observable for part of each year.
constexpr double kRubinCadence        = 3.0;   // days
constexpr double kRubinVisibleDays    = 250.0; // observable window per year
constexpr int    kRubinYears          = 10;

// Representative 1-sigma per-epoch uncertainties. Flat stand-ins for errlsstM / errRomanM,
// which need file-loaded survey data this fixture deliberately avoids.
constexpr double kErrMagRubin  = 0.020; // [mag]
constexpr double kErrMagRoman  = 0.005; // [mag] -- space-based, no atmosphere
constexpr double kErrAstRubin  = 0.500; // [mas]
// NOTE: Roman has no real astrometric error model yet (see OPEN_ITEMS.md). This is a
// placeholder standing in for one, consistent with what Bulge_LSST.cpp currently does.
constexpr double kErrAstRoman  = 0.050; // [mas]

// ---------------------------------------------------------------------------------------------
// The test events. Chosen to span tE and, at each tE, to contrast an event peaking inside a
// Roman season against one peaking in a Roman gap -- the contrast the whole joint-fit science
// case rests on.
// ---------------------------------------------------------------------------------------------
struct FixtureEvent {
    const char* name;
    double tE;    // Einstein crossing time [days]
    double t0;    // time of peak magnification [days since survey start]
    double u0;    // impact parameter [Einstein radii]
    double piE;   // microlensing parallax []
};

// The t0 values carry a deliberate +0.37 d offset so that no epoch of either cadence grid
// lands exactly on t0. That is not cosmetic: at timh == t0 the mus1/mus2 terms in the modelled
// source position vanish identically, perturbing them changes nothing, and FisherM's
// astrometric CHECK aborts the program. Real t0 is drawn from a continuous distribution so the
// coincidence is measure-zero there, but a regular test grid hits it easily.
const FixtureEvent kEvents[] = {
    {"short_inseason",  5.0,   36.37, 0.30, 0.30},
    {"short_ingap",     5.0,  900.37, 0.30, 0.30},
    {"mid_inseason",   25.0, 1636.37, 0.30, 0.15},
    {"mid_ingap",      25.0, 2500.37, 0.30, 0.15},
    {"long_inseason", 100.0, 1636.37, 0.30, 0.10},
    {"long_ingap",    100.0, 2500.37, 0.30, 0.10},
    {"verylong",      900.0, 1800.37, 0.30, 0.05},
};

// Parameters held fixed across all events, so that differences in the printed sigmas are
// attributable to tE / t0 / piE and the cadence, not to incidental source-star differences.
constexpr double kXi        = 0.70;   // trajectory angle [radian]; kept away from 0 and pi/2,
                                      // where ErrorCal's tan(xi) and 1/tan(xi) terms blow up
constexpr double kTetE      = 0.50;   // angular Einstein radius [mas]
constexpr double kDs        = 8.00;   // source distance [kpc]
constexpr double kDl        = 4.00;   // lens distance [kpc]
constexpr double kMbsRubin  = 20.00;  // Rubin baseline magnitude (source+blend)
constexpr double kMbsRoman  = 19.00;  // Roman F146 baseline magnitude
constexpr double kFbRubin   = 0.50;   // Rubin blend fraction
constexpr double kFbRoman   = 0.80;   // Roman blend fraction -- higher, because F146's much
                                      // smaller PSF deblends what Rubin's seeing cannot
constexpr double kMus1      = +0.008; // source proper motion [mas/day] (~2.9 mas/yr)
constexpr double kMus2      = -0.005;
constexpr double kMul1      = +0.003; // lens proper motion [mas/day]
constexpr double kMul2      = +0.002;

// Fill in the geometry and stellar quantities that lightcurve() and ErrorCal() read but that
// are not what this fixture varies.
void setupStatic(source& s, lens& l)
{
    s.lon = kLon;
    s.lat = kLat;
    s.TET = (360.0 - s.lon) / RAa;  // same convention as Bulge_LSST.cpp main()
    s.FI  = s.lat / RAa;

    // Reproduce Lensing.cpp's deltao convention for this sight line.
    double tetd = s.TET;
    if (s.TET > pi) tetd = s.TET - 2.0 * pi;
    l.deltao = pi - std::fabs(tetd);
    if (tetd < 0.0) l.deltao = -1.0 * l.deltao;

    s.Ds  = kDs;
    l.Dl  = kDl;
    s.xi  = kXi;
    l.tetE = kTetE;

    s.mbs[0] = kMbsRubin;  s.fb[0] = kFbRubin;
    s.mbs[1] = kMbsRoman;  s.fb[1] = kFbRoman;

    s.mus1 = kMus1;  s.mus2 = kMus2;
    l.mul1 = kMul1;  l.mul2 = kMul2;
    s.mus  = std::sqrt(s.mus1 * s.mus1 + s.mus2 * s.mus2);
    l.mul  = std::sqrt(l.mul1 * l.mul1 + l.mul2 * l.mul2);
    l.murel = std::sqrt((s.mus1 - l.mul1) * (s.mus1 - l.mul1)
                      + (s.mus2 - l.mul2) * (s.mus2 - l.mul2));
}

// True (noiseless) model magnitude for telescope tt at the current lightcurve() state.
double modelMag(const source& s, int tt, double Astar)
{
    return s.mbs[tt] - 2.5 * std::log10(Astar * s.fb[tt] + 1.0 - s.fb[tt]);
}

// Append one epoch to the light-curve buffers. The stored magnitude and astrometric position
// are the NOISELESS model values at the true parameters -- which is exactly what a Fisher
// forecast wants. FisherM differences the perturbed model against these, and for the one-sided
// stencil (sig2, used for tE and piE) the stored value does not cancel, so it must be the true
// model rather than a noise realisation.
void appendEpoch(source& s, lens& l, astromet& as, double tim, int tele, int& ndw)
{
    lightcurve(s, l, as, tim);
    s.Astar = (s.ut * s.ut + 2.0) / std::sqrt(s.ut * s.ut * (s.ut * s.ut + 4.0));

    l.timn[ndw] = tim;
    l.magn[ndw] = modelMag(s, tele, s.Astar);
    l.errm[ndw] = (tele == 0) ? kErrMagRubin : kErrMagRoman;
    l.soux[ndw] = s.pos1c;
    l.souy[ndw] = s.pos2c;
    l.erra[ndw] = (tele == 0) ? kErrAstRubin : kErrAstRoman;
    l.tele[ndw] = tele;
    ++ndw;
}

// Half-width of the fit window around t0, in days. Epochs outside it are dropped.
//
// This is partly physical and partly defensive.
//
// Physical: a real light-curve fit uses a window scaled to the event duration. Epochs at
// |t - t0| >> tE carry only baseline information, and no analyst fits a 5-day event using ten
// years of data either side of it.
//
// Defensive: the astrometric branch of FisherM asserts
//     CHECK(!(s.pos1c == l.soux[i] && s.pos2c == l.souy[i]));
// which fires when a parameter perturbation changes the modelled source position by less than
// a double can represent. That happens for epochs enormously far from peak, where the
// astrometric deflection def1c/def2c has decayed as 1/u^2 into the floating-point noise floor,
// so a piE perturbation returns bit-identical coordinates. Without this window, the short-tE
// events below trip that assert and abort the fixture. This is a latent hazard in the
// production code too (see OPEN_ITEMS.md); the window sidesteps it rather than papering over
// it, and the fixture is deliberately left able to surface it again if the window is widened.
double fitHalfWindow(double tE)
{
    double w = 20.0 * tE;
    if (w < 100.0)  w = 100.0;   // floor: always keep some baseline either side
    if (w > Tobs)   w = Tobs;
    return w;
}

// Build the full merged light curve for one event. Returns total epoch count; reports the
// per-survey split through nL / nR.
int buildLightCurve(source& s, lens& l, astromet& as, int& nL, int& nR)
{
    int ndw = 0;
    nL = 0;
    nR = 0;

    const double halfWin = fitHalfWindow(l.tE);
    const double tMin    = l.t0 - halfWin;
    const double tMax    = l.t0 + halfWin;

    // Roman: dense sampling inside each season.
    for (double start : kRomanSeasonStarts) {
        for (double t = start; t <= start + kRomanSeasonLength; t += kRomanCadence) {
            if (t < 0.0 || t > Tobs) continue;
            if (t < tMin || t > tMax) continue;
            appendEpoch(s, l, as, t, 1, ndw);
            ++nR;
        }
    }

    // Rubin: sparse sampling across the visible window of each year.
    for (int y = 0; y < kRubinYears; ++y) {
        const double yearStart = y * year;
        for (double t = yearStart; t <= yearStart + kRubinVisibleDays; t += kRubinCadence) {
            if (t < 0.0 || t > Tobs) continue;
            if (t < tMin || t > tMax) continue;
            appendEpoch(s, l, as, t, 0, ndw);
            ++nL;
        }
    }
    return ndw;
}

} // namespace

namespace {

// Recompute F[SRUBIN] + F[SROMAN] and compare against F[SJOINT] element by element.
//
// This is the sharpest available check that the partitioning is correct. Fisher information is
// a sum of independent per-epoch terms, and every epoch feeds the joint matrix and exactly one
// single-survey matrix, so the identity is exact up to floating-point summation order.
bool checkAdditivity(const covarian& co, const char* evName, int dim, bool photometric)
{
    bool ok = true;
    double worst = 0.0;
    for (int j = 0; j < dim; ++j) {
        for (int k = 0; k < dim; ++k) {
            const double joint = photometric
                ? gsl_matrix_get(co.inputA[SJOINT].get(), j, k)
                : gsl_matrix_get(co.inputB[SJOINT].get(), j, k);
            const double parts = photometric
                ? gsl_matrix_get(co.inputA[SRUBIN].get(), j, k)
                    + gsl_matrix_get(co.inputA[SROMAN].get(), j, k)
                : gsl_matrix_get(co.inputB[SRUBIN].get(), j, k)
                    + gsl_matrix_get(co.inputB[SROMAN].get(), j, k);
            const double scale = std::fabs(joint) + std::fabs(parts) + 1e-300;
            const double rel   = std::fabs(joint - parts) / scale;
            if (rel > worst) worst = rel;
            if (rel > 1e-9) ok = false;
        }
    }
    if (!ok) {
        std::cerr << "FAIL [" << evName << "] " << (photometric ? "photometric" : "astrometric")
                  << " matrix: F[JOINT] != F[RUBIN] + F[ROMAN]; worst relative mismatch "
                  << worst << "\n";
    }
    return ok;
}

// Joint forecast versus single-survey forecast, for the SHARED parameters only.
//
// Restricted to the geometric parameters {u0, tE, piE, xi, t0} on purpose. The naive statement
// "adding data can never increase Fisher information" applies to a fixed parameter set. Here the
// partitions deliberately have different parameter sets: the joint fit marginalizes over BOTH
// telescopes' flux parameters, while a single-survey fit carries only its own. For a
// telescope-specific flux parameter that makes the comparison meaningless -- sigma(fb1) is
// legitimately larger in the joint fit, because the joint fit must also solve for Rubin's flux
// scale while Rubin's epochs contribute nothing to fb1 itself.
//
// Even on the shared parameters this is an empirical expectation rather than a theorem: by the
// Schur complement, the joint fit's effective information for the shared block is
// F_rubin + F_roman(shared) minus a positive semi-definite penalty for marginalizing over the
// extra flux parameters. In practice the added data dominates by a wide margin, so a violation
// here is a strong signal that something is wrong and worth investigating -- which is why it is
// still a hard failure.
bool checkJointNoWorse(const covarian& co, const char* evName,
                       const char* const* pnames, int dim, bool photometric)
{
    bool ok = true;
    for (int q = 1; q < NSURV; ++q) {
        const bool okJoint = photometric ? co.okA[SJOINT] : co.okB[SJOINT];
        const bool okPart  = photometric ? co.okA[q]      : co.okB[q];
        if (!okJoint || !okPart) continue;  // nothing to compare against
        for (int k = 0; k < dim; ++k) {
            // Shared geometric parameters only (photometric indices 0,1,3,4,5). Index 2 is fb0,
            // 6 mbs0, 7 fb1, 8 mbs1 -- all telescope-specific. The astrometric set has no
            // telescope-specific parameters, so all of it is comparable.
            if (photometric && (k == 2 || k >= 6)) continue;
            const double sJoint = photometric ? co.Era[SJOINT][k] : co.Erb[SJOINT][k];
            const double sPart  = photometric ? co.Era[q][k]      : co.Erb[q][k];
            // A negative sigma means the parameter is not in this partition's active subset
            // (Rubin's matrix says nothing about Roman's flux scale). Nothing to compare.
            if (sJoint < 0.0 || sPart < 0.0) continue;
            // 1e-9 relative slack absorbs floating-point noise in the inversion.
            if (sJoint > sPart * (1.0 + 1e-9)) {
                std::cerr << "FAIL [" << evName << "] sigma(" << pnames[k] << ") joint="
                          << sJoint << " > " << (q == SRUBIN ? "rubin=" : "roman=") << sPart
                          << "  -- adding data made the forecast worse\n";
                ok = false;
            }
        }
    }
    return ok;
}

// Step C1's acceptance criterion, checked on every event and every survey partition.
//
// The pre-C1 parameter set {u0, tE, fb, piE, xi} is exactly the leading 5x5 submatrix of the
// current 6-parameter information matrix -- dropping t0's row and column is the same thing as
// asserting perfect knowledge of when the peak occurred. So both numbers come from the SAME
// accumulated information on the SAME event: a perfectly paired comparison, no second run and
// no rebuild.
//
// Cramer-Rao: marginalizing over a genuinely free parameter can only loosen the bound on the
// others. sigma(tE) with t0 free must therefore be >= sigma(tE) with t0 held fixed. A decrease
// anywhere means the extra parameter is somehow adding information, which is impossible.
//
// This works only because invert_matrix operates on a scratch copy and leaves inputA intact.
bool checkT0Marginalization(const covarian& co, const char* evName, bool verbose)
{
    if (Nx < 6) return true;  // nothing to compare against
    bool ok = true;

    for (int q = 0; q < NSURV; ++q) {
        if (!co.okA[q]) continue;

        // Pre-C1 parameter set {u0, tE, fb, piE, xi} = indices 0-4, intersected with this
        // partition's active subset (Roman's matrix has no fb0 at index 2).
        std::vector<int> sub;
        for (int k : activePhotParams(q, co.nepochA[SRUBIN], co.nepochA[SROMAN])) if (k <= 4) sub.push_back(k);
        const int nsub = static_cast<int>(sub.size());

        gsl_matrix      *F5 = gsl_matrix_alloc(nsub, nsub);
        gsl_matrix      *I5 = gsl_matrix_alloc(nsub, nsub);
        gsl_permutation *p5 = gsl_permutation_alloc(nsub);
        for (int a = 0; a < nsub; ++a)
            for (int b = 0; b < nsub; ++b)
                gsl_matrix_set(F5, a, b, gsl_matrix_get(co.inputA[q].get(), sub[a], sub[b]));

        int s5;
        gsl_linalg_LU_decomp(F5, p5, &s5);
        const double det5 = gsl_linalg_LU_det(F5, s5);

        if (det5 != 0.0 && std::isfinite(det5)) {
            gsl_linalg_LU_invert(F5, p5, I5);
            const double v = gsl_matrix_get(I5, 1, 1); //index 1 is tE in every subset
            if (std::isfinite(v) && v >= 0.0) {
                const double sigFixed = std::sqrt(v);
                const double sigFree  = co.Era[q][1];
                const char*  qn = (q == SJOINT) ? "joint" : (q == SRUBIN ? "rubin" : "roman");
                if (sigFree + 1e-12 < sigFixed) {
                    std::cerr << "FAIL [" << evName << "/" << qn << "] sigma(tE) DECREASED when "
                              << "t0 was added: fixed=" << sigFixed << " free=" << sigFree << "\n";
                    ok = false;
                } else if (verbose) {
                    std::cout << "#   t0-marginalization " << std::setw(5) << qn
                              << ": sigma(tE) fixed=" << std::scientific << std::setprecision(4)
                              << sigFixed << "  free=" << sigFree
                              << "  ratio=" << std::fixed << std::setprecision(4)
                              << sigFree / sigFixed << "\n";
                }
            }
        }
        gsl_permutation_free(p5);
        gsl_matrix_free(I5);
        gsl_matrix_free(F5);
    }
    return ok;
}

// The direct acceptance test for Step C4: F * F^-1 must be the identity.
//
// This validates the entire normalize / invert / rescale round trip in one shot -- if the
// diagonal scaling were applied inconsistently, or undone in the wrong order, the product would
// not come back as I even though every individual step looked reasonable. It only works because
// invert_matrix operates on scratch copies and leaves the accumulated information matrix intact.
//
// Tolerance is scaled by the condition number: losing about log10(cond) digits is expected and
// unavoidable, so a fixed tolerance would either be vacuous for well-conditioned matrices or
// spuriously fail for legitimately degenerate ones.
bool checkInverseRoundTrip(const covarian& co, const char* evName, int dim, bool photometric)
{
    bool ok = true;
    for (int q = 0; q < NSURV; ++q) {
        const bool valid = photometric ? co.okA[q] : co.okB[q];
        if (!valid) continue;

        const gsl_matrix* F  = photometric ? co.inputA[q].get() : co.inputB[q].get();
        const gsl_matrix* Fi = photometric ? co.inverA[q].get() : co.inverB[q].get();
        const double cond    = photometric ? co.condA[q] : co.condB[q];

        // Only the active submatrix was inverted, so only it can round-trip to the identity.
        static const std::vector<int> kAllAst = {0, 1, 2, 3};
        const std::vector<int> actP = photometric
            ? activePhotParams(q, co.nepochA[SRUBIN], co.nepochA[SROMAN]) : std::vector<int>();
        const std::vector<int>& act = photometric ? actP : kAllAst;
        (void)dim;
        const int n = static_cast<int>(act.size());

        // Measure the residual in the NORMALIZED basis, which is the one the inversion actually
        // worked in. Evaluating F * F^-1 in raw units means multiplying entries that span many
        // orders of magnitude, and the cancellation in that product swamps the answer for reasons
        // that have nothing to do with whether the inverse is correct. With
        // s_i = 1/sqrt(F_ii), Ftilde = D F D and Ftilde^-1 = D^-1 F^-1 D^-1, every term of the
        // product is O(1) and the residual measures what we actually care about.
        std::vector<double> s(n);
        for (int a = 0; a < n; ++a) {
            const double d = gsl_matrix_get(F, act[a], act[a]);
            s[a] = (d > 0.0) ? 1.0 / std::sqrt(d) : 0.0;
        }
        double worst = 0.0;
        for (int a = 0; a < n; ++a) {
            for (int b = 0; b < n; ++b) {
                double sum = 0.0;
                for (int c = 0; c < n; ++c) {
                    const double Ft  = gsl_matrix_get(F,  act[a], act[c]) * s[a] * s[c];
                    const double Fti = (s[c] > 0.0 && s[b] > 0.0)
                        ? gsl_matrix_get(Fi, act[c], act[b]) / (s[c] * s[b]) : 0.0;
                    sum += Ft * Fti;
                }
                const double target = (a == b) ? 1.0 : 0.0;
                const double err = std::fabs(sum - target);
                if (err > worst) worst = err;
            }
        }
        // Losing roughly log10(cond) digits from a ~1e-16 base is expected and unavoidable.
        const double tol = std::max(1e-9, 1e-14 * (cond > 0.0 ? cond : 1.0));
        if (worst > tol) {
            const char* qn = (q == SJOINT) ? "joint" : (q == SRUBIN ? "rubin" : "roman");
            std::cerr << "FAIL [" << evName << "/" << qn << "] "
                      << (photometric ? "photometric" : "astrometric")
                      << " F*Finv deviates from identity by " << worst
                      << " (tol " << tol << ", cond " << cond << ")\n";
            ok = false;
        }
    }
    return ok;
}

const char* kPhotNames[] = {"u0", "tE", "fb0", "piE", "xi", "t0", "mbs0", "fb1", "mbs1"};
const char* kAstNames[]  = {"tetE", "mus1", "mus2", "piE"};

void printRow(const char* label, int nep, int ok, double cond,
              const std::vector<double>& era, const std::vector<double>& erb)
{
    std::cout << std::left << std::setw(10) << label
              << std::right << std::setw(7) << nep
              << std::setw(5) << ok
              << std::scientific << std::setprecision(2) << std::setw(11) << cond;
    std::cout << std::scientific << std::setprecision(3);
    for (int k = 0; k < Nx; ++k) std::cout << std::setw(12) << era[k];
    std::cout << std::setw(12) << erb[0] << "\n";
}

} // namespace

int main()
{
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();

    int failures = 0;

    std::cout << "# Fisher-matrix fixture -- synthetic events, no data files required\n"
              << "# Nx=" << Nx << " (photometric)  Ny=" << Ny << " (astrometric)\n"
              << "# One block per event, one row per survey partition (Step C5).\n"
              << "# ok=1 usable inverse, ok=0 not characterizable (sigmas print as -1).\n"
              << "# Asserted: F[joint] == F[rubin] + F[roman] exactly, sigma_joint <= both, and\n"
              << "# sigma(tE) never decreases when t0 is marginalized over (Step C1 acceptance).\n";

    for (const auto& ev : kEvents) {
        setupStatic(*s, *l);
        l->tE  = ev.tE;
        l->t0  = ev.t0;
        l->u0  = ev.u0;
        l->piE = ev.piE;

        int nL = 0, nR = 0;
        const int ndw = buildLightCurve(*s, *l, *as, nL, nR);

        // buildLightCurve left lightcurve() evaluated at the final epoch; FisherM perturbs from
        // whatever the structs currently hold, so restore the true parameter point first.
        l->tE  = ev.tE;
        l->t0  = ev.t0;
        l->u0  = ev.u0;
        l->piE = ev.piE;
        s->xi  = kXi;
        s->fb[0] = kFbRubin;
        s->fb[1] = kFbRoman;

        // See OPEN_ITEMS.md: an epoch exactly at t0 makes the mus1/mus2 perturbations no-ops and
        // trips FisherM's astrometric CHECK. The t0 values are offset off-grid to avoid it.
        for (int i = 0; i < ndw; ++i) {
            if (l->timn[i] == l->t0) {
                std::cerr << "fixture error: event '" << ev.name << "' has an epoch exactly at "
                          << "t0=" << l->t0 << "; nudge t0 off the cadence grid.\n";
                return 1;
            }
        }

        FisherM(*s, *l, *as, *co, ndw);
        ErrorCal(*co, *l, *s);

        std::cout << "\n# " << ev.name << "  tE=" << std::fixed << std::setprecision(1) << ev.tE
                  << "  t0=" << ev.t0 << "  ndw=" << ndw
                  << "  (rubin=" << nL << ", roman=" << nR << ")\n";
        std::cout << std::left << std::setw(10) << "# survey"
                  << std::right << std::setw(7) << "nep" << std::setw(5) << "ok" << std::setw(11) << "cond";
        for (int k = 0; k < Nx; ++k) std::cout << std::setw(12) << kPhotNames[k];
        std::cout << std::setw(12) << kAstNames[0] << "\n";

        printRow("joint", co->nepochA[SJOINT], co->okA[SJOINT], co->condA[SJOINT], co->Era[SJOINT], co->Erb[SJOINT]);
        printRow("rubin", co->nepochA[SRUBIN], co->okA[SRUBIN], co->condA[SRUBIN], co->Era[SRUBIN], co->Erb[SRUBIN]);
        printRow("roman", co->nepochA[SROMAN], co->okA[SROMAN], co->condA[SROMAN], co->Era[SROMAN], co->Erb[SROMAN]);

        {
            static const char* kSyn[] = {"none", "both-alone", "rubin-only-alone",
                                         "roman-only-alone", "joint-only-RESCUE"};
            const int sc = synergyClass(*co);
            std::cout << "#   synergy: " << kSyn[sc];
            // Quantify what the joint fit bought over the best either survey managed alone.
            double best = -1.0;
            if (co->okA[SRUBIN]) best = co->Era[SRUBIN][1];
            if (co->okA[SROMAN] && (best < 0.0 || co->Era[SROMAN][1] < best))
                best = co->Era[SROMAN][1];
            if (co->okA[SJOINT] && best > 0.0) {
                std::cout << "   sigma(tE) joint/best-single = "
                          << std::fixed << std::setprecision(4)
                          << co->Era[SJOINT][1] / best;
            } else if (sc == SYN_JOINT_ONLY) {
                std::cout << "   (neither survey alone could characterize this event)";
            }
            std::cout << "\n";
        }

        // ---- assertions ----
        if (!checkAdditivity(*co, ev.name, Nx, true))  ++failures;
        if (!checkAdditivity(*co, ev.name, Ny, false)) ++failures;
        if (!checkJointNoWorse(*co, ev.name, kPhotNames, Nx, true))  ++failures;
        if (!checkJointNoWorse(*co, ev.name, kAstNames,  Ny, false)) ++failures;
        if (!checkT0Marginalization(*co, ev.name, true)) ++failures;
        if (!checkInverseRoundTrip(*co, ev.name, Nx, true))  ++failures;
        if (!checkInverseRoundTrip(*co, ev.name, Ny, false)) ++failures;

        // Condition numbers must be finite and >= 1 wherever a partition was accepted.
        for (int q = 0; q < NSURV; ++q) {
            if (co->okA[q] && !(std::isfinite(co->condA[q]) && co->condA[q] >= 1.0)) {
                std::cerr << "FAIL [" << ev.name << "] accepted photometric partition " << q
                          << " has bad condition number " << co->condA[q] << "\n";
                ++failures;
            }
            // A rejected partition must not carry a condition number at all. Without this, a
            // stale value from the previous event survives and looks like a real measurement.
            // Parameters outside a partition's active subset must never carry a number.
            if (co->okA[q]) {
                const auto act = activePhotParams(q, co->nepochA[SRUBIN], co->nepochA[SROMAN]);
                for (int k = 0; k < Nx; ++k) {
                    const bool isActive = std::find(act.begin(), act.end(), k) != act.end();
                    if (!isActive && co->Era[q][k] >= 0.0) {
                        std::cerr << "FAIL [" << ev.name << "] partition " << q << " reports "
                                  << "sigma(" << kPhotNames[k] << ")=" << co->Era[q][k]
                                  << " for a parameter it carries no information about\n";
                        ++failures;
                    }
                }
            }
            if (!co->okA[q] && co->condA[q] != -1.0) {
                std::cerr << "FAIL [" << ev.name << "] rejected photometric partition " << q
                          << " carries a stale condition number " << co->condA[q] << "\n";
                ++failures;
            }
        }

        // A partition with no epochs must be flagged not-characterizable, never given numbers.
        if (nR == 0 && co->okA[SROMAN] != 0) {
            std::cerr << "FAIL [" << ev.name << "] Roman has zero epochs but okA[SROMAN]="
                      << co->okA[SROMAN] << "; a partition with no data must not be usable\n";
            ++failures;
        }
        if (nL == 0 && co->okA[SRUBIN] != 0) {
            std::cerr << "FAIL [" << ev.name << "] Rubin has zero epochs but okA[SRUBIN]="
                      << co->okA[SRUBIN] << "\n";
            ++failures;
        }
        // Epoch counts must partition exactly.
        if (co->nepochA[SJOINT] != co->nepochA[SRUBIN] + co->nepochA[SROMAN]) {
            std::cerr << "FAIL [" << ev.name << "] epoch counts do not partition\n";
            ++failures;
        }
    }

    std::cout << "\n#\n";
    if (failures == 0) {
        std::cout << "# PASS -- all assertions held\n";
        return 0;
    }
    std::cout << "# FAIL -- " << failures << " assertion failure(s)\n";
    return 1;
}
