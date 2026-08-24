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

// The joint forecast can never be worse than a single-survey forecast, for EVERY parameter that
// survey constrains -- including its own telescope-specific flux parameters.
//
// This is a theorem, not an expectation, and it is worth spelling out because an earlier version
// of this comment got it wrong. Partition the joint parameters into A (the single survey's active
// set) and B (the other survey's flux parameters). The other survey contributes nothing to B, so
// by the Schur complement the joint fit's effective information on A is
//
//     S = F_thisSurvey[A,A] + ( F_other[A,A] - F_other[A,B] F_other[B,B]^-1 F_other[B,A] )
//
// and the bracketed term is the Schur complement of the other survey's own information matrix,
// which is positive semi-definite. Hence S >= F_thisSurvey[A,A] in the Loewner order, so every
// diagonal element of the inverse can only shrink: sigma_joint <= sigma_single, always.
//
// A violation therefore means a bug, never a physical effect. It caught exactly that once: fb0
// was still perturbing s.fb[tt] rather than s.fb[0], so on Roman epochs parameters 2 and 7 both
// moved Roman's blend fraction and the joint matrix carried a duplicated direction.
bool checkJointNoWorse(const covarian& co, const char* evName,
                       const char* const* pnames, int dim, bool photometric)
{
    bool ok = true;
    for (int q = 1; q < NSURV; ++q) {
        const bool okJoint = photometric ? co.okA[SJOINT] : co.okB[SJOINT];
        const bool okPart  = photometric ? co.okA[q]      : co.okB[q];
        if (!okJoint || !okPart) continue;  // nothing to compare against
        for (int k = 0; k < dim; ++k) {
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

// ---------------------------------------------------------------------------------------------
// Step C3: finite-difference step-size convergence sweep.
//
// For each event, each photometric parameter and each survey partition, vary that parameter's
// step over about two decades with everything else at default, and record the recovered sigma.
//
// What to expect. Finite-difference error is U-shaped in the step h: truncation error grows with
// h (the model is not linear over the step), round-off error grows as 1/h (subtracting two nearly
// equal model values and dividing by a small number). The flat bottom is the plateau, and a
// trustworthy step sits in the middle of it.
//
// mbs (indices 6 and 8) is the control: the model magnitude is exactly linear in baseline
// magnitude, so its finite difference is exact at any step and its curve must be perfectly flat.
// Structure there means the sweep itself is broken, not the physics.
int runSweep()
{
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();

    // Twelve decades, geometric, centred so that BOTH failure modes are visible either side of
    // the production step (scale 1, after Step C3's kFDStepScale retune).
    //
    // The grid must span both walls of the U or the plot proves nothing -- a curve that is flat
    // everywhere you looked only means you did not look far enough. That was the original
    // mistake: the first sweep ran 0.1-10 around the LEGACY steps, saw no plateau anywhere, and
    // it took widening the window to discover the legacy steps were ~1e4 times too large and the
    // whole window sat on the truncation branch.
    //
    //   scale << 1  -- round-off. The step drives the model-magnitude difference down toward
    //                  double precision; subtracting two values agreeing to ~1e-9 leaves few
    //                  significant digits, amplified by the 1/h division. sigma collapses.
    //   scale >> 1  -- truncation. The secant stops matching the tangent; the derivative picks
    //                  up curvature and the Fisher matrix gains information that is not there.
    //                  scale 1e4 recovers the legacy steps this project used before Step C3.
    static const double kScales[] = {1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5,
                                     1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1,
                                     1.0,
                                     1.0e+1, 1.0e+2, 1.0e+3, 1.0e+4};

    std::cout << "event,tE,survey,param,param_idx,scale,sigma,cond,ok\n";

    for (const auto& ev : kEvents) {
        for (int pidx = 0; pidx < Nx; ++pidx) {
            for (double sc : kScales) {
                setupStatic(*s, *l);
                l->tE  = ev.tE;  l->t0 = ev.t0;
                l->u0  = ev.u0;  l->piE = ev.piE;

                int nL = 0, nR = 0;
                const int ndw = buildLightCurve(*s, *l, *as, nL, nR);

                l->tE  = ev.tE;  l->t0 = ev.t0;
                l->u0  = ev.u0;  l->piE = ev.piE;
                s->xi  = kXi;
                s->fb[0] = kFbRubin;  s->fb[1] = kFbRoman;
                s->mbs[0] = kMbsRubin; s->mbs[1] = kMbsRoman;

                for (int q = 0; q < Nx; ++q) co->deltaScale[q] = 1.0;
                co->deltaScale[pidx] = sc;

                FisherM(*s, *l, *as, *co, ndw);
                ErrorCal(*co, *l, *s);

                for (int q = 0; q < NSURV; ++q) {
                    const char* qn = (q == SJOINT) ? "joint" : (q == SRUBIN ? "rubin" : "roman");
                    std::cout << ev.name << ',' << ev.tE << ',' << qn << ','
                              << kPhotNames[pidx] << ',' << pidx << ','
                              << sc << ','
                              << std::scientific << std::setprecision(10) << co->Era[q][pidx]
                              << ',' << co->condA[q] << ',' << co->okA[q]
                              << std::defaultfloat << '\n';
                }
            }
        }
    }
    return 0;
}

// ---------------------------------------------------------------------------------------------
// Diagnostic: eigen-structure of the NORMALIZED photometric Fisher matrix.
//
// The condition number invert_matrix reports is lambda_max/lambda_min of exactly this matrix.
// This mode prints the whole spectrum plus the eigenvector belonging to lambda_min, i.e. the
// parameter combination the data constrains least. That vector is the physical content behind
// a large condition number: it names WHICH degeneracy the event suffers from.
int runEigen(double scale)
{
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();
    for (int q = 0; q < Nx; ++q) co->deltaScale[q] = scale;
    std::cout << "# all finite-difference steps scaled by " << scale << "\n";

    for (const auto& ev : kEvents) {
        setupStatic(*s, *l);
        l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
        int nL = 0, nR = 0;
        const int ndw = buildLightCurve(*s, *l, *as, nL, nR);
        l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
        s->xi = kXi; s->fb[0] = kFbRubin; s->fb[1] = kFbRoman;
        s->mbs[0] = kMbsRubin; s->mbs[1] = kMbsRoman;

        FisherM(*s, *l, *as, *co, ndw);

        const auto act = activePhotParams(SJOINT, co->nepochA[SRUBIN], co->nepochA[SROMAN]);
        const int dim = int(act.size());
        gsl_matrix* F = co->inputA[SJOINT].get();

        // D F D with D = diag(1/sqrt(F_ii)): unit diagonal, so every remaining off-diagonal
        // entry is a correlation coefficient and the units cancel out entirely.
        std::vector<double> sc(dim);
        bool bad = false;
        for (int i = 0; i < dim; ++i) {
            const double d = gsl_matrix_get(F, act[i], act[i]);
            if (!(d > 0.0)) { bad = true; break; }
            sc[i] = 1.0 / std::sqrt(d);
        }
        if (bad) { std::cout << "\n" << ev.name << ": a parameter has zero information\n"; continue; }

        gsl_matrix* Ft = gsl_matrix_alloc(dim, dim);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                gsl_matrix_set(Ft, i, j, gsl_matrix_get(F, act[i], act[j]) * sc[i] * sc[j]);

        gsl_vector* ew = gsl_vector_alloc(dim);
        gsl_matrix* evec = gsl_matrix_alloc(dim, dim);
        gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(dim);
        gsl_eigen_symmv(Ft, ew, evec, w);
        gsl_eigen_symmv_sort(ew, evec, GSL_EIGEN_SORT_VAL_DESC);
        gsl_eigen_symmv_free(w);

        std::cout << "\n# " << ev.name << "  tE=" << ev.tE << "  joint, dim=" << dim << "\n";
        std::cout << "  eigenvalues: ";
        for (int i = 0; i < dim; ++i)
            std::cout << std::scientific << std::setprecision(2) << gsl_vector_get(ew, i) << " ";
        std::cout << "\n  cond = lmax/lmin = " << std::scientific << std::setprecision(3)
                  << gsl_vector_get(ew, 0) / gsl_vector_get(ew, dim - 1) << "\n";
        std::cout << "  worst-constrained direction (eigenvector of lmin):\n    ";
        for (int i = 0; i < dim; ++i) {
            const double c = gsl_matrix_get(evec, i, dim - 1);
            if (std::fabs(c) > 0.15)
                std::cout << std::showpos << std::fixed << std::setprecision(2) << c
                          << std::noshowpos << "*" << kPhotNames[act[i]] << "  ";
        }
        std::cout << "\n";
        gsl_matrix_free(Ft); gsl_matrix_free(evec); gsl_vector_free(ew);
    }
    return 0;
}

// ---------------------------------------------------------------------------------------------
// Sentinel discipline in ErrorCal (Step D1).
//
// Era[]/Erb[] use -1.0 to mean "this partition could not measure this parameter". That is a
// SENTINEL, not a small error bar, and every consumer has to test for it before doing
// arithmetic. ErrorCal did not: it took the better of the two independent parallax routes with
// MIN(photometric, astrometric), which prefers -1 over any real sigma the moment the astrometric
// matrix is singular. The result was a negative sigma(piE) propagating into the mass and distance
// -- see DEVIATIONS entries 8 and 19.3.
//
// The live stub run that verified Step D1 did not contain a single event in that state, so this
// exercises it directly instead of waiting for one to turn up. Each case sets the partition flags
// and sigmas by hand and checks what ErrorCal makes of them.
//
// The discriminating assertion is case 2: with photometry good and astrometry singular, resu[3]
// must come out POSITIVE. Under the old MIN it was negative, every time.
// ---------------------------------------------------------------------------------------------
bool checkSentinelDiscipline()
{
    int fails = 0;

    // Case: {name, okA, okB, Era[3] (photometric piE), Erb[3] (astrometric piE), Erb[0] (tetE)}
    struct Case {
        const char* name;
        int    okA, okB;
        double eraPiE, erbPiE, erbTetE;
        bool   wantPiE;   //resu[3] should be a real measurement
        bool   wantMass;  //resu[9] / relMl should be a real measurement
    };
    const Case cases[] = {
        {"both routes available",      1, 1,  0.010, 0.020, 0.05, true,  true },
        {"astrometry singular",        1, 0,  0.010,  -1.0, -1.0, true,  false},
        {"photometry singular",        0, 1,   -1.0, 0.020, 0.05, true,  true },
        {"neither available",          0, 0,   -1.0,  -1.0, -1.0, false, false},
    };

    std::cout << "\n# --- ErrorCal sentinel discipline (Step D1) ---\n";
    std::cout << "# " << std::left << std::setw(24) << "case"
              << std::setw(12) << "resu[3]" << std::setw(12) << "resu[9]"
              << std::setw(12) << "relMl[J]" << "\n";

    for (const auto& c : cases) {
        auto s  = std::make_unique<source>();
        auto l  = std::make_unique<lens>();
        auto co = std::make_unique<covarian>();
        setupStatic(*s, *l);
        l->u0 = 0.30; l->tE = 40.0; l->piE = 0.12;

        // ErrorCal does NOT read Era/Erb -- it RECOMPUTES them from the inverse covariance
        // matrices, as sqrt of the diagonal. So the inverses are what has to be set up; writing
        // Era/Erb directly here would be silently overwritten. Diagonal matrices also make the
        // tE-xi correlation term ErrorCal forms from the (1,4) element exactly zero, which keeps
        // everything asserted below independent of it.
        for (int q = 0; q < NSURV; ++q) {
            co->okA[q] = c.okA;
            co->okB[q] = c.okB;
            co->nepochA[q] = 100;

            gsl_matrix_set_zero(co->inverA[q].get());
            for (int k = 0; k < Nx; ++k) gsl_matrix_set(co->inverA[q].get(), k, k, 0.01 * 0.01);
            if (c.eraPiE > 0.0)
                gsl_matrix_set(co->inverA[q].get(), 3, 3, c.eraPiE * c.eraPiE);

            gsl_matrix_set_zero(co->inverB[q].get());
            for (int k = 0; k < Ny; ++k) gsl_matrix_set(co->inverB[q].get(), k, k, 0.01 * 0.01);
            if (c.erbTetE > 0.0)
                gsl_matrix_set(co->inverB[q].get(), 0, 0, c.erbTetE * c.erbTetE);
            if (c.erbPiE > 0.0)
                gsl_matrix_set(co->inverB[q].get(), 3, 3, c.erbPiE * c.erbPiE);
        }

        ErrorCal(*co, *l, *s);

        std::cout << "# " << std::left << std::setw(24) << c.name
                  << std::setw(12) << co->resu[3]
                  << std::setw(12) << co->resu[9]
                  << std::setw(12) << co->relMl[SJOINT] << "\n";

        // A sigma is either a positive measurement or exactly the -1 sentinel. Never anything
        // in between, and never a negative number that is not the sentinel.
        if (co->resu[3] < 0.0 && co->resu[3] != -1.0) {
            std::cerr << "FAIL [sentinel/" << c.name << "] resu[3] = " << co->resu[3]
                      << " -- a negative sigma that is not the -1 sentinel. This is the"
                      << " unguarded MIN(photometric, astrometric) bug.\n";
            ++fails;
        }
        if (c.wantPiE && !(co->resu[3] > 0.0)) {
            std::cerr << "FAIL [sentinel/" << c.name << "] resu[3] = " << co->resu[3]
                      << " but at least one parallax route was available\n";
            ++fails;
        }
        if (!c.wantPiE && co->resu[3] != -1.0) {
            std::cerr << "FAIL [sentinel/" << c.name << "] resu[3] = " << co->resu[3]
                      << " but neither parallax route was available\n";
            ++fails;
        }
        // Case 1: MIN must pick the SMALLER of the two available routes, not merely one of them.
        if (c.okA && c.okB && std::fabs(co->resu[3] - c.eraPiE / std::fabs(l->piE)) > 1e-12) {
            std::cerr << "FAIL [sentinel/" << c.name << "] resu[3] = " << co->resu[3]
                      << " did not take the better of the two parallax routes\n";
            ++fails;
        }
        // The mass needs BOTH ingredients. Missing either means -1, never a number built
        // from a sentinel.
        const bool massOK = (co->resu[9] > 0.0);
        if (massOK != c.wantMass) {
            std::cerr << "FAIL [sentinel/" << c.name << "] resu[9] = " << co->resu[9]
                      << ", expected " << (c.wantMass ? "a measurement" : "the -1 sentinel")
                      << "\n";
            ++fails;
        }
        // relMl is the same physical quantity as resu[9] computed independently; they must agree.
        if (c.wantMass && std::fabs(co->relMl[SJOINT] - co->resu[9]) > 1e-12) {
            std::cerr << "FAIL [sentinel/" << c.name << "] relMl[SJOINT] = "
                      << co->relMl[SJOINT] << " disagrees with resu[9] = " << co->resu[9] << "\n";
            ++fails;
        }
        if (!c.wantMass && co->relMl[SJOINT] != -1.0) {
            std::cerr << "FAIL [sentinel/" << c.name << "] relMl[SJOINT] = "
                      << co->relMl[SJOINT] << ", expected the -1 sentinel\n";
            ++fails;
        }
    }
    return fails == 0;
}

// ---------------------------------------------------------------------------------------------
// Roman season clustering (Step D1).
//
// buildRomanSchedule recovers the observing seasons from the epoch times themselves rather than
// restating the generator's constants, so it depends on one assumption: that within-season epoch
// spacing and between-season gaps are cleanly separated by SEASON_GAP_MIN_DAYS. main() guards
// that assumption and refuses to run when it fails, because dt_edge and t0zone would still look
// entirely reasonable while meaning nothing.
//
// Both halves are exercised here on SYNTHETIC schedules, so this needs no data files: a healthy
// one that must cluster correctly, and a pathological one whose in-season cadence exceeds the
// threshold and which the guard must catch.
// ---------------------------------------------------------------------------------------------
namespace {
// Fill ro.tim with nSeasons seasons of length seasonLen sampled every cadence days, starting at
// day 730 and separated by gap days. The vector is longer than the synthetic schedule needs, so
// the epochs repeat cyclically -- buildRomanSchedule de-duplicates, exactly as it does for the
// real file where every field repeats every epoch.
void fillSynthetic(roman& ro, int nSeasons, double seasonLen, double cadence, double gap)
{
    std::vector<double> t;
    double start = 730.0;
    for (int i = 0; i < nSeasons; ++i) {
        for (double u = 0.0; u <= seasonLen + 1e-9; u += cadence) t.push_back(start + u);
        start += seasonLen + gap;
    }
    for (size_t i = 0; i < ro.tim.size(); ++i) ro.tim[i] = t[i % t.size()];
}
} //namespace

bool checkSeasonClustering()
{
    int fails = 0;
    std::cout << "\n# --- Roman season clustering (Step D1) ---\n";

    // Guard predicate, kept identical in form to the one in Bulge_LSST.cpp main().
    auto guardTrips = [](const RomanSchedule& sc) {
        return sc.seasons.size() < 2
            or sc.maxInSeasonSpacing >= SEASON_GAP_MIN_DAYS
            or sc.minSeasonGap       <= SEASON_GAP_MIN_DAYS
            or sc.minSeasonLength    <= 0.0;
    };

    {   // Healthy: 3 seasons of 70 d sampled every 5 d, separated by 110 d gaps.
        auto ro = std::make_unique<roman>();
        fillSynthetic(*ro, 3, 70.0, 5.0, 110.0);
        const RomanSchedule sc = buildRomanSchedule(*ro);
        std::cout << "# healthy   : " << sc.seasons.size() << " seasons, max in-season "
                  << sc.maxInSeasonSpacing << " d, min gap " << sc.minSeasonGap
                  << " d, shortest season " << sc.minSeasonLength
                  << " d, guard " << (guardTrips(sc) ? "TRIPS" : "passes") << "\n";

        if (sc.seasons.size() != 3) {
            std::cerr << "FAIL [seasons/healthy] recovered " << sc.seasons.size()
                      << " seasons, expected 3\n";
            ++fails;
        }
        if (guardTrips(sc)) {
            std::cerr << "FAIL [seasons/healthy] the guard rejected a well-separated schedule\n";
            ++fails;
        }
        if (!sc.seasons.empty()) {
            // Season 0 spans 730 -> 800; season 1 starts 70+110 = 180 d after season 0's start.
            if (std::fabs(sc.seasons[0].first - 730.0) > 1e-9
                or std::fabs(sc.seasons[0].second - 800.0) > 1e-9
                or std::fabs(sc.seasons[1].first - 910.0) > 1e-9) {
                std::cerr << "FAIL [seasons/healthy] wrong boundaries: season 0 ["
                          << sc.seasons[0].first << ", " << sc.seasons[0].second
                          << "], season 1 starts " << sc.seasons[1].first << "\n";
                ++fails;
            }
            // dt_edge sign convention: negative inside a season, positive outside, and
            // off-mission must be distinguishable from a mid-mission gap.
            struct P { double t0; int zone; bool negative; };
            const P probes[] = {
                {700.0,  T0_OFF_MISSION, false},  //before the first epoch
                {735.0,  T0_IN_SEASON,   true },  //5 d into season 0
                {855.0,  T0_IN_GAP,      false},  //mid-gap between seasons 0 and 1
                {9000.0, T0_OFF_MISSION, false},  //after the last epoch
            };
            for (const auto& p : probes) {
                const int    z  = sc.zone(p.t0);
                const double dt = sc.dtToSeasonEdge(p.t0);
                if (z != p.zone) {
                    std::cerr << "FAIL [seasons/healthy] t0=" << p.t0 << " zone=" << z
                              << ", expected " << p.zone << "\n";
                    ++fails;
                }
                if ((dt < 0.0) != p.negative) {
                    std::cerr << "FAIL [seasons/healthy] t0=" << p.t0 << " dt_edge=" << dt
                              << " has the wrong sign; negative means in-season\n";
                    ++fails;
                }
            }
            // Mid-gap is equidistant from both edges: 855 is 55 d from 800 and 55 d from 910.
            if (std::fabs(sc.dtToSeasonEdge(855.0) - 55.0) > 1e-9) {
                std::cerr << "FAIL [seasons/healthy] mid-gap dt_edge = "
                          << sc.dtToSeasonEdge(855.0) << ", expected 55\n";
                ++fails;
            }
        }
    }

    {   // Pathological: in-season cadence of 25 d exceeds SEASON_GAP_MIN_DAYS, so every epoch
        // looks like a season boundary and each "season" ends up holding exactly one epoch.
        //
        // This case is why minSeasonLength exists. The first two margins both look HEALTHY here:
        // maxInSeasonSpacing stays 0 (no spacing was ever classified as in-season, so it is never
        // updated) and minSeasonGap is 25 d, comfortably above the threshold. Only the zero-length
        // seasons give it away. Writing this test is what found that hole in the guard.
        auto ro = std::make_unique<roman>();
        fillSynthetic(*ro, 3, 70.0, 25.0, 110.0);
        const RomanSchedule sc = buildRomanSchedule(*ro);
        std::cout << "# degenerate: " << sc.seasons.size() << " seasons, max in-season "
                  << sc.maxInSeasonSpacing << " d, min gap " << sc.minSeasonGap
                  << " d, shortest season " << sc.minSeasonLength
                  << " d, guard " << (guardTrips(sc) ? "TRIPS" : "passes") << "\n";

        if (!guardTrips(sc)) {
            std::cerr << "FAIL [seasons/degenerate] a schedule whose in-season cadence exceeds "
                      << SEASON_GAP_MIN_DAYS << " d was accepted. dt_edge and t0zone would be "
                      << "fiction and nothing would say so.\n";
            ++fails;
        }
    }

    return fails == 0;
}

int main(int argc, char** argv)
{
    if (argc > 1 && std::string(argv[1]) == "--sweep") return runSweep();
    if (argc > 1 && std::string(argv[1]) == "--eigen")
        return runEigen(argc > 2 ? std::atof(argv[2]) : 1.0);

    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();

    // --scale S: multiply EVERY finite-difference step by S before running the normal fixture.
    //
    // The sweep varies one parameter's step at a time, which is right for locating each
    // parameter's plateau but cannot answer the question that actually matters: does the
    // headline joint-vs-single-survey comparison survive moving all the steps at once? sigma_p
    // comes from an inverse, so it depends on every row of the matrix, not just row p. This
    // mode moves them together so the full table -- including the joint/best-single sigma(tE)
    // ratio the thesis claim rests on -- can be diffed between step choices:
    //     ./fishertest > prod.txt && ./fishertest --scale 1e-3 > tuned.txt && diff prod.txt tuned.txt
    double allScale = 1.0;
    if (argc > 2 && std::string(argv[1]) == "--scale") {
        allScale = std::atof(argv[2]);
        if (!(allScale > 0.0)) {
            std::cerr << "fixture error: --scale needs a positive number, got '"
                      << argv[2] << "'\n";
            return 2;
        }
        for (int q = 0; q < Nx; ++q) co->deltaScale[q] = allScale;
        std::cout << "# ALL finite-difference steps scaled by " << allScale << "\n";
    }

    int failures = 0;

    // Data-free unit checks of the Step D1 logic. Run first so a failure here is seen before
    // the event table, which is long.
    if (!checkSentinelDiscipline()) ++failures;
    if (!checkSeasonClustering())   ++failures;

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
