// Fisher-matrix regression fixture.
//
// Hand-builds a small set of microlensing events with known parameters, synthesises their light
// curves on representative Roman and Rubin cadences, runs the real FisherM and ErrorCal, and
// prints the sigmas in a stable, diffable format. Runs in seconds, bypassing the Monte Carlo,
// whose detection efficiency makes a Fisher-scored event cost tens of minutes.
//
// Reads NO data files (those are absent after a fresh clone): source magnitudes, blend fractions,
// cadences and photometric errors are hardcoded representative values.
//
// IS:     a regression harness -- run before and after a change to FisherM and diff the output.
// IS NOT: a validation of absolute precision. The events are hand-picked and the error model is a
//         flat stand-in for errlsstM / errRomanM.
//
// Usage: make fishertest && ./fishertest   (modes: --scale S, --sweep, --eigen [S],
//        --astro-variants, --sweep-astro, --finite-source, --lens-motion)
//
// The light-curve model's source is a point (pspl, the default here: the baseline every change to the
// Fisher matrices is diffed against) or a uniform disc (espl); --source-model pspl|espl picks it, in any
// mode. Every event of the table has the same source radius, kRho; --finite-source is the mode that puts the
// source on the lens, compares the two models and fits the source radius (ln rho) where the data see it. Under
// espl, --sweep runs the finite-source events (kFiniteEvents) and adds ln rho to the parameters swept.

#include <random>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_integration.h>
#include "common.h"
#include "types.h"
#include "events/lightcurve.h"
#include "events/vbm_model.h"
#include "surveys/visits.h"
#include "surveys/schedule.h"
#include "fisher/fisher.h"

namespace {

// ---------------------------------------------------------------------------------------------
// Sight line. The centre of the GBTDS fields. The fixture never runs the sightline setup, so main()
// points the light-curve model at it (vbmSetSightline): the Earth's orbit seen from here is the
// parallax of every event, and (kLon, kLat) fixes the axes (e_l, e_b) the event's xi is measured on.
// ---------------------------------------------------------------------------------------------
constexpr double kLon = 1.0;   // Galactic longitude [deg]
constexpr double kLat = -1.5;  // Galactic latitude  [deg]

// ---------------------------------------------------------------------------------------------
// Cadence model. Simplified stand-ins for the two surveys' observing patterns -- enough to
// reproduce the structural feature that matters (Roman dense-but-seasonal, Rubin sparse-but-
// continuous), without needing the real visit lists.
// ---------------------------------------------------------------------------------------------

// Roman GBTDS: 6 high-cadence seasons of ~72 days around a multi-year mid-mission gap. The real
// F146 cadence is ~12-15 min; sampled here at 4 h for speed, preserving the seasonal structure.
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
// NOTE: Roman has no real astrometric error model yet; this is a placeholder consistent with
// src/sim/characterize.cpp.
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

// The t0 values carry a +0.37 d offset so no epoch of either cadence grid lands exactly on t0:
// there the mus1/mus2 terms vanish identically and FisherM's astrometric CHECK aborts.
const FixtureEvent kEvents[] = {
    {"short_inseason",  5.0,   36.37, 0.30, 0.30},
    {"short_ingap",     5.0,  900.37, 0.30, 0.30},
    {"mid_inseason",   25.0, 1636.37, 0.30, 0.15},
    {"mid_ingap",      25.0, 2500.37, 0.30, 0.15},
    {"long_inseason", 100.0, 1636.37, 0.30, 0.10},
    {"long_ingap",    100.0, 2500.37, 0.30, 0.10},
    {"verylong",      900.0, 1800.37, 0.30, 0.05},
};

// The finite-source events of --sweep under a disc: a tE 25 d event with piE 0.15 peaking inside a Roman
// season, the source of radius rho passing the lens at u0 = z rho, with z and rho such that the data measure ln rho:
// a dwarf, a mid-sized star and a giant, near and outside the disc. The sweep drops the prior on ln rho, so that
// sigma(ln rho) is the data's, and keeps it wherever the matrix inverts (sigLnRhoMax infinite): the step, not the
// adoption, is what it checks. The radii are away from the
// rows of the library's disc table (log-spaced; 1e-3 and 1e-2 lie on one, where a central difference in ln rho
// straddles a row and averages the two cells beside it). Inside the disc the table's cells in u/rho (0.01) set
// the resolution of the u0 derivative, which the u0 step (absolute) cannot follow for a dwarf. kEvents' own
// source (kRho, 250 radii from the lens) never fits ln rho.
struct FiniteEvent { FixtureEvent ev; double rho; };
const FiniteEvent kFiniteEvents[] = {
    {{"fs_dwarf_z2",    25.0, 1636.37, 2.0 * 1.5e-3, 0.15}, 1.5e-3},
    {{"fs_mid_z0.8",    25.0, 1636.37, 0.8 * 1.5e-2, 0.15}, 1.5e-2},
    {{"fs_mid_z1.5",    25.0, 1636.37, 1.5 * 1.5e-2, 0.15}, 1.5e-2},
    {{"fs_mid_z2",      25.0, 1636.37, 2.0 * 1.5e-2, 0.15}, 1.5e-2},
    {{"fs_mid_z4",      25.0, 1636.37, 4.0 * 1.5e-2, 0.15}, 1.5e-2},
    {{"fs_giant_z1.5",  25.0, 1636.37, 1.5 * 5.0e-2, 0.15}, 5.0e-2},
    {{"fs_giant_z2",    25.0, 1636.37, 2.0 * 5.0e-2, 0.15}, 5.0e-2},
    {{"fs_giant_z3",    25.0, 1636.37, 3.0 * 5.0e-2, 0.15}, 5.0e-2},
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
constexpr double kRho       = 1.2e-3; // source radius in Einstein radii: a 1 Rsun dwarf at 8 kpc, theta* = 5.8e-4
                                      // mas, against thetaE = 0.5 mas. Read by the disc model only; the events
                                      // pass the source at u0 = 0.3, 250 radii away, where it changes nothing

// Fill in the geometry and stellar quantities that the model (evaluateTrack, evaluateModel) and
// ErrorCal() read.
void setupStatic(source& s, lens& l)
{
    s.lon = kLon;
    s.lat = kLat;

    s.Ds  = kDs;
    l.Dl  = kDl;
    s.xi  = kXi;
    l.tetE = kTetE;
    s.rho = kRho;

    s.mbs[0] = kMbsRubin;  s.fb[0] = kFbRubin;
    s.mbs[1] = kMbsRoman;  s.fb[1] = kFbRoman;

    s.mus1 = kMus1;  s.mus2 = kMus2;
    l.mul1 = kMul1;  l.mul2 = kMul2;
    s.mus  = std::sqrt(s.mus1 * s.mus1 + s.mus2 * s.mus2);
    l.mul  = std::sqrt(l.mul1 * l.mul1 + l.mul2 * l.mul2);
    l.murel = std::sqrt((s.mus1 - l.mul1) * (s.mus1 - l.mul1)
                      + (s.mus2 - l.mul2) * (s.mus2 - l.mul2));
}

// True (noiseless) model magnitude for telescope tt at magnification Astar.
double modelMag(const source& s, int tt, double Astar)
{
    return s.mbs[tt] - 2.5 * std::log10(Astar * s.fb[tt] + 1.0 - s.fb[tt]);
}

// Append one epoch to the light-curve buffers. The stored magnitude and astrometric position are
// the NOISELESS model values at the true parameters, which is what a Fisher forecast wants (and
// what the one-sided stencils difference against).
void appendEpoch(source& s, lens& l, astromet& as, double tim, int tele, int& ndw, int season = -1)
{
    // Roman observes from L2, so its observer displacement differs from Rubin's.
    TrackPoint p;
    evaluateTrack(s, l, as, tele, &tim, 1, true, &p);

    l.timn[ndw] = tim;
    l.magn[ndw] = modelMag(s, tele, p.A);
    l.errm[ndw] = (tele == 0) ? kErrMagRubin : kErrMagRoman;
    l.soux[ndw] = p.pos1c;
    l.souy[ndw] = p.pos2c;
    l.erra[ndw] = (tele == 0) ? kErrAstRubin : kErrAstRoman;
    l.tele[ndw] = tele;
    // Season and roll of a Roman exposure; seasons alternate spring / autumn from the first.
    l.rseas[ndw] = (tele == 1) ? season : -1;
    l.rroll[ndw] = (tele == 1) ? season % 2 : -1;
    ++ndw;
}

// Half-width of the fit window around t0, in days; epochs outside it are dropped. A real fit uses
// a window scaled to the event duration. The window also avoids FisherM's astrometric
// CHECK(!(s.pos1c == l.soux[i] && s.pos2c == l.souy[i])), which fires for epochs so far from peak
// that the deflection (decaying as 1/u^2) is below double precision and a piE perturbation returns
// bit-identical coordinates.
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
    for (int si = 0; si < int(sizeof(kRomanSeasonStarts) / sizeof(double)); ++si) {
        const double start = kRomanSeasonStarts[si];
        for (double t = start; t <= start + kRomanSeasonLength; t += kRomanCadence) {
            if (t < 0.0 || t > Tobs) continue;
            if (t < tMin || t > tMax) continue;
            appendEpoch(s, l, as, t, 1, ndw, si);
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

// Recompute F[SRUBIN] + F[SROMAN] and compare against F[SJOINT] element by element; exact up to
// floating-point summation order, since each epoch feeds the joint and exactly one survey matrix.
bool checkAdditivity(const covarian& co, const char* evName, int dim, bool photometric)
{
    bool ok = true;
    double worst = 0.0;
    for (int j = 0; j < dim; ++j) {
        for (int k = 0; k < dim; ++k) {
            // The photometric matrices carry the prior on ln rho once each (F[9][9]); the data parts are what add.
            const auto prior = [&](int q) {
                return (photometric and j == 9 and k == 9 and co.rhoFit[q] and std::isfinite(co.rhoPriorSigLn))
                    ? 1.0 / (co.rhoPriorSigLn * co.rhoPriorSigLn) : 0.0;
            };
            const double joint = photometric
                ? gsl_matrix_get(co.inputA[SJOINT].get(), j, k) - prior(SJOINT)
                : gsl_matrix_get(co.inputB[SJOINT].get(), j, k);
            const double parts = photometric
                ? gsl_matrix_get(co.inputA[SRUBIN].get(), j, k) - prior(SRUBIN)
                    + gsl_matrix_get(co.inputA[SROMAN].get(), j, k) - prior(SROMAN)
                : gsl_matrix_get(co.inputB[SRUBIN].get(), j, k)
                    + gsl_matrix_get(co.inputB[SROMAN].get(), j, k);
            // (The prior is added and taken off again, so an F[9][9] far below it is only known to its round-off.)
            const double scale = std::fabs(joint) + std::fabs(parts) + prior(SJOINT) + prior(SRUBIN) + prior(SROMAN) + 1e-300;
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
// survey constrains. Partition the joint parameters into A (the single survey's active set) and B
// (the other survey's flux parameters). The other survey contributes nothing to B, so the joint
// information on A is
//     S = F_this[A,A] + ( F_other[A,A] - F_other[A,B] F_other[B,B]^-1 F_other[B,A] ),
// and the bracket is a Schur complement of a PSD matrix, hence PSD. So S >= F_this[A,A] in the
// Loewner order and sigma_joint <= sigma_single. A violation is a bug, never physics (it once
// caught fb0 perturbing s.fb[tt] instead of s.fb[0]).
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
            // A negative sigma means the parameter is outside this partition's active subset.
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

// Adding t0 as a free parameter can only loosen the bound on the others (Cramer-Rao). The
// {u0, tE, fb, piE, xi} set is the leading submatrix of the information matrix, so both sigmas
// come from the same accumulated information on the same event. sigma(tE) with t0 free must be
// >= sigma(tE) with t0 fixed. Needs invert_matrix to operate on a scratch copy.
bool checkT0Marginalization(const covarian& co, const char* evName, bool verbose)
{
    if (Nx < 6) return true;  // nothing to compare against
    bool ok = true;

    for (int q = 0; q < NSURV; ++q) {
        if (!co.okA[q]) continue;

        // {u0, tE, fb, piE, xi} = indices 0-4, intersected with this partition's active subset.
        std::vector<int> sub;
        for (int k : activePhotParams(co, q)) if (k <= 4) sub.push_back(k);
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

// F * F^-1 must be the identity, validating the normalize / invert / rescale round trip.
// Needs invert_matrix to leave the accumulated information matrix intact. The tolerance scales with
// the condition number, since about log10(cond) digits are lost.
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
            ? activePhotParams(co, q) : std::vector<int>();
        const std::vector<int>& act = photometric ? actP : kAllAst;
        (void)dim;
        const int n = static_cast<int>(act.size());

        // Measure the residual in the NORMALIZED basis the inversion worked in, where every term
        // of the product is O(1): Ftilde = D F D, Ftilde^-1 = D^-1 F^-1 D^-1, s_i = 1/sqrt(F_ii).
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
        // Losing roughly log10(cond) digits from a ~1e-16 base is expected.
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

const char* kPhotNames[] = {"u0", "tE", "fb0", "piE", "xi", "t0", "mbs0", "fb1", "mbs1", "lnrho"};
const char* kAstNames[]  = {"tetE", "mus1", "mus2", "piE"};

// How many photometric parameters the tables show. ln rho (the last) exists only for a finite source, so a
// point-source run shows the nine it has: its tables are those of a model without the tenth.
int nPhotShown() { return (vbmSourceModel() == SourceModel::ESPL) ? Nx : Nx - 1; }

void printRow(const char* label, int nep, int ok, double cond,
              const std::vector<double>& era, const std::vector<double>& erb)
{
    std::cout << std::left << std::setw(10) << label
              << std::right << std::setw(7) << nep
              << std::setw(5) << ok
              << std::scientific << std::setprecision(2) << std::setw(11) << cond;
    std::cout << std::scientific << std::setprecision(3);
    for (int k = 0; k < nPhotShown(); ++k) std::cout << std::setw(12) << era[k];
    std::cout << std::setw(12) << erb[0] << "\n";
}

} // namespace

// ---------------------------------------------------------------------------------------------
// Finite-difference step-size convergence sweep.
//
// For each event, photometric parameter and survey partition, vary that parameter's step with
// everything else at default and record sigma. Truncation error grows with the step and round-off
// as 1/step, so a trustworthy step sits on the flat bottom of the U. mbs (indices 6 and 8) is the
// control: the model is exactly linear in it, so its curve must be perfectly flat.
int runSweep()
{
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();

    // Thirteen scales, geometric, spanning both walls of the U around the production step
    // (scale 1): round-off for scale << 1 (the difference of nearly equal model magnitudes is
    // amplified by 1/h) and truncation for scale >> 1 (the derivative picks up curvature and the
    // Fisher matrix gains spurious information). Scale 1e4 recovers the legacy steps.
    static const double kScales[] = {1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5,
                                     1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1,
                                     1.0,
                                     1.0e+1, 1.0e+2, 1.0e+3, 1.0e+4};

    std::cout << "event,tE,survey,param,param_idx,scale,sigma,cond,ok\n";

    // The table's events, each with the source of setupStatic; under a disc, the finite-source ones, which
    // fit ln rho (index 9, the last parameter, absent from a point-source sweep).
    std::vector<std::pair<FixtureEvent, double>> events;
    if (vbmSourceModel() == SourceModel::ESPL) for (const auto& f : kFiniteEvents) events.push_back({f.ev, f.rho});
    else                                       for (const auto& e : kEvents) events.push_back({e, kRho});

    co->sigLnRhoMax = std::numeric_limits<double>::infinity();
    co->rhoPriorSigLn = std::numeric_limits<double>::infinity();
    for (const auto& [ev, rho] : events) {
        for (int pidx = 0; pidx < nPhotShown(); ++pidx) {
            for (double sc : kScales) {
                // u0 is stepped by an absolute amount: a step reaching the lens would put the source on its
                // other side (u0 < 0), which the model does not take.
                if (pidx == 0 and FD_STEP_U0 * kFDStepScale * sc >= ev.u0) continue;
                setupStatic(*s, *l);
                s->rho = rho;
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
// Diagnostic: eigen-structure of the NORMALIZED photometric Fisher matrix, whose
// lambda_max/lambda_min is the condition number invert_matrix reports. Prints the spectrum and the
// eigenvector of lambda_min, i.e. WHICH parameter combination is degenerate.
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

        const auto act = activePhotParams(*co, SJOINT);
        const int dim = int(act.size());
        gsl_matrix* F = co->inputA[SJOINT].get();

        // D F D with D = diag(1/sqrt(F_ii)): unit diagonal, so off-diagonals are correlations.
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
// Sentinel discipline in ErrorCal.
//
// Era[]/Erb[] use -1.0 for "this partition could not measure this parameter"; every consumer must
// test for it before arithmetic. ErrorCal once took MIN(photometric, astrometric) piE, which
// prefers -1 over any real sigma when the astrometric matrix is singular. Each case sets the
// partition flags and sigmas by hand and checks what ErrorCal makes of them. The discriminating
// case: photometry good and astrometry singular must give a POSITIVE resu[3].
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

    std::cout << "\n# --- ErrorCal sentinel discipline ---\n";
    std::cout << "# " << std::left << std::setw(24) << "case"
              << std::setw(12) << "resu[3]" << std::setw(12) << "resu[9]"
              << std::setw(12) << "relMl[J]" << "\n";

    for (const auto& c : cases) {
        auto s  = std::make_unique<source>();
        auto l  = std::make_unique<lens>();
        auto co = std::make_unique<covarian>();
        setupStatic(*s, *l);
        l->u0 = 0.30; l->tE = 40.0; l->piE = 0.12;

        // ErrorCal RECOMPUTES Era/Erb as sqrt of the inverse-matrix diagonals, so the inverses
        // are what must be set up. Diagonal matrices make the tE-xi correlation term exactly zero.
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

        // A sigma is either a positive measurement or exactly the -1 sentinel.
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
        // The mass needs BOTH ingredients; missing either means -1.
        const bool massOK = (co->resu[9] > 0.0);
        if (massOK != c.wantMass) {
            std::cerr << "FAIL [sentinel/" << c.name << "] resu[9] = " << co->resu[9]
                      << ", expected " << (c.wantMass ? "a measurement" : "the -1 sentinel")
                      << "\n";
            ++fails;
        }
        // relMl is the same quantity as resu[9] computed independently; they must agree.
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
// Roman season clustering.
//
// buildRomanSchedule recovers the observing seasons from the epoch times, assuming within-season
// spacing and between-season gaps are cleanly separated by SEASON_GAP_MIN_DAYS; main() guards that
// assumption. Both halves are exercised on SYNTHETIC schedules (no data files): a healthy one that
// must cluster correctly, and a pathological one the guard must catch.
// ---------------------------------------------------------------------------------------------
namespace {
// Fill ro.tim with nSeasons seasons of length seasonLen sampled every cadence days, starting at
// day 730 and separated by gap days. The epochs repeat cyclically; buildRomanSchedule
// de-duplicates, as it does for the real file.
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
    std::cout << "\n# --- Roman season clustering ---\n";

    // Guard predicate, same form as in buildRomanSeasons (src/run/inputs.cpp).
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
            // dt_edge sign convention: negative in-season, positive outside; off-mission is
            // distinguishable from a mid-mission gap.
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

    {   // Pathological: an in-season cadence of 25 d exceeds SEASON_GAP_MIN_DAYS, so every epoch
        // looks like a season boundary and each "season" holds one epoch. maxInSeasonSpacing stays
        // 0 and minSeasonGap is 25 d, so only the zero-length seasons (minSeasonLength) trip the guard.
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

// ---------------------------------------------------------------------------------------------
// Astrometric noise variants W / N / P (AST_SIGC in config/parameters.h).
//
// Prints, per event and partition, sigma(tetE), the astrometric sigma(piE) and relMl under each
// variant, and ASSERTS sigma_W <= sigma_N <= sigma_P on tetE. N adds roll offsets and noise
// correlated within a day to W, and P adds season offsets and a larger sigma_c; neither extra
// parameters nor positively correlated noise can add information, so a violation is a bug. Also asserted:
// relMlV of every variant against the lens-mass error formed here from an independent inversion of the variant's
// matrix (route, thetaE-piE or ln rho-piE covariance; see lnCovTetPi), and relMlV[W] == relMl.
// ---------------------------------------------------------------------------------------------

// One stray epoch from a telescope must not make the joint photometric matrix singular. Take an
// event with many Roman epochs, keep exactly ONE Rubin epoch, and require that the joint matrix
// inverts and that sigma_joint <= sigma_Roman for tE.
bool checkFewEpochTelescope()
{
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();
    const auto& ev = kEvents[0];                       // short_inseason
    setupStatic(*s, *l);
    l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
    int nL = 0, nR = 0;
    int ndw = buildLightCurve(*s, *l, *as, nL, nR);
    // Keep the Roman epochs and the first Rubin epoch only, compacting the buffers in place.
    int w = 0; bool keptRubin = false;
    for (int i = 0; i < ndw; ++i) {
        if (l->tele[i] == 0) { if (keptRubin) continue; keptRubin = true; }
        l->timn[w] = l->timn[i]; l->magn[w] = l->magn[i]; l->errm[w] = l->errm[i];
        l->soux[w] = l->soux[i]; l->souy[w] = l->souy[i]; l->erra[w] = l->erra[i];
        l->tele[w] = l->tele[i]; l->rseas[w] = l->rseas[i]; l->rroll[w] = l->rroll[i];
        ++w;
    }
    ndw = w;
    l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
    s->xi = kXi; s->fb[0] = kFbRubin; s->fb[1] = kFbRoman; s->mbs[0] = kMbsRubin; s->mbs[1] = kMbsRoman;
    FisherM(*s, *l, *as, *co, ndw);
    ErrorCal(*co, *l, *s);
    const bool ok = co->okA[SJOINT] == 1 && co->okA[SROMAN] == 1
                 && co->Era[SJOINT][1] <= co->Era[SROMAN][1] * (1.0 + 1e-9);
    std::cout << "# --- one stray Rubin epoch ---\n# okA_J " << co->okA[SJOINT]
              << "  okA_R " << co->okA[SROMAN] << "  sigma_tE J/R " << co->Era[SJOINT][1] << " / "
              << co->Era[SROMAN][1] << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// Step-size sweep of the ASTROMETRIC derivatives, the counterpart of --sweep. One parameter's step
// is scaled at a time; prints CSV of sigma (variant W) per event, partition, parameter and scale.
int runSweepAstro()
{
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();
    static const double kScales[] = {1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 2.0};
    const char* pn[Ny] = {"tetE", "mus1", "mus2", "piE"};
    std::cout << "event,tE,survey,param,param_idx,scale,sigma,cond,ok\n";
    for (const auto& ev : kEvents) {
        for (int j = 0; j < Ny; ++j) {
            for (double sc : kScales) {
                setupStatic(*s, *l);
                l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
                int nL = 0, nR = 0;
                const int ndw = buildLightCurve(*s, *l, *as, nL, nR);
                l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
                s->xi = kXi; s->fb[0] = kFbRubin; s->fb[1] = kFbRoman;
                s->mbs[0] = kMbsRubin; s->mbs[1] = kMbsRoman;
                for (int k = 0; k < Ny; ++k) co->deltaScaleB[k] = 1.0;
                co->deltaScaleB[j] = sc;
                FisherM(*s, *l, *as, *co, ndw);
                ErrorCal(*co, *l, *s);
                for (int q = 0; q < NSURV; ++q) {
                    const char* qn = (q == SJOINT) ? "joint" : (q == SRUBIN ? "rubin" : "roman");
                    std::cout << ev.name << ',' << ev.tE << ',' << qn << ',' << pn[j] << ',' << j << ','
                              << sc << ',' << std::scientific << std::setprecision(10)
                              << co->Erb[q][j] << ',' << co->condB[q] << ',' << co->okB[q]
                              << std::defaultfloat << '\n';
                }
            }
        }
    }
    return 0;
}

int runAstroVariants()
{
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();
    const char* vn[NAVAR] = {"W", "N", "P"};
    int bad = 0, nMassChecked = 0, nMassAst = 0;
    std::cout << "# Step 3c astrometric noise variants: W white | N roll offsets + "
              << AST_SIGC_N << " mas/day | P roll offsets + " << AST_SIGC_P << " mas/day\n"
              << "# event,tE,partition,variant,ok,sig_tetE,sig_piE_ast,relMl,ratio_tetE_to_W\n";
    for (const auto& ev : kEvents) {
        setupStatic(*s, *l);
        l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
        int nL = 0, nR = 0;
        const int ndw = buildLightCurve(*s, *l, *as, nL, nR);
        l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
        s->xi = kXi; s->fb[0] = kFbRubin; s->fb[1] = kFbRoman;
        s->mbs[0] = kMbsRubin; s->mbs[1] = kMbsRoman;
        FisherM(*s, *l, *as, *co, ndw);
        ErrorCal(*co, *l, *s);
        for (int q = 0; q < NSURV; ++q) {
            const char* qn = (q == SJOINT) ? "joint" : (q == SRUBIN ? "rubin" : "roman");
            const double w0 = co->ErbV[AV_W][q][0];
            for (int v = 0; v < NAVAR; ++v)
                std::cout << ev.name << ',' << ev.tE << ',' << qn << ',' << vn[v] << ','
                          << co->okBV[v][q] << ',' << std::scientific << std::setprecision(4)
                          << co->ErbV[v][q][0] << ',' << co->ErbV[v][q][3] << ','
                          << co->relMlV[v][q] << ',' << std::fixed << std::setprecision(3)
                          << ((w0 > 0.0 && co->ErbV[v][q][0] > 0.0) ? co->ErbV[v][q][0] / w0 : -1.0)
                          << std::defaultfloat << '\n';
            // Independent check of the lens-mass error of every variant: invert FBV[v][q] here (GSL LU, not the
            // production invertNormalized), take rt (the astrometric thetaE of this variant combined with the
            // photometric one by inverse variance) and the better piE route, form c by the rule (photometric piE:
            // -wp cov(ln rho, ln piE); astrometric piE: (1 - wp) cov(thetaE, piE) / (thetaE piE)), and compare
            // sqrt(rt^2 + rpi^2 - 2c) with relMlV. Variant W is the main result: relMlV[AV_W] == relMl exactly.
            if (co->relMlV[AV_W][q] != co->relMl[q]) {
                std::cout << "#   FAIL " << ev.name << " " << qn << ": relMlV[W] " << co->relMlV[AV_W][q]
                          << " != relMl " << co->relMl[q] << "\n";
                ++bad;
            }
            for (int v = 0; v < NAVAR; ++v) {
                if (!co->okBV[v][q]) continue;
                gsl_matrix* M = gsl_matrix_alloc(Ny, Ny);
                gsl_matrix* Minv = gsl_matrix_alloc(Ny, Ny);
                gsl_permutation* perm = gsl_permutation_alloc(Ny);
                for (int j = 0; j < Ny; ++j)
                    for (int m = 0; m < Ny; ++m) gsl_matrix_set(M, j, m, co->FBV[v][q][j * Ny + m]);
                int sgn = 0;
                gsl_linalg_LU_decomp(M, perm, &sgn);
                gsl_linalg_LU_invert(M, perm, Minv);
                const double sT = std::sqrt(gsl_matrix_get(Minv, 0, 0)), sP = std::sqrt(gsl_matrix_get(Minv, 3, 3));
                const double cTP = gsl_matrix_get(Minv, 0, 3);
                gsl_matrix_free(M); gsl_matrix_free(Minv); gsl_permutation_free(perm);
                const double ra = sT / l->tetE, rph = co->relTetEph[q];
                const double rt = (rph >= 0.0) ? ra * rph / std::hypot(ra, rph) : ra;
                const double wp = (rph >= 0.0) ? rt * rt / (rph * rph) : 0.0;
                const double rp1 = (co->Era[q][3] >= 0.0) ? co->Era[q][3] / l->piE : -1.0, rp2 = sP / l->piE;
                const bool   phot = (rp1 >= 0.0 and rp1 < rp2);
                const double rpi = phot ? rp1 : rp2;
                const double c = phot ? -wp * gsl_matrix_get(co->inverA[q].get(), 9, 3) / l->piE * double(co->rhoMeas[q])
                                      : (1.0 - wp) * cTP / (l->tetE * l->piE);
                const double want = std::sqrt(rt * rt + rpi * rpi - 2.0 * c);
                if (!(std::fabs(co->relMlV[v][q] / want - 1.0) <= 1e-6)) {
                    std::cout << "#   FAIL " << ev.name << " " << qn << " " << vn[v] << ": relMlV " << co->relMlV[v][q]
                              << " vs independent " << want << (phot ? " (phot piE)" : " (ast piE)") << "\n";
                    ++bad;
                }
                ++nMassChecked;
                if (!phot) ++nMassAst;
            }
            for (int v = 1; v < NAVAR; ++v) {
                const double a0 = co->ErbV[v - 1][q][0], a1 = co->ErbV[v][q][0];
                if (co->okBV[v - 1][q] && co->okBV[v][q] && a1 < a0 * (1.0 - 1e-9)) {
                    std::cout << "#   FAIL " << ev.name << " " << qn << ": sigma_tetE " << vn[v]
                              << " < " << vn[v - 1] << "\n";
                    ++bad;
                }
            }
        }
    }
    std::cout << "# relMlV checked against an independent inversion: " << nMassChecked << " cases, " << nMassAst
              << " on the astrometric piE route\n";
    std::cout << (bad ? "# FAIL -- variant ordering violated\n" : "# PASS -- sigma_W <= sigma_N <= sigma_P everywhere\n");
    return bad ? 1 : 0;
}

// ---------------------------------------------------------------------------------------------
// The finite source. One event (tE 25 d, peaking inside a Roman season, piE 0.15) with a source of radius
// rho passing the lens at u0 = z rho: a dwarf (rho 1e-3, theta* = 5e-4 mas), a mid-sized star (1e-2) and a
// giant (5e-2, 43 Rsun at 8 kpc against thetaE 0.5 mas, a luminous bulge giant; 1e-2 is 8.6 Rsun), at z from
// the lens inside the disc out to where the disc has no effect. Each is forecast three ways, every one on its
// own noiseless data:
//   point  the point-source model,
//   fixed  the uniform disc with rho held at its true value, ln rho not a parameter (rhoFixed),
//   prod   the uniform disc as production forecasts it: ln rho a parameter of every partition, with the prior
//          RHO_PRIOR_SIGLN; rho is measured where the DATA-only sigma(ln rho) is at most SIGLNRHO_MAX,
//   free   the same without the prior (rhoPriorSigLn infinite), where its matrix inverts.
// Printed per partition: the peak magnification of the disc, dchiFS, the data-only and posterior sigma(ln rho) of
// production, whether rho is measured ("meas"), sigma(ln rho) of the prior-free fit, the fractional error of the
// astrometric (rel_a), photometric (rel_ph) and combined (rel_c) thetaE, the correlation of ln rho and piE, and
// the ratios sigma(fixed)/sigma(point), sigma(prod)/sigma(point) and sigma(free)/sigma(point) of u0, tE and piE.
// '-' = not measured.
// Asserted: (1) for u0 >= 30 rho every production sigma of u0, tE, piE is within 1% of the fixed-rho one; (2) where
// the prior-free matrix is well conditioned (sigma(ln rho) < 10) the data-only sigma of production equals it to
// 1e-6 (the Schur identity); (3) the posterior sigma(ln rho) is at most the prior's; (4) the joint data-only
// sigma(ln rho) is no worse than a single survey's; (5) F_joint = F_Rubin + F_Roman on the data part; fitting ln
// rho only loosens the other sigmas; a partition that does not measure rho reports no photometric thetaE and its
// combined thetaE is the astrometric one; the combined thetaE is no worse than either route; the point source's
// magnification is the disc model's own track through magnifOf; no partition falls back to rho fixed.
struct FsRun {
    double peak = 0.0;
    double sig[NSURV][3];                         // u0, tE, piE; -1 = not measured
    double tetE[NSURV], dchi[NSURV], slrData[NSURV], slrPost[NSURV], relPh[NSURV], relC[NSURV], corr[NSURV];
    int    fit[NSURV], meas[NSURV];
    bool   additive = true;
};

namespace {

enum class FsMode { Point, Fixed, Prod, Free };

FsRun fsForecast(FsMode mode, double rho, double z, source& s, lens& l, astromet& as, covarian& co)
{
    static const int kP[] = {0, 1, 3};     // u0, tE, piE
    vbmSetSourceModel(mode == FsMode::Point ? SourceModel::PSPL : SourceModel::ESPL);
    setupStatic(s, l);
    s.rho = rho;
    l.tE = 25.0;  l.t0 = 1636.37;  l.u0 = z * rho;  l.piE = 0.15;
    int nL = 0, nR = 0;
    const int ndw = buildLightCurve(s, l, as, nL, nR);
    l.tE = 25.0;  l.t0 = 1636.37;  l.u0 = z * rho;  l.piE = 0.15;
    s.xi = kXi;  s.fb[0] = kFbRubin;  s.fb[1] = kFbRoman;
    s.mbs[0] = kMbsRubin;  s.mbs[1] = kMbsRoman;
    TrackPoint p;
    evaluateTrack(s, l, as, 0, &l.t0, 1, true, &p);
    FsRun r;
    r.peak = p.A;
    const double kInf = std::numeric_limits<double>::infinity();
    co.rhoFixed      = (mode == FsMode::Fixed);
    co.rhoPriorSigLn = (mode == FsMode::Free) ? kInf : RHO_PRIOR_SIGLN;
    co.sigLnRhoMax   = (mode == FsMode::Free) ? kInf : SIGLNRHO_MAX;
    FisherM(s, l, as, co, ndw);
    ErrorCal(co, l, s);
    r.additive = checkAdditivity(co, "finite-source", Nx, true);
    for (int q = 0; q < NSURV; ++q) {
        for (int a = 0; a < 3; ++a) r.sig[q][a] = co.okA[q] ? co.Era[q][kP[a]] : -1.0;
        r.tetE[q] = co.okB[q] ? co.Erb[q][0] : -1.0;
        r.dchi[q] = co.dchiFS[q];
        r.fit[q]  = co.rhoFit[q];
        r.meas[q] = co.rhoMeas[q];
        r.slrData[q] = co.sigLnRhoData[q];
        r.slrPost[q] = (co.okA[q] and co.rhoFit[q]) ? co.Era[q][9] : -1.0;
        r.relPh[q] = co.relTetEph[q];
        r.relC[q]  = co.relTetEc[q];
        r.corr[q]  = co.corrRhoPiE[q];
    }
    co.rhoFixed = false;  co.rhoPriorSigLn = RHO_PRIOR_SIGLN;  co.sigLnRhoMax = SIGLNRHO_MAX;
    return r;
}

// The disc model's track is the point source's: at every epoch of the event, both observers, the model with
// the source radius set to zero (a point source) has the same u and the magnification magnifOf(u).
bool checkPointIsTrack(double rho, double z, source& s, lens& l, astromet& as)
{
    vbmSetSourceModel(SourceModel::ESPL);
    setupStatic(s, l);
    s.rho = rho;
    l.tE = 25.0;  l.t0 = 1636.37;  l.u0 = z * rho;  l.piE = 0.15;
    int nL = 0, nR = 0;
    const int ndw = buildLightCurve(s, l, as, nL, nR);
    l.tE = 25.0;  l.t0 = 1636.37;  l.u0 = z * rho;  l.piE = 0.15;
    s.xi = kXi;
    bool ok = true;
    for (int tt = 0; tt < 2; ++tt) {
        std::vector<double> t;
        for (int i = 0; i < ndw; ++i) if (l.tele[i] == tt) t.push_back(l.timn[i]);
        if (t.empty()) continue;
        std::vector<ModelPoint> disc(t.size()), pt(t.size());
        evaluateModel(s, l, as, tt, t.data(), int(t.size()), false, disc.data());
        s.rho = 0.0;
        evaluateModel(s, l, as, tt, t.data(), int(t.size()), false, pt.data());
        s.rho = rho;
        for (size_t i = 0; i < t.size(); ++i)
            if (pt[i].u != disc[i].u or pt[i].A != magnifOf(disc[i].u)) {
                std::cerr << "FAIL [finite/track] rho " << rho << " z " << z << " telescope " << tt << " epoch " << i
                          << ": u " << pt[i].u << " vs " << disc[i].u << ", A " << pt[i].A << " vs "
                          << magnifOf(disc[i].u) << "\n";
                ok = false;
                break;
            }
    }
    return ok;
}

// ---- The derivative of the disc magnification with respect to ln rho, three ways ----
// Exact uniform disc, independent of the library. The magnification of a disc of radius rho at distance u from a
// point lens is the mean of the point-lens magnification over the disc, A = (1/pi rho^2) int A_ps dS, and
// d A / d rho = (2 / rho)(B - A) with B the mean of A_ps over the RIM of the disc, so
//     d A / d ln rho = 2 (B - A).
// In polar coordinates about the lens, g(d) = A_ps(d) d = (d^2 + 2) / sqrt(d^2 + 4) integrates in closed form,
// int_0^x g = x sqrt(x^2 + 4) / 2, which gives the circles inside the disc entirely (d < rho - u); the others are
// partly inside, by the angle 2 atan2(rho sin psi, u + rho cos psi), psi the rim point's angle about the disc's
// centre, and d(psi)^2 = u^2 + rho^2 + 2 u rho cos psi. Both integrands are smooth in psi (the rim average
// has an integrable 1/d at u = rho). GSL adaptive Gauss-Kronrod, 1e-13.
double gDisc(double d) { return (d * d + 2.0) / std::sqrt(d * d + 4.0); }
double aPoint(double d) { return gDisc(d) / d; }

struct DiscPar { double u, rho; };

double rimIntegrand(double psi, void* p)
{
    const DiscPar& q = *static_cast<DiscPar*>(p);
    const double d2 = q.u * q.u + q.rho * q.rho + 2.0 * q.u * q.rho * std::cos(psi);
    return aPoint(std::sqrt(std::max(d2, 1e-300)));
}
double partialIntegrand(double psi, void* p)
{
    const DiscPar& q = *static_cast<DiscPar*>(p);
    const double d = std::sqrt(std::max(q.u * q.u + q.rho * q.rho + 2.0 * q.u * q.rho * std::cos(psi), 1e-300));
    const double ang = std::atan2(q.rho * std::sin(psi), q.u + q.rho * std::cos(psi));
    return gDisc(d) * 2.0 * ang * (q.u * q.rho * std::sin(psi) / d);
}

double integrate(double (*f)(double, void*), DiscPar& par)
{
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(2000);
    gsl_function F;  F.function = f;  F.params = &par;
    double res = 0.0, err = 0.0;
    gsl_set_error_handler_off();
    gsl_integration_qags(&F, 0.0, M_PI, 0.0, 1e-13, 2000, w, &res, &err);
    gsl_integration_workspace_free(w);
    return res;
}

double exactDisc(double u, double rho)
{
    DiscPar par{u, rho};
    const double x = std::max(rho - u, 0.0);
    return (2.0 * M_PI * 0.5 * x * std::sqrt(x * x + 4.0) + integrate(partialIntegrand, par)) / (M_PI * rho * rho);
}
double exactDerivative(double u, double rho)
{
    DiscPar par{u, rho};
    const double B = integrate(rimIntegrand, par) / M_PI;
    return 2.0 * (B - exactDisc(u, rho));
}

// The library's own finite difference at the production step, from one-epoch evaluations of its disc magnification.
double vbmDerivative(double u, double rho)
{
    ModelEvent ev{};
    ev.source = SourceModel::ESPL;
    const double step = FD_STEP_LNRHO * kFDStepScale;
    ev.rho = rho * std::exp(+step);  const double ap = vbmMagnification(ev, u);
    ev.rho = rho * std::exp(-step);  const double am = vbmMagnification(ev, u);
    return (ap - am) / (2.0 * step);
}

// u at which the library leaves its table for the quadrupole, quad * 10 = VBM_TOL (bisection).
double vbmSwitch(double rho)
{
    double lo = 1e-3 * rho, hi = 1e4;
    for (int i = 0; i < 200; ++i) {
        const double mid = std::sqrt(lo * hi);
        (discQuadrupole(mid, rho) * 10.0 < VBM_TOL ? hi : lo) = mid;
    }
    return hi;
}

int runRhoDerivative()
{
    vbmSetSourceModel(SourceModel::ESPL);
    int failures = 0;
    static const double kRho[] = {1e-5, 3e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 1.0};
    static const double kZ[]   = {0.5, 0.8, 1.5, 2, 3, 5, 7, 10, 15, 20, 30, 50, 100, 300, 1000};
    constexpr double kResolved = 1e-9;   // |exact dA/d ln rho| below this is lost in the exact integral's round-off
    std::cout << "\n# d A / d ln rho at u = (u/rho) rho: exact uniform disc (independent integration), the analytic quadrupole 2 quad,\n"
              << "# and the library's finite difference at the production step. err = relative to exact; '*' = the library is on its\n"
              << "# quadrupole branch there. Rows with |exact| < " << kResolved << " are unresolved by the exact integral and not judged.\n"
              << "# " << std::setw(8) << "rho" << std::setw(8) << "u/rho" << std::setw(13) << "exact" << std::setw(13) << "quad"
              << std::setw(10) << "err(quad)" << std::setw(13) << "vbm FD" << std::setw(10) << "err(FD)" << "\n";
    for (double rho : kRho) {
        double lastBadQuad = 0.0;   // largest u/rho at which the quadrupole is off by more than 1%
        const double us = vbmSwitch(rho);
        for (double z : kZ) {
            const double u = z * rho;
            const double ex = exactDerivative(u, rho), qd = 2.0 * discQuadrupole(u, rho), fd = vbmDerivative(u, rho);
            const bool resolved = std::fabs(ex) >= kResolved;
            const bool analytic = (u >= rhoDerivAnalyticU(rho));   // the branch production takes
            const double eq = resolved ? std::fabs(qd / ex - 1.0) : 0.0, ef = resolved ? std::fabs(fd / ex - 1.0) : 0.0;
            if (resolved and eq > 0.01) lastBadQuad = z;
            std::cout << "  " << std::scientific << std::setprecision(1) << std::setw(8) << rho << std::defaultfloat << std::setprecision(6) << std::setw(8) << z
                      << std::scientific << std::setprecision(4) << std::setw(13) << ex << std::setw(13) << qd
                      << std::setprecision(1) << std::setw(10) << eq << std::setprecision(4) << std::setw(13) << fd
                      << std::setprecision(1) << std::setw(10) << ef << (u > us ? " *" : "") << std::defaultfloat << "\n";
            // Sources larger than 0.1 Einstein radii (a few in 1e4 draws) are held to 6%: the expansion in rho and the
            // library's table are both coarser there.
            const double tolQ = (rho > 0.1) ? 0.06 : 0.01, tolF = (rho > 0.1) ? 0.06 : 0.02;
            if (resolved and analytic and eq > tolQ) {
                std::cerr << "FAIL [rho-deriv] rho " << rho << " u/rho " << z << ": quadrupole off by " << eq << "\n";
                ++failures;
            }
            if (resolved and not analytic and ef > tolF) {
                std::cerr << "FAIL [rho-deriv] rho " << rho << " u/rho " << z << ": library finite difference off by " << ef << "\n";
                ++failures;
            }
        }
        std::cout << "#   rho " << rho << ": quadrupole within 1% of exact beyond u/rho " << lastBadQuad
                  << "; the library's switch to its quadrupole is at u/rho " << us / rho
                  << "; production takes the analytic term from u/rho " << rhoDerivAnalyticU(rho) / rho << "\n";
        if (rhoDerivAnalyticU(rho) >= us) { std::cerr << "FAIL [rho-deriv] finite differences reach the library's switch at rho " << rho << "\n"; ++failures; }
    }

    // The spike: a finite difference that straddles the switch.
    const double rho = 6e-4, us = vbmSwitch(rho);
    double maxOld = 0.0, maxNew = 0.0, typical = 2.0 * discQuadrupole(us, rho);
    for (int k = -100; k <= 100; ++k) {
        const double u = us * (1.0 + k * 2e-9);
        const double oldD = vbmDerivative(u, rho);
        const double newD = (u >= rhoDerivAnalyticU(rho)) ? 2.0 * discQuadrupole(u, rho) : vbmDerivative(u, rho);
        maxOld = std::max(maxOld, std::fabs(oldD));
        maxNew = std::max(maxNew, std::fabs(newD));
    }
    std::cout << "# switch scan, rho " << rho << ", u = " << us << " (1 + k 2e-9), k = -100..100: true dA/dln rho " << typical
              << ", max |old finite difference| " << maxOld << " (x" << maxOld / typical << "), max |new| " << maxNew
              << " (x" << maxNew / typical << ")\n";
    if (!(maxOld > 100.0 * typical)) std::cout << "#   note: the old difference did not spike on this scan\n";
    if (!(maxNew < 1.5 * typical)) { std::cerr << "FAIL [rho-deriv] the derivative spikes across the library's switch\n"; ++failures; }
    return failures;
}

}  // namespace

int runFiniteSource()
{
    struct Grid { double rho; std::vector<double> z; };
    static const Grid kGrid[] = {
        {1.0e-3, {1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 10.0, 30.0}},
        {1.0e-2, {0.3, 0.8, 2.0, 3.0, 4.0, 4.5, 5.0, 6.0, 30.0}},
        {5.0e-2, {0.3, 0.8, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 30.0}},
    };
    const SourceModel keep = vbmSourceModel();
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();
    int failures = 0;
    const long fallbacks0 = rhoFallbackCount();

    std::cout << "# Finite source, tE 25 d, t0 1636.37 (in a Roman season), piE 0.15. ln rho is a parameter of every partition, with the\n"
              << "# prior sigma(ln rho) = " << RHO_PRIOR_SIGLN << "; rho is measured (meas) where the data-only sigma(ln rho) <= " << SIGLNRHO_MAX << ".\n"
              << "# sigma(u0), sigma(tE), sigma(piE): ratio to the point-source forecast, for rho fixed | production | prior-free fit.\n"
              << "# s_dat, s_post: data-only and posterior sigma(ln rho) of production; s_free: the prior-free fit's. '-' = not measured.\n"
              << "# rel_a, rel_ph, rel_c: fractional error of thetaE from astrometry, from photometry (theta*/rho,\n"
              << "# sigma(ln theta*) = " << THETASTAR_SIGLN << "), and combined (production).\n"
              << "# " << std::left << std::setw(7) << "rho" << std::setw(7) << "u0/rho" << std::setw(7) << "survey" << std::right
              << std::setw(9) << "Apk,disc" << std::setw(10) << "dchiFS" << std::setw(5) << "meas"
              << std::setw(10) << "s_dat" << std::setw(10) << "s_post" << std::setw(10) << "s_free" << std::setw(9) << "rel_a"
              << std::setw(9) << "rel_ph" << std::setw(9) << "rel_c"
              << std::setw(6) << "corr" << " |" << std::setw(28) << "u0 fix|prod|free" << std::setw(28) << "tE fix|prod|free"
              << std::setw(28) << "piE fix|prod|free" << "\n";

    auto ratio = [](double a, double b) {
        std::ostringstream o;
        if (a > 0.0 and b > 0.0) o << std::fixed << std::setprecision(3) << a / b; else o << "-";
        return o.str();
    };
    auto sci = [](double v, int prec = 2) {
        std::ostringstream o;
        if (v > 0.0) o << std::scientific << std::setprecision(prec) << v; else o << "-";
        return o.str();
    };

    for (const Grid& g : kGrid) {
        for (double z : g.z) {
            if (!checkPointIsTrack(g.rho, z, *s, *l, *as)) ++failures;
            const FsRun point = fsForecast(FsMode::Point, g.rho, z, *s, *l, *as, *co);
            const FsRun fixed = fsForecast(FsMode::Fixed, g.rho, z, *s, *l, *as, *co);
            const FsRun prod  = fsForecast(FsMode::Prod,  g.rho, z, *s, *l, *as, *co);
            const FsRun freeF = fsForecast(FsMode::Free,  g.rho, z, *s, *l, *as, *co);
            for (int q = 0; q < NSURV; ++q) {
                const char* qn = (q == SJOINT) ? "joint" : (q == SRUBIN ? "rubin" : "roman");
                const double relA = (prod.tetE[q] > 0.0) ? prod.tetE[q] / kTetE : -1.0;
                std::cout << "  " << std::left << std::setw(7) << g.rho << std::setw(7) << z << std::setw(7) << qn << std::right
                          << std::setprecision(4) << std::fixed << std::setw(9) << prod.peak
                          << std::setw(10) << std::setprecision(1) << prod.dchi[q] << std::setw(5) << prod.meas[q]
                          << std::setw(10) << sci(prod.slrData[q]) << std::setw(10) << sci(prod.slrPost[q])
                          << std::setw(10) << sci(freeF.slrPost[q]) << std::setw(9) << sci(relA) << std::setw(9) << sci(prod.relPh[q])
                          << std::setw(9) << sci(prod.relC[q]);
                if (prod.meas[q]) std::cout << std::setw(6) << std::setprecision(2) << std::fixed << prod.corr[q];
                else              std::cout << std::setw(6) << "-";
                std::cout << " |";
                for (int a = 0; a < 3; ++a)
                    std::cout << std::setw(9) << ratio(fixed.sig[q][a], point.sig[q][a]) << std::setw(9)
                              << ratio(prod.sig[q][a], point.sig[q][a]) << std::setw(10) << ratio(freeF.sig[q][a], point.sig[q][a]);
                std::cout << std::defaultfloat << std::setprecision(6) << "\n";

                // ---- assertions ----
                const std::string tag = "[finite/" + std::to_string(g.rho) + "," + std::to_string(z) + "," + qn + "] ";
                if (fixed.fit[q] or fixed.slrPost[q] != -1.0 or fixed.meas[q]) {
                    std::cerr << "FAIL " << tag << "rho fixed but ln rho is reported\n";
                    ++failures;
                }
                if (!prod.additive or !freeF.additive) ++failures;                                            // (5)
                if (prod.fit[q]) {
                    if (!(prod.slrPost[q] > 0.0 and prod.slrPost[q] <= RHO_PRIOR_SIGLN * (1.0 + 1e-12))) {    // (3)
                        std::cerr << "FAIL " << tag << "posterior sigma(ln rho) " << prod.slrPost[q] << " above the prior\n";
                        ++failures;
                    }
                    if (freeF.fit[q] and freeF.slrPost[q] > 0.0 and freeF.slrPost[q] < 10.0) {                 // (2)
                        if (!(std::fabs(prod.slrData[q] / freeF.slrPost[q] - 1.0) <= 1e-6)) {
                            std::cerr << "FAIL " << tag << "data-only sigma(ln rho) " << prod.slrData[q] << " vs prior-free " << freeF.slrPost[q] << "\n";
                            ++failures;
                        }
                    }
                    if (z >= 30.0)                                                                             // (1)
                        for (int a = 0; a < 3; ++a)
                            if (fixed.sig[q][a] > 0.0 and !(std::fabs(prod.sig[q][a] / fixed.sig[q][a] - 1.0) <= 0.01)) {
                                std::cerr << "FAIL " << tag << "u0 >= 30 rho: production sigma " << a << " " << prod.sig[q][a]
                                          << " vs rho fixed " << fixed.sig[q][a] << "\n";
                                ++failures;
                            }
                    // Fitting a parameter, even with a prior, cannot tighten the others.
                    for (int a = 0; a < 3; ++a)
                        if (fixed.sig[q][a] > 0.0 and prod.sig[q][a] < fixed.sig[q][a] * (1.0 - 1e-9)) {
                            std::cerr << "FAIL " << tag << "sigma(" << kPhotNames[a == 0 ? 0 : (a == 1 ? 1 : 3)]
                                      << ") fell when ln rho was fitted: " << fixed.sig[q][a] << " -> " << prod.sig[q][a] << "\n";
                            ++failures;
                        }
                }
                if (prod.meas[q]) {
                    if (!(prod.slrData[q] > 0.0 and prod.slrData[q] <= SIGLNRHO_MAX) or !(prod.relPh[q] >= THETASTAR_SIGLN) or
                        !(std::fabs(prod.corr[q]) <= 1.0)) {
                        std::cerr << "FAIL " << tag << "sigma_data " << prod.slrData[q] << ", rel_ph " << prod.relPh[q]
                                  << ", corr " << prod.corr[q] << "\n";
                        ++failures;
                    }
                    const double ra = relA, rp = prod.relPh[q], rc = prod.relC[q];
                    if (!(rc <= rp * (1.0 + 1e-12)) or (ra > 0.0 and !(rc <= ra * (1.0 + 1e-12)))) {
                        std::cerr << "FAIL " << tag << "combined thetaE error " << rc << " is worse than a route (" << ra << ", " << rp << ")\n";
                        ++failures;
                    }
                } else {
                    if (prod.relPh[q] != -1.0 or (relA > 0.0 and std::fabs(prod.relC[q] / relA - 1.0) > 1e-12)) {
                        std::cerr << "FAIL " << tag << "rho not measured but a photometric thetaE is used: rel_ph " << prod.relPh[q]
                                  << ", rel_c " << prod.relC[q] << ", rel_a " << relA << "\n";
                        ++failures;
                    }
                }
            }
            // (4) the joint data-only sigma(ln rho) is no worse than a single survey's.
            for (int q = 1; q < NSURV; ++q)
                if (prod.slrData[SJOINT] > 0.0 and prod.slrData[q] > 0.0 and prod.slrData[SJOINT] > prod.slrData[q] * (1.0 + 1e-9)) {
                    std::cerr << "FAIL [finite/" << g.rho << "," << z << "] data-only sigma(ln rho) joint " << prod.slrData[SJOINT]
                              << " > single survey " << prod.slrData[q] << "\n";
                    ++failures;
                }
        }
    }
    if (rhoFallbackCount() != fallbacks0) {
        std::cerr << "FAIL [finite] " << rhoFallbackCount() - fallbacks0 << " partition(s) fell back to rho fixed\n";
        ++failures;
    }
    failures += runRhoDerivative();
    vbmSetSourceModel(keep);
    std::cout << (failures ? "# FAIL -- finite-source assertions failed\n" : "# PASS -- all finite-source assertions held\n");
    return failures ? 1 : 0;
}

// ---------------------------------------------------------------------------------------------
// The lens's motion on the sky and its error (--lens-motion; ErrorCal, lensMotion in events/vbm_model.h). The
// fixture's lens is not self-consistent as it stands (its piE, tE and xi are set by hand, not drawn from a Dl and
// the proper motions), so here each event gets the lens its parameters imply: pi_rel = piE thetaE fixes Dl
// (1 / Dl = pi_rel + 1 / Ds), and mul is built from the geocentric motion tE and xi describe and the Earth's
// velocity at t0 -- by this test's own inverse of geocentricMotion, which the production function then has to
// undo: geocentricMotion(Dl, mus, mul) returns the event's tE and xi, and lensMotion returns mul and v_L. Neither
// touches the light curve, which depends on tE, xi, piE, thetaE and mus only.
// Per event and partition: the true |mu_L| [mas/day] and |v_L| [km/s], the errors relMuL = sigma(|mu_L|) / |mu_L|
// and sigVL = sqrt(lambda_max) of the 2x2 covariance of v_L [km/s] (the error along the worst-measured direction)
// by linear propagation, and a Monte Carlo of the same: kLmDraws parameter vectors
// from the covariance of the six parameters, pushed through the closed-form function; "mc/lin" is the sample's
// standard deviation of |mu_L| over the linear one, and sqrt of the largest eigenvalue of the SAMPLE covariance of
// the v_L draws over sigVL. The linear propagation is exact to first order, so
// expect ratios of 1 within a few percent where the relative errors are below ~0.1, and departures (the function
// is nonlinear: |.|, 1/pi_rel, the angle) above it. Then, the joint relMuL next to the previous formula's
// resu[13] (old) and the old formula with the sign of the tE-xi cross term of the mu_1 term corrected (old+),
// which show what the sign and the thetaE-mus correlations change.
// Asserted: the identity above (to 1e-9); relMuL and sigVL are the -1 sentinel exactly where a partition lacks a
// matrix or relMl; Sigma(thetaE) / thetaE^2 is relTetEc^2; on BOTH piE routes the variance of ln thetaE - ln piE is
// relMl^2 and, for the joint partition, that of ln thetaE + ln piE times (1 - D_L/D_S)^2 is rel_Dl^2 (resu[10]), to
// 1e-9 (relMass, lnCovTetPi: c = cov(ln thetaE, ln piE) of the combined thetaE with the route's piE). The same two
// are also sampled in the Monte Carlo (the columns "ln M", "ln pi_rel", mc/lin at the full and the scaled errors; the
// scaled ones asserted within 3%); the route and the correlation c / (rt rpi) are printed per row, and at least one
// row must take the astrometric route with |c / (rt rpi)| above 0.05, so the astrometric term is exercised.
// ---------------------------------------------------------------------------------------------
constexpr int kLmDraws = 200000;

int runLensMotion()
{
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();
    int failures = 0, nAstCorr = 0;
    std::mt19937_64 rng(20261012ULL);
    std::normal_distribution<double> gauss(0.0, 1.0);

    std::vector<std::pair<FixtureEvent, double>> events;   // the event, and its source radius
    for (const auto& ev : kEvents) events.push_back({ev, kRho});
    if (vbmSourceModel() == SourceModel::ESPL)
        for (const auto& fe : kFiniteEvents) events.push_back({fe.ev, fe.rho});

    std::cout << "# Lens motion: |mu_L| [mas/day], |v_L| [km/s] (heliocentric), relMuL = sigma(|mu_L|)/|mu_L|, sigVL = sqrt(lambda_max of cov(v_L)) [km/s], frame-independent.\n"
              << "# mc/lin: sample sd of " << kLmDraws << " draws from the parameter covariance over the linear propagation, for |mu_L| then, for v_L, sqrt(largest eigenvalue of the sample covariance) over sigVL\n"
              << "# ln M, ln pi_rel: the sample sd of ln thetaE -/+ ln piE (over the draws with thetaE, piE > 0) over the linear one, at the full and at the scaled errors; c/rt.rpi: the correlation of ln thetaE and ln piE used.\n"
              << "# (rel_v = sigVL/|v_L|, heliocentric); the last two at the covariance scaled in sigma by 'scale' to bring relMuL, rel_v (sigVL/|v_L|) and every parameter's error (xi in rad; mus against |mu_L|) to at most 2%.\n"
              << "# 'old', 'old+': the previous joint resu[13] and the same with the mu_1 tE-xi cross term's sign corrected. '-' = not measured.\n";
    for (const auto& [ev, rho] : events) {
        setupStatic(*s, *l);
        s->rho = rho;
        l->tE = ev.tE;  l->t0 = ev.t0;  l->u0 = ev.u0;  l->piE = ev.piE;
        int nL = 0, nR = 0;
        const int ndw = buildLightCurve(*s, *l, *as, nL, nR);
        l->tE = ev.tE;  l->t0 = ev.t0;  l->u0 = ev.u0;  l->piE = ev.piE;
        s->xi = kXi;  s->fb[0] = kFbRubin;  s->fb[1] = kFbRoman;

        // The self-consistent lens.
        const SkyFrame& F = vbmSkyFrame();
        l->pirel = l->piE * l->tetE;
        l->Dl = 1.0 / (l->pirel + 1.0 / s->Ds);
        double vE[2];
        earthVelocityNE(l->t0, vE);
        const double gl = -l->tetE / l->tE * std::cos(s->xi), gb = -l->tetE / l->tE * std::sin(s->xi);
        const double hN = F.toNE[0][0] * gl + F.toNE[0][1] * gb + vE[0] * l->pirel;
        const double hE = F.toNE[1][0] * gl + F.toNE[1][1] * gb + vE[1] * l->pirel;
        l->mul1 = s->mus1 + F.toNE[0][0] * hN + F.toNE[1][0] * hE;
        l->mul2 = s->mus2 + F.toNE[0][1] * hN + F.toNE[1][1] * hE;
        l->mul  = std::hypot(l->mul1, l->mul2);
        {
            const GeocentricMotion g = geocentricMotion(l->t0, l->tetE, l->pirel, s->mus1, s->mus2, l->mul1, l->mul2);
            if (std::fabs(g.tE / l->tE - 1.0) > 1e-9 or std::fabs(g.xi - s->xi) > 1e-9) {
                std::cerr << "FAIL [lens-motion/" << ev.name << "] the constructed lens does not give tE, xi back: tE "
                          << g.tE << " vs " << l->tE << ", xi " << g.xi << " vs " << s->xi << "\n";
                ++failures;
            }
        }

        FisherM(*s, *l, *as, *co, ndw);
        ErrorCal(*co, *l, *s);

        const MotionContext mc = motionContext(l->t0, 1.0 / s->Ds);
        const double p0[NLM] = {l->tE, s->xi, l->piE, l->tetE, s->mus1, s->mus2};
        const LensMotion m0 = lensMotion(mc, p0[0], p0[1], p0[2], p0[3], p0[4], p0[5]);
        const double mu0 = std::hypot(m0.mu[0], m0.mu[1]), v0 = std::hypot(m0.v[0], m0.v[1]);
        const double vconv = l->Dl * AU / (1000.0 * 24.0 * 3600.0);
        if (!(co->lmCheck < 1e-9) or std::fabs(v0 / (l->mul * vconv) - 1.0) > 1e-9 or std::fabs(m0.DL / l->Dl - 1.0) > 1e-9) {
            std::cerr << "FAIL [lens-motion/" << ev.name << "] closed form does not reproduce the draw: lmCheck "
                      << co->lmCheck << ", |v_L| " << v0 << " vs " << l->mul * vconv << ", D_L " << m0.DL << " vs " << l->Dl << "\n";
            ++failures;
        }

        std::cout << "\n# " << ev.name << "  tE=" << std::fixed << std::setprecision(1) << ev.tE << "  rho=" << std::scientific
                  << std::setprecision(2) << rho << "  |mu_L|=" << std::setprecision(5) << mu0 << " mas/day  |v_L|="
                  << std::fixed << std::setprecision(2) << v0 << " km/s  Dl=" << l->Dl << " kpc  lmCheck=" << std::scientific
                  << std::setprecision(1) << co->lmCheck << "\n";
        std::cout << std::left << std::setw(10) << "# survey" << std::right << std::setw(4) << "ok"
                  << std::setw(11) << "relMuL" << std::setw(11) << "mc/lin" << std::setw(11) << "sigVL" << std::setw(11) << "mc/lin"
                  << std::setw(10) << "rel_v" << std::setw(8) << "scale" << std::setw(9) << "mc/lin" << std::setw(9) << "mc/lin"
                  << std::setw(8) << "ln rho" << std::setw(8) << "piE" << std::setw(13) << "dlnMl-relMl"
                  << std::setw(9) << "c/rt.rpi" << std::setw(9) << "lnM" << std::setw(9) << "lnPirel" << std::setw(9) << "lnM" << std::setw(9) << "lnPirel" << "\n";
        static const char* kName[NSURV] = {"joint", "rubin", "roman"};
        for (int q = 0; q < NSURV; ++q) {
            double S[NLM][NLM];
            const bool have = lensMotionCovariance(*co, *l, q, S);
            const bool want = co->okA[q] and co->okB[q] and co->relMl[q] >= 0.0;
            if (have != want or (have != (co->relMuL[q] >= 0.0)) or (have != (co->sigVL[q] >= 0.0)) or
                (!have and (co->relMuL[q] != -1.0 or co->sigVL[q] != -1.0))) {
                std::cerr << "FAIL [lens-motion/" << ev.name << "] partition " << q << " sentinel: relMuL " << co->relMuL[q]
                          << ", sigVL " << co->sigVL[q] << ", relMl " << co->relMl[q] << "\n";
                ++failures;
            }
            if (!have) {
                std::cout << std::left << std::setw(10) << kName[q] << std::right << std::setw(4) << 0 << "         -          -          -          -\n";
                continue;
            }
            // Linear propagation, independently: relMuL / sigVL from ErrorCal against the same from S here.
            double Smu[2][2], Sv[2][2];
            lensMotionPropagate(mc, p0, S, Smu, Sv);
            const double ux = m0.mu[0] / mu0, uy = m0.mu[1] / mu0;
            const double linMu = std::sqrt(ux * ux * Smu[0][0] + 2 * ux * uy * Smu[0][1] + uy * uy * Smu[1][1]);
            const double linV  = std::sqrt(0.5 * (Sv[0][0] + Sv[1][1]) + std::sqrt(0.25 * (Sv[0][0] - Sv[1][1]) * (Sv[0][0] - Sv[1][1]) + Sv[0][1] * Sv[0][1]));

            // thetaE and ln Ml against relTetEc and relMl.
            const double tt = l->tetE, pp = l->piE;
            const double rp1 = (co->Era[q][3] >= 0.0) ? co->Era[q][3] / (std::fabs(l->piE) + eps) : -1.0;
            const double rp2 = (co->Erb[q][3] >= 0.0) ? co->Erb[q][3] / (std::fabs(l->piE) + eps) : -1.0;
            const bool photPi = (rp1 >= 0.0 and (rp2 < 0.0 or rp1 < rp2));
            const double vTet = S[LM_TETE][LM_TETE] / (tt * tt);
            const double vLnMl = vTet + S[LM_PIE][LM_PIE] / (pp * pp) - 2.0 * S[LM_TETE][LM_PIE] / (tt * pp);
            if (std::fabs(vTet / (co->relTetEc[q] * co->relTetEc[q]) - 1.0) > 1e-9) {
                std::cerr << "FAIL [lens-motion/" << ev.name << "] partition " << q << ": Sigma(thetaE)/thetaE^2 = " << vTet
                          << " vs relTetEc^2 = " << co->relTetEc[q] * co->relTetEc[q] << "\n";
                ++failures;
            }
            const double vLnPi = vTet + S[LM_PIE][LM_PIE] / (pp * pp) + 2.0 * S[LM_TETE][LM_PIE] / (tt * pp);
            const double cRho = S[LM_TETE][LM_PIE] / (tt * pp) / std::sqrt(vTet * S[LM_PIE][LM_PIE] / (pp * pp));   // c / (rt rpi)
            if (!photPi and std::fabs(cRho) > 0.05) ++nAstCorr;
            const double dMl = vLnMl / (co->relMl[q] * co->relMl[q]) - 1.0;
            if (q == SJOINT) {
                const double rd = (1.0 - l->Dl / s->Ds) * (1.0 - l->Dl / s->Ds) * vLnPi / (co->resu[10] * co->resu[10]) - 1.0;
                if (!(std::fabs(rd) <= 1e-9)) {
                    std::cerr << "FAIL [lens-motion/" << ev.name << "] joint: (1 - Dl/Ds)^2 var(ln thetaE + ln piE) = "
                              << (1.0 - l->Dl / s->Ds) * (1.0 - l->Dl / s->Ds) * vLnPi << " vs rel_Dl^2 = " << co->resu[10] * co->resu[10] << "\n";
                    ++failures;
                }
            }
            if (!(std::fabs(dMl) <= 1e-9)) {
                std::cerr << "FAIL [lens-motion/" << ev.name << "] partition " << q << ": var(ln thetaE - ln piE) = " << vLnMl
                          << " vs relMl^2 = " << co->relMl[q] * co->relMl[q] << "\n";
                ++failures;
            }
            if (std::fabs(linMu / mu0 / co->relMuL[q] - 1.0) > 1e-9 or std::fabs(linV / co->sigVL[q] - 1.0) > 1e-9) {
                std::cerr << "FAIL [lens-motion/" << ev.name << "] partition " << q << ": ErrorCal's relMuL, sigVL differ from the propagation here\n";
                ++failures;
            }

            // Monte Carlo from the covariance f^2 S: the correlation matrix's eigenvectors, scaled by the sigmas.
            // Returns the sample standard deviations of |mu_L| and |v_L| over the linear ones of the same covariance.
            double sg[NLM];
            for (int i = 0; i < NLM; ++i) sg[i] = std::sqrt(S[i][i]);
            gsl_matrix* C = gsl_matrix_alloc(NLM, NLM);
            for (int i = 0; i < NLM; ++i)
                for (int j = 0; j < NLM; ++j) gsl_matrix_set(C, i, j, S[i][j] / (sg[i] * sg[j]));
            gsl_vector* ev_ = gsl_vector_alloc(NLM);
            gsl_matrix* V = gsl_matrix_alloc(NLM, NLM);
            gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(NLM);
            gsl_eigen_symmv(C, ev_, V, w);
            double A[NLM][NLM];                                  // dp = f A z, A = sg * V sqrt(lambda)
            for (int i = 0; i < NLM; ++i)
                for (int k = 0; k < NLM; ++k) A[i][k] = sg[i] * gsl_matrix_get(V, i, k) * std::sqrt(std::max(gsl_vector_get(ev_, k), 0.0));
            gsl_eigen_symmv_free(w);  gsl_vector_free(ev_);  gsl_matrix_free(V);  gsl_matrix_free(C);
            auto monteCarlo = [&](double f, double& rMu, double& rV, double& rM, double& rP) {
                double n = 0, mM = 0, sM = 0, mV[2] = {0, 0}, cV[2][2] = {{0, 0}, {0, 0}};    // Welford
                double nl = 0, mLm = 0, sLm = 0, mLp = 0, sLp = 0;        // ln thetaE -/+ ln piE, over the draws where both are > 0
                for (int d = 0; d < kLmDraws; ++d) {
                    double z[NLM], p[NLM];
                    for (int k = 0; k < NLM; ++k) z[k] = gauss(rng);
                    for (int i = 0; i < NLM; ++i) {
                        p[i] = p0[i];
                        for (int k = 0; k < NLM; ++k) p[i] += f * A[i][k] * z[k];
                    }
                    const LensMotion m = lensMotion(mc, p[0], p[1], p[2], p[3], p[4], p[5]);
                    if (p[LM_TETE] > 0.0 and p[LM_PIE] > 0.0) {
                        const double lm = std::log(p[LM_TETE]) - std::log(p[LM_PIE]), lp = std::log(p[LM_TETE]) + std::log(p[LM_PIE]);
                        nl += 1.0;
                        const double dm = lm - mLm, dp = lp - mLp;
                        mLm += dm / nl;  sLm += dm * (lm - mLm);
                        mLp += dp / nl;  sLp += dp * (lp - mLp);
                    }
                    const double a = std::hypot(m.mu[0], m.mu[1]);
                    n += 1.0;
                    const double dM = a - mM;
                    mM += dM / n;  sM += dM * (a - mM);
                    const double d0 = m.v[0] - mV[0], d1 = m.v[1] - mV[1];
                    mV[0] += d0 / n;  mV[1] += d1 / n;
                    cV[0][0] += d0 * (m.v[0] - mV[0]);  cV[1][1] += d1 * (m.v[1] - mV[1]);  cV[0][1] += d0 * (m.v[1] - mV[1]);
                }
                rMu = std::sqrt(sM / (n - 1.0)) / (f * linMu);
                const double va = cV[0][0] / (n - 1.0), vd = cV[1][1] / (n - 1.0), vb = cV[0][1] / (n - 1.0);
                rV  = std::sqrt(0.5 * (va + vd) + std::hypot(0.5 * (va - vd), vb)) / (f * linV);
                rM  = std::sqrt(sLm / (nl - 1.0)) / (f * std::sqrt(vLnMl));
                rP  = std::sqrt(sLp / (nl - 1.0)) / (f * std::sqrt(vLnPi));
            };
            double rMu1, rV1, rMu2, rV2, rM1, rP1, rM2, rP2;
            monteCarlo(1.0, rMu1, rV1, rM1, rP1);
            // The same at a covariance scaled down until every error is at most 2%, if they are not already:
            // the regime where the linear propagation has to hold, whatever the fixture's errors are.
            // The condition is on the size of every parameter's error, not only on the two along the vectors: an error
            // across them (xi, or the proper motions against |mu_L|) is nonlinear in its own right.
            double worst = std::max(co->relMuL[q], co->sigVL[q] / v0);
            for (int i = 0; i < NLM; ++i) {
                const double ref = (i == LM_XI) ? 1.0 : (i == LM_MUS1 or i == LM_MUS2) ? mu0 : std::fabs(p0[i]);
                worst = std::max(worst, sg[i] / ref);
            }
            const double fSmall = std::min(1.0, 0.02 / worst);
            monteCarlo(fSmall, rMu2, rV2, rM2, rP2);
            if (std::fabs(rMu2 - 1.0) > 0.03 or std::fabs(rV2 - 1.0) > 0.03 or std::fabs(rM2 - 1.0) > 0.03 or std::fabs(rP2 - 1.0) > 0.03) {
                std::cerr << "FAIL [lens-motion/" << ev.name << "] partition " << q << ": at errors scaled by " << fSmall
                          << " the Monte Carlo/linear ratios are " << rMu2 << " (|mu_L|), " << rV2 << " (|v_L|), " << rM2
                          << " (ln M), " << rP2 << " (ln pi_rel)\n";
                ++failures;
            }

            std::cout << std::left << std::setw(10) << kName[q] << std::right << std::setw(4) << 1 << std::scientific << std::setprecision(3)
                      << std::setw(11) << co->relMuL[q] << std::fixed << std::setprecision(3) << std::setw(11) << rMu1
                      << std::scientific << std::setprecision(3) << std::setw(11) << co->sigVL[q]
                      << std::fixed << std::setprecision(3) << std::setw(11) << rV1
                      << std::setprecision(3) << std::setw(10) << co->sigVL[q] / v0
                      << std::setprecision(4) << std::setw(8) << fSmall << std::setprecision(3) << std::setw(9) << rMu2 << std::setw(9) << rV2
                      << std::setw(8) << (co->rhoMeas[q] ? "yes" : "no") << std::setw(8) << (photPi ? "phot" : "ast")
                      << std::scientific << std::setprecision(2) << std::setw(13) << dMl << std::fixed << std::setprecision(3)
                      << std::setw(9) << cRho << std::setw(9) << rM1 << std::setw(9) << rP1 << std::setw(9) << rM2 << std::setw(9) << rP2 << "\n";
        }

        // The previous joint formula, recomputed: its f1, f2 and sigmul from the joint matrices.
        {
            const auto& Era = co->Era[SJOINT];
            const auto& Erb = co->Erb[SJOINT];
            const double rtE = Era[1] / (std::fabs(l->tE) + eps), rTet = Erb[0] / (std::fabs(l->tetE) + eps);
            const double corr1 = gsl_matrix_get(co->inverA[SJOINT].get(), 1, 4) / std::fabs(Era[1] * Era[4]);
            const double tn = std::tan(s->xi), muGeo = l->tetE / l->tE;
            auto oldRel = [&](double sign1) {
                const double f1 = rTet * rTet + rtE * rtE + (Era[4] * tn) * (Era[4] * tn) + sign1 * 2.0 * rtE * Era[4] * tn * corr1;
                const double f2 = rTet * rTet + rtE * rtE + (Era[4] / tn) * (Era[4] / tn) - 2.0 * rtE * Era[4] / tn * corr1;
                const double s1 = std::sqrt(Erb[1] * Erb[1] + muGeo * muGeo * std::cos(s->xi) * std::cos(s->xi) * std::fabs(f1));
                const double s2 = std::sqrt(Erb[2] * Erb[2] + muGeo * muGeo * std::sin(s->xi) * std::sin(s->xi) * std::fabs(f2));
                return std::sqrt(s1 * s1 + s2 * s2) / (std::fabs(l->mul) + eps);
            };
            if (co->okA[SJOINT] and co->okB[SJOINT])
                std::cout << "#   joint relMuL: new " << std::scientific << std::setprecision(3) << co->relMuL[SJOINT]
                          << "   old " << oldRel(-1.0) << "   old+ " << oldRel(+1.0) << "   (resu[13] = " << co->resu[13] << ")\n";
            if (co->resu[13] != co->relMuL[SJOINT]) {
                std::cerr << "FAIL [lens-motion/" << ev.name << "] resu[13] is not the joint relMuL\n";
                ++failures;
            }
        }
    }
    std::cout << "\n# rows on the astrometric piE route with |c / (rt rpi)| > 0.05: " << nAstCorr << "\n";
    if (nAstCorr == 0) {
        std::cerr << "FAIL [lens-motion] no fixture row takes the astrometric route with a non-negligible c: the astrometric term is untested\n";
        ++failures;
    }
    std::cout << "\n#\n" << (failures ? "# FAIL -- lens-motion assertions failed\n" : "# PASS -- all lens-motion assertions held\n");
    return failures ? 1 : 0;
}

int main(int argc, char** argv)
{
    // --source-model pspl|espl in any mode, taken out of the argument list before the modes read it.
    SourceModel model = SourceModel::PSPL;
    {
        int kept = 1;
        for (int a = 1; a < argc; ++a) {
            if (std::string(argv[a]) != "--source-model") { argv[kept++] = argv[a]; continue; }
            const std::string v = (a + 1 < argc) ? argv[++a] : "";
            if      (v == "pspl") model = SourceModel::PSPL;
            else if (v == "espl") model = SourceModel::ESPL;
            else {
                std::cerr << "fixture error: --source-model needs pspl or espl, got '" << v << "'\n";
                return 2;
            }
        }
        argc = kept;
    }
    vbmSetSourceModel(model);
    vbmSetSightline(kLon, kLat);
    if (argc > 1 && std::string(argv[1]) == "--finite-source") return runFiniteSource();
    if (argc > 1 && std::string(argv[1]) == "--lens-motion") return runLensMotion();
    if (argc > 1 && std::string(argv[1]) == "--astro-variants") return runAstroVariants();
    if (argc > 1 && std::string(argv[1]) == "--sweep-astro") return runSweepAstro();
    if (argc > 1 && std::string(argv[1]) == "--sweep") return runSweep();
    if (argc > 1 && std::string(argv[1]) == "--eigen")
        return runEigen(argc > 2 ? std::atof(argv[2]) : 1.0);

    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();

    // --scale S: multiply EVERY finite-difference step by S before running the normal fixture, to
    // check that the joint-vs-single-survey comparison survives moving all the steps at once
    // (the sweep moves one at a time):
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

    // Data-free unit checks, run first so a failure is seen before the long event table.
    if (!checkSentinelDiscipline()) ++failures;
    if (!checkSeasonClustering())   ++failures;
    if (!checkFewEpochTelescope())  ++failures;

    if (model == SourceModel::ESPL)
        std::cout << "# SOURCE MODEL espl: a uniform disc of radius rho = " << kRho << " Einstein radii\n";
    std::cout << "# Fisher-matrix fixture -- synthetic events, no data files required\n"
              << "# Nx=" << nPhotShown() << " (photometric)  Ny=" << Ny << " (astrometric)\n"
              << "# One block per event, one row per survey partition.\n"
              << "# ok=1 usable inverse, ok=0 not characterizable (sigmas print as -1).\n"
              << "# Asserted: F[joint] == F[rubin] + F[roman] exactly, sigma_joint <= both, and\n"
              << "# sigma(tE) never decreases when t0 is marginalized over.\n";

    for (const auto& ev : kEvents) {
        setupStatic(*s, *l);
        l->tE  = ev.tE;
        l->t0  = ev.t0;
        l->u0  = ev.u0;
        l->piE = ev.piE;

        int nL = 0, nR = 0;
        const int ndw = buildLightCurve(*s, *l, *as, nL, nR);

        // Restore the true parameter point, since FisherM perturbs from whatever the structs hold.
        l->tE  = ev.tE;
        l->t0  = ev.t0;
        l->u0  = ev.u0;
        l->piE = ev.piE;
        s->xi  = kXi;
        s->fb[0] = kFbRubin;
        s->fb[1] = kFbRoman;

        // An epoch exactly at t0 makes the mus1/mus2 perturbations no-ops and trips FisherM's
        // astrometric CHECK; the t0 values are offset off-grid to avoid it.
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
        for (int k = 0; k < nPhotShown(); ++k) std::cout << std::setw(12) << kPhotNames[k];
        std::cout << std::setw(12) << kAstNames[0] << "\n";

        printRow("joint", co->nepochA[SJOINT], co->okA[SJOINT], co->condA[SJOINT], co->Era[SJOINT], co->Erb[SJOINT]);
        printRow("rubin", co->nepochA[SRUBIN], co->okA[SRUBIN], co->condA[SRUBIN], co->Era[SRUBIN], co->Erb[SRUBIN]);
        printRow("roman", co->nepochA[SROMAN], co->okA[SROMAN], co->condA[SROMAN], co->Era[SROMAN], co->Erb[SROMAN]);

        {
            static const char* kSyn[] = {"none", "both-alone", "rubin-only-alone",
                                         "roman-only-alone", "joint-only-RESCUE"};
            const int sc = synergyClass(*co);
            std::cout << "#   synergy: " << kSyn[sc];
            // What the joint fit bought over the best either survey managed alone.
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
            // A rejected partition must not carry a condition number (a stale value from the
            // previous event would look real); inactive parameters must not carry a sigma.
            if (co->okA[q]) {
                const auto act = activePhotParams(*co, q);
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

        // The finite source's results: a point source fits no ln rho and computes no dchiFS; a disc fits it in every
        // characterised partition, with a posterior sigma no larger than the prior's; the data-only sigma obeys the
        // Schur identity; rho is measured where that is at most sigLnRhoMax, and only there is there a photometric
        // thetaE or a correlation to read.
        for (int q = 0; q < NSURV; ++q) {
            const bool point = (model == SourceModel::PSPL);
            const bool fit = co->rhoFit[q] != 0, meas = co->rhoMeas[q] != 0;
            bool bad = false;
            if (point and (co->dchiFS[q] != -1.0 or fit or meas or co->sigLnRhoData[q] != -1.0)) bad = true;
            if (not point and co->okA[q] and not fit and not co->rhoFallback[q]) bad = true;
            if (fit) {
                const double sp = co->Era[q][9], sd = co->sigLnRhoData[q];
                if (not (co->okA[q] and sp > 0.0 and sp <= co->rhoPriorSigLn * (1.0 + 1e-12))) bad = true;
                if (sd > 0.0 and std::fabs(1.0 / (sd * sd) + 1.0 / (co->rhoPriorSigLn * co->rhoPriorSigLn) - 1.0 / (sp * sp)) > 1e-6 / (sp * sp)) bad = true;
                if (meas != (sd > 0.0 and sd <= co->sigLnRhoMax)) bad = true;
            } else if (co->Era[q][9] != -1.0 or co->sigLnRhoData[q] != -1.0 or meas) bad = true;
            if (not meas and (co->relTetEph[q] != -1.0 or not std::isnan(co->corrRhoPiE[q]))) bad = true;
            if (meas and not (co->relTetEph[q] >= THETASTAR_SIGLN)) bad = true;
            if (bad) {
                std::cerr << "FAIL [" << ev.name << "] partition " << q << " finite-source results inconsistent: dchiFS "
                          << co->dchiFS[q] << ", ln rho " << (fit ? "fitted" : "not fitted") << ", measured " << meas
                          << ", sigma(ln rho) posterior " << co->Era[q][9] << " data " << co->sigLnRhoData[q] << "\n";
                ++failures;
            }
        }

        // A partition with no epochs must be not-characterizable.
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
