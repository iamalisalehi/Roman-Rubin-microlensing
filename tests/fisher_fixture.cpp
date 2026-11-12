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
//        --astro-variants, --sweep-astro, --finite-source, --lens-motion, --blend)
//
// Rubin's bands: the default mode also checks the per-band fit (a blend fraction and baseline per Rubin band, indices
// 10-19; sigma(f_S) = FB_PRIOR_SIG on every active band's blend fraction, F146's too): prior off reproduces the
// single-band forecast bit for bit, a band at baseline or with fewer than kMinTeleEpochs visits changes no physical
// sigma, and a chromatic r/i/z event keeps the joint no worse than either survey. --sweep adds that event.
//
// The light-curve model's source is a point (pspl, the default here: the baseline every change to the
// Fisher matrices is diffed against) or a uniform disc (espl); --source-model pspl|espl picks it, in any
// mode. Every event of the table has the same source radius, kRho; --finite-source is the mode that puts the
// source on the lens, compares the two models and fits the source radius (ln rho) where the data see it. Under
// espl, --sweep runs the finite-source events (kFiniteEvents) and adds ln rho to the parameters swept.

#include <random>
#include <cstring>
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
    // Per-filter copies (index 2 = Rubin r, 6 = Roman F146), as func_source fills them
    s.blend[2] = kFbRubin;  s.magb[2] = kMbsRubin;
    s.blend[6] = kFbRoman;  s.magb[6] = kMbsRoman;
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
void appendEpoch(source& s, lens& l, astromet& as, double tim, int tele, int& ndw, int season = -1, int filt = -1)
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
    l.filt[ndw] = (filt >= 0) ? filt : ((tele == 0) ? 2 : 6); // default: Rubin r, Roman F146
    // Season and roll of a Roman exposure; seasons alternate spring / autumn from the first.
    l.rseas[ndw] = (tele == 1) ? season : -1;
    l.rroll[ndw] = (tele == 1) ? season % 2 : -1;
    ++ndw;
}

// The fixture's own 10-year window, fixed so that the regression table does not move when the simulated
// window (config's Tobs) is redefined.
constexpr double kFixtureTobs = 10.0 * year;

// Half-width of the fit window around t0, in days; epochs outside it are dropped. A real fit uses
// a window scaled to the event duration. The window also avoids FisherM's astrometric
// CHECK(!(s.pos1c == l.soux[i] && s.pos2c == l.souy[i])), which fires for epochs so far from peak
// that the deflection (decaying as 1/u^2) is below double precision and a piE perturbation returns
// bit-identical coordinates.
double fitHalfWindow(double tE)
{
    double w = 20.0 * tE;
    if (w < 100.0)  w = 100.0;   // floor: always keep some baseline either side
    if (w > kFixtureTobs)   w = kFixtureTobs;
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
            if (t < 0.0 || t > kFixtureTobs) continue;
            if (t < tMin || t > tMax) continue;
            appendEpoch(s, l, as, t, 1, ndw, si);
            ++nR;
        }
    }

    // Rubin: sparse sampling across the visible window of each year.
    for (int y = 0; y < kRubinYears; ++y) {
        const double yearStart = y * year;
        for (double t = yearStart; t <= yearStart + kRubinVisibleDays; t += kRubinCadence) {
            if (t < 0.0 || t > kFixtureTobs) continue;
            if (t < tMin || t > tMax) continue;
            appendEpoch(s, l, as, t, 0, ndw);
            ++nL;
        }
    }
    return ndw;
}

// Make the event chromatic: Rubin's visits cycle r, i, z (in stored order), each band with its own blend fraction and
// baseline magnitude -- f_S = 0.3, 0.5, 0.6, as the Rubin bands of a real blended source differ -- and the per-telescope
// copies of r's pair set to match, as func_source does.
void makeChromatic(source& s, lens& l, int ndw)
{
    s.blend[2] = 0.30;  s.magb[2] = 20.3;
    s.blend[3] = 0.50;  s.magb[3] = 19.9;
    s.blend[4] = 0.60;  s.magb[4] = 19.7;
    s.fb[0] = s.blend[2];  s.mbs[0] = s.magb[2];
    static const int kCycle[3] = {2, 3, 4};
    int c = 0;
    for (int i = 0; i < ndw; ++i) if (l.tele[i] == 0) l.filt[i] = kCycle[c++ % 3];
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
            // The photometric matrices carry the prior on ln rho once each (F[9][9]), and each matrix that holds a band the prior
            // on its blend fraction (F[fb_b][fb_b]): the joint and Rubin ones for Rubin's bands, the joint and Roman ones for
            // F146; the data parts are what add.
            const auto prior = [&](int q) {
                double p = 0.0;
                if (photometric and j == 9 and k == 9 and co.rhoFit[q] and std::isfinite(co.rhoPriorSigLn))
                    p += 1.0 / (co.rhoPriorSigLn * co.rhoPriorSigLn);
                if (photometric and j == k and q != SROMAN and std::isfinite(co.fbPriorSig))
                    for (int b = 0; b < kRubinBands; ++b)
                        if (co.nRubinBand[b] > 0 and j == fbIndexOfFilter(b)) p += 1.0 / (co.fbPriorSig * co.fbPriorSig);
                if (photometric and j == k and q != SRUBIN and std::isfinite(co.fbPriorSig) and j == fbIndexOfFilter(6)
                    and co.nepochA[SROMAN] >= kMinTeleEpochs)
                    p += 1.0 / (co.fbPriorSig * co.fbPriorSig);                  //Roman's F146, in the joint and Roman matrices
                return p;
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

const char* kPhotNames[] = {"u0", "tE", "fb0", "piE", "xi", "t0", "mbs0", "fb1", "mbs1", "lnrho",
                            "fb_u", "mbs_u", "fb_g", "mbs_g", "fb_i", "mbs_i", "fb_z", "mbs_z", "fb_y", "mbs_y"};
const char* kAstNames[]  = {"tetE", "mus1", "mus2", "piE"};

// How many photometric parameters the tables show: indices 0-9, the single-band ones. ln rho (9) exists only for a finite
// source, so a point-source run shows the nine it has: its tables are those of a model without the tenth. Rubin's other
// bands (10-19) get a line of their own, only when a band is in the fit (printRow).
int nPhotShown() { return (vbmSourceModel() == SourceModel::ESPL) ? 10 : 9; }

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
    // Rubin's other bands, when the partition fits any: their sigmas on a line of their own.
    bool bands = false;
    for (int k = 10; k < Nx; ++k) bands = bands or era[k] > 0.0;
    if (bands) {
        std::cout << "#   bands  ";
        for (int k = 10; k < Nx; ++k) if (era[k] > 0.0) std::cout << " " << kPhotNames[k] << "=" << era[k];
        std::cout << "\n";
    }
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
    // fit ln rho (index 9, the last single-band parameter, absent from a point-source sweep). The point-source sweep
    // adds a chromatic event (mid_inseason with Rubin's visits cycling r, i, z; makeChromatic), the only one with
    // Rubin bands beyond r: it sweeps the blend fraction and baseline of i and z (14-17) as well, mbs_i and mbs_z being
    // controls like mbs0 (exactly linear, so flat).
    struct SweepEvent { FixtureEvent ev; double rho; bool chromatic; };
    std::vector<SweepEvent> events;
    if (vbmSourceModel() == SourceModel::ESPL) for (const auto& f : kFiniteEvents) events.push_back({f.ev, f.rho, false});
    else {
        for (const auto& e : kEvents) events.push_back({e, kRho, false});
        events.push_back({{"mid_chromatic", kEvents[2].tE, kEvents[2].t0, kEvents[2].u0, kEvents[2].piE}, kRho, true});
    }

    co->sigLnRhoMax = std::numeric_limits<double>::infinity();
    co->rhoPriorSigLn = std::numeric_limits<double>::infinity();
    // The blend prior is dropped as well: with it a Rubin-poor partition's sigma(fb) is held by the prior and would pass
    // the plateau check whatever the step. A partition that does not invert without it is skipped (ok = 0 rows).
    co->fbPriorSig = std::numeric_limits<double>::infinity();
    for (const auto& [ev, rho, chromatic] : events) {
        std::vector<int> params;
        for (int pidx = 0; pidx < nPhotShown(); ++pidx) params.push_back(pidx);
        if (chromatic) for (const int pidx : {14, 15, 16, 17}) params.push_back(pidx);
        for (const int pidx : params) {
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
                if (chromatic) makeChromatic(*s, *l, ndw);

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
        co->nRubinBand[2] = 100;   // Rubin's r band is in the fit, so the joint and Rubin sets hold fb0, mbs0
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
        l->tele[w] = l->tele[i]; l->filt[w] = l->filt[i]; l->rseas[w] = l->rseas[i]; l->rroll[w] = l->rroll[i];
        ++w;
    }
    ndw = w;
    l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
    s->xi = kXi; s->fb[0] = kFbRubin; s->fb[1] = kFbRoman; s->mbs[0] = kMbsRubin; s->mbs[1] = kMbsRoman;
    // FisherM reads the per-filter blend and baseline (r is 2, F146 is 6): reset those with the per-telescope ones
    s->blend[2] = kFbRubin; s->blend[6] = kFbRoman; s->magb[2] = kMbsRubin; s->magb[6] = kMbsRoman;
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

// ---------------------------------------------------------------------------------------------
// The blend-fraction errors (--blend; the event table's sigfb0_J sigfb0_L sigfb1_J sigfb1_R columns are
// co.Era[q][2] and co.Era[q][7]). Per fixture event and partition: sigma(fb0) (Rubin's f_S) and sigma(fb1) (Roman's),
// asserted > 0 where the partition is characterised and the parameter is in activePhotParams, and the -1 sentinel
// otherwise; and the joint error no larger than either single survey's (sigma_J <= sigma_L for fb0, <= sigma_R for fb1,
// relative slack 1e-6). Then the accuracy at small f_S, where the new selection sends many more very blended Rubin
// events to the Fisher matrix: copies of one Rubin-covered event with Rubin's f_S set to 1e-2 ... 1e-5, the
// finite-difference steps scaled by 0.1, 1 and 10 (as --scale does); sigma(tE) and sigma(fb0) of the joint and of the
// Rubin matrix must agree across the scales within 1% (the Rubin matrix only for f_S >= 1e-4). The table is printed; the exit status is non-zero on a failure.
int runBlend()
{
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();
    int failures = 0;

    // One forecast at the event ev with Rubin's f_S = fRubin and every step scaled by stepScale.
    const auto forecast = [&](const FixtureEvent& ev, double fRubin, double stepScale) {
        for (int q = 0; q < Nx; ++q) co->deltaScale[q] = stepScale;
        setupStatic(*s, *l);
        s->fb[0] = fRubin;  s->blend[2] = fRubin;
        l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
        int nL = 0, nR = 0;
        const int ndw = buildLightCurve(*s, *l, *as, nL, nR);
        l->tE = ev.tE; l->t0 = ev.t0; l->u0 = ev.u0; l->piE = ev.piE;
        s->xi = kXi; s->fb[0] = fRubin; s->fb[1] = kFbRoman; s->blend[2] = fRubin;
        FisherM(*s, *l, *as, *co, ndw);
        ErrorCal(*co, *l, *s);
    };

    static const char* kPart[] = {"joint", "rubin", "roman"};
    std::cout << "# sigma(fb0), sigma(fb1) per fixture event and partition (-1 = not characterised)\n"
              << std::left << std::setw(16) << "# event" << std::setw(8) << "survey" << std::right << std::setw(4) << "ok"
              << std::setw(14) << "sigma(fb0)" << std::setw(14) << "sigma(fb1)" << "\n";
    std::cout << std::scientific << std::setprecision(5);
    for (const auto& ev : kEvents) {
        forecast(ev, kFbRubin, 1.0);
        for (int q = 0; q < NSURV; ++q) {
            std::cout << std::left << std::setw(16) << ev.name << std::setw(8) << kPart[q] << std::right << std::setw(4)
                      << int(co->okA[q]) << std::setw(14) << co->Era[q][2] << std::setw(14) << co->Era[q][7] << "\n";
            const auto act = activePhotParams(*co, q);
            for (const int k : {2, 7}) {
                const bool want = co->okA[q] && std::find(act.begin(), act.end(), k) != act.end();
                const double v = co->Era[q][k];
                if (want ? !(v > 0.0) : v != -1.0) {
                    std::cerr << "FAIL [" << ev.name << "] partition " << kPart[q] << ": sigma(" << kPhotNames[k] << ")=" << v
                              << (want ? " where it is characterised (want > 0)\n" : " where it is not (want the -1 sentinel)\n");
                    ++failures;
                }
            }
        }
        const double slack = 1e-6;
        if (co->okA[SJOINT] && co->okA[SRUBIN] && co->Era[SJOINT][2] > 0.0 && co->Era[SRUBIN][2] > 0.0
            && co->Era[SJOINT][2] > co->Era[SRUBIN][2] * (1.0 + slack)) {
            std::cerr << "FAIL [" << ev.name << "] sigma(fb0) joint " << co->Era[SJOINT][2] << " > Rubin " << co->Era[SRUBIN][2] << "\n";
            ++failures;
        }
        if (co->okA[SJOINT] && co->okA[SROMAN] && co->Era[SJOINT][7] > 0.0 && co->Era[SROMAN][7] > 0.0
            && co->Era[SJOINT][7] > co->Era[SROMAN][7] * (1.0 + slack)) {
            std::cerr << "FAIL [" << ev.name << "] sigma(fb1) joint " << co->Era[SJOINT][7] << " > Roman " << co->Era[SROMAN][7] << "\n";
            ++failures;
        }
    }

    // Small f_S in Rubin.
    const FixtureEvent& ev = kEvents[2];   // mid_inseason: Rubin and Roman both cover it
    static const double kFs[]     = {1e-2, 1e-3, 1e-4, 1e-5};
    static const double kScales[] = {0.1, 1.0, 10.0};
    std::cout << "\n# Rubin f_S small, event " << ev.name << ": sigma(tE) [d] and sigma(fb0) at step scales 0.1 / 1 / 10; spread = (max-min)/min\n"
              << std::left << std::setw(10) << "# f_S" << std::setw(8) << "survey";
    for (const char* p : {"tE", "fb0"}) {
        for (double sc : kScales) { std::ostringstream h; h << p << "@" << sc; std::cout << std::right << std::setw(13) << h.str(); }
        std::cout << std::setw(11) << (std::string(p) + " spread");
    }
    std::cout << "\n";
    for (const double f : kFs) {
        double v[2][2][3];   // [partition: joint, rubin][param: tE, fb0][scale]
        int okq[2][3];
        for (int is = 0; is < 3; ++is) {
            forecast(ev, f, kScales[is]);
            for (int ip = 0; ip < 2; ++ip) {
                const int q = ip == 0 ? SJOINT : SRUBIN;
                okq[ip][is] = co->okA[q];
                v[ip][0][is] = co->Era[q][1];
                v[ip][1][is] = co->Era[q][2];
            }
        }
        for (int ip = 0; ip < 2; ++ip) {
            const bool nearSingular = (ip == 1 && f < 5e-5);
            std::cout << std::left << std::scientific << std::setprecision(0) << std::setw(10) << f
                      << std::setw(8) << (ip == 0 ? "joint" : "rubin") << std::setprecision(5);
            for (int ik = 0; ik < 2; ++ik) {
                double lo = v[ip][ik][0], hi = v[ip][ik][0];
                for (int is = 0; is < 3; ++is) {
                    std::cout << std::right << std::setw(13) << v[ip][ik][is];
                    lo = std::min(lo, v[ip][ik][is]); hi = std::max(hi, v[ip][ik][is]);
                }
                const bool good = okq[ip][0] && okq[ip][1] && okq[ip][2] && lo > 0.0;
                const double spread = good ? (hi - lo) / lo : -1.0;
                std::cout << std::setw(11) << std::setprecision(3) << spread << std::setprecision(5);
                // Rubin alone at f_S = 1e-5 is near-singular (sigma(tE) ~ 1/f_S: no information on tE from Rubin alone), so its
                // step spread is not held to 1%; it is asserted to be far from characterised instead (sigma(tE) / tE > 100).
                if (nearSingular) {
                    if (ik == 0 && !(good && lo / ev.tE > 100.0)) {
                        std::cerr << "FAIL [f_S=" << f << "] rubin sigma(tE)/tE not > 100 at every step scale\n";
                        ++failures;
                    }
                } else if (!good || spread > 0.01) {
                    std::cerr << "FAIL [f_S=" << f << "] " << (ip == 0 ? "joint" : "rubin") << " sigma(" << (ik ? "fb0" : "tE")
                              << ") " << (good ? "differs across the step scales by " : "not characterised at every scale; spread ") << spread << "\n";
                    ++failures;
                }
            }
            std::cout << "\n";
        }
    }

    std::cout << "\n#\n" << (failures ? "# FAIL -- blend assertions failed\n" : "# PASS -- all blend assertions held\n");
    return failures ? 1 : 0;
}

// ---------------------------------------------------------------------------------------------
// Rubin's bands (a blend fraction and a baseline magnitude per band, indices 10-19; sigma(f_S) = FB_PRIOR_SIG on
// every active band's blend fraction, Roman's F146 included).
//
// The light-curve magnification is the same in every band; only m_b = mbs_b - 2.5 log10(fb_b A + 1 - fb_b) differs. The
// checks below hold the per-band machinery to what it must do: leave the single-band (r) forecast of an event as it was
// before the other bands existed when the new prior is switched off, ignore a band with too few visits, and not let a
// band that never sees the bump (an empty fb row, held invertible by the prior) move any physical parameter.
// ---------------------------------------------------------------------------------------------
namespace {

// The forecast (Era, okA, condA of indices 0-9) the f9e9938 build printed for every kEvents event and partition, at full
// precision (hexadecimal floats, exact), under the point source (espl 0) and the uniform disc (espl 1). Generated once
// from that build; never regenerated to make a check pass.
struct PriorOffRef { int espl; const char* name; int q; int ok; double cond; double era[10]; };
const PriorOffRef kPriorOffRef[] = {
    {0, "short_inseason", 0, 1, 0x1.78c60cd6937a6p+13, {0x1.371feaad9edcap-6, 0x1.7d900dccb5102p-5, 0x1.144befe4a49c3p-5, 0x1.a756ef50c1716p+1, 0x1.031400ff13a23p+2, 0x1.1f3e68a7100d8p-6, 0x1.9eb5c168e70fcp-9, 0x1.de5a580d7ff5cp-7, 0x1.542b44caa3e36p-12, -0x1p+0}},
    {0, "short_inseason", 1, 1, 0x1.961e7df161309p+15, {0x1.9321abdc73a5ep-3, 0x1.7469948d57f02p+3, 0x1.add3c0dd61f2bp-2, 0x1.8aa39ccc0b1eep+11, 0x1.aa019bd6dd984p+11, 0x1.2fe9b49d47322p-3, 0x1.3d3c5d58120a1p-8, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "short_inseason", 2, 1, 0x1.102f775e793eap+22, {0x1.7039736b864a4p-2, 0x1.789c33ceefa69p-3, -0x1p+0, 0x1.fd199915a1933p+5, 0x1.1348992fceb49p+6, 0x1.8078eddab7ce5p-5, -0x1p+0, 0x1.0ba8fd08093f5p-6, 0x1.a2e9f9e025429p-12, -0x1p+0}},
    {0, "short_ingap", 0, 1, 0x1.d127a2761a53p+17, {0x1.945274625c9fdp-2, 0x1.85439fa3df33bp+4, 0x1.5f15672165d1ep-1, 0x1.e484e8a3e6576p+11, 0x1.368bc753ec801p+11, 0x1.2117aaba3b0aep-3, 0x1.0706939a634efp-8, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "short_ingap", 1, 1, 0x1.d127a2761a53p+17, {0x1.945274625c9fdp-2, 0x1.85439fa3df33bp+4, 0x1.5f15672165d1ep-1, 0x1.e484e8a3e6576p+11, 0x1.368bc753ec801p+11, 0x1.2117aaba3b0aep-3, 0x1.0706939a634efp-8, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "short_ingap", 2, 0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "mid_inseason", 0, 1, 0x1.5501912a50c26p+14, {0x1.b143d4ba31e8cp-8, 0x1.9ee83d5757593p-1, 0x1.51fb4d7a36965p-7, 0x1.ce5b2930e473bp-1, 0x1.2847b9cefb787p+0, 0x1.1912d2f075d09p-7, 0x1.672c5f0b749d3p-10, 0x1.9d8d192e04359p-7, 0x1.070daf7482a0cp-12, -0x1p+0}},
    {0, "mid_inseason", 1, 1, 0x1.3c0fa19761426p+13, {0x1.2c394b62c495bp-4, 0x1.854429162452bp+3, 0x1.3ec87322ac4f4p-3, 0x1.842bf3deabdc3p+3, 0x1.ebc57c19fcb4cp+3, 0x1.7bb2c64b963aap-3, 0x1.9af0273a36429p-10, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "mid_inseason", 2, 1, 0x1.97d3353c64177p+16, {0x1.c3dc452ce149fp-7, 0x1.c2cd8481885d2p+0, -0x1p+0, 0x1.fae4f3e365e7p+0, 0x1.44638247d11cp+1, 0x1.752b44fc18868p-7, -0x1p+0, 0x1.7af7aa9c5afbep-6, 0x1.2456a828383edp-12, -0x1p+0}},
    {0, "mid_ingap", 0, 1, 0x1.6ac75d859088ep+20, {0x1.63d5335ffa8c2p+7, 0x1.0326d46f66daep+10, 0x1.1a8e75397d673p+6, 0x1.46bd007de24e9p+7, 0x1.fd5e38c80bc43p+6, 0x1.ac033f1818853p+5, 0x1.0bb360198aeeep-9, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "mid_ingap", 1, 1, 0x1.6ac75d859088ep+20, {0x1.63d5335ffa8c2p+7, 0x1.0326d46f66daep+10, 0x1.1a8e75397d673p+6, 0x1.46bd007de24e9p+7, 0x1.fd5e38c80bc43p+6, 0x1.ac033f1818853p+5, 0x1.0bb360198aeeep-9, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "mid_ingap", 2, 0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "long_inseason", 0, 1, 0x1.49f213ec1314p+13, {0x1.db11cd4fe906dp-9, 0x1.17dc148fd79abp+0, 0x1.2029b2efd00e2p-7, 0x1.07a0c48b5a9c2p-5, 0x1.0ea214716ca43p-4, 0x1.f22130a1e32a7p-6, 0x1.7a476b12efe8fp-11, 0x1.82367acfc91b1p-7, 0x1.f8601fb7819e1p-14, -0x1p+0}},
    {0, "long_inseason", 1, 1, 0x1.6e2085c41522bp+11, {0x1.0a822471eb76bp-5, 0x1.8d390cb432c42p+3, 0x1.19363b14b7514p-4, 0x1.8ec7227b7125ep-3, 0x1.07dfd2144bb07p-2, 0x1.dd09bcc126159p-2, 0x1.a4d6ea30efedep-11, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "long_inseason", 2, 1, 0x1.b85ebe6475ebfp+13, {0x1.123f52032ee26p-8, 0x1.1ba31ee3fb9a3p+0, -0x1p+0, 0x1.2088af8242508p-5, 0x1.3c7d33116b63ep-4, 0x1.09d2ed478a712p-5, -0x1p+0, 0x1.bdf6b516ca537p-7, 0x1.f8832a6c7ddf1p-14, -0x1p+0}},
    {0, "long_ingap", 0, 1, 0x1.7291971501f8ep+13, {0x1.a8b464277423p-3, 0x1.c1302b6acd443p+3, 0x1.303783ac90074p-2, 0x1.13ca5b4056536p-1, 0x1.02d5c5d5354ap-1, 0x1.f45582fbca219p+1, 0x1.cbafda6ffb86bp-11, 0x1.56f213202b1cdp-1, 0x1.01937dbc20fd2p-12, -0x1p+0}},
    {0, "long_ingap", 1, 1, 0x1.725217f7adf17p+13, {0x1.a9368f4909341p-3, 0x1.c178738687ef7p+3, 0x1.30d89af48730dp-2, 0x1.14a18e871031ep-1, 0x1.03a16a2a68d2cp-1, 0x1.f55435a6aefa8p+1, 0x1.cbafeed072abcp-11, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "long_ingap", 2, 1, 0x1.330b8e1844664p+26, {0x1.0f35291f8ea9ap+8, 0x1.f81007f129cbbp+13, -0x1p+0, 0x1.2624f0f1c1672p+5, 0x1.ee1ef5253bd94p+4, 0x1.058a2a4080f64p+11, -0x1p+0, 0x1.0b80e482bdf7ep+9, 0x1.ae08aa3bb8ee6p-9, -0x1p+0}},
    {0, "verylong", 0, 1, 0x1.2f9b4fe06f6f9p+16, {0x1.461d74f5a9174p-8, 0x1.77f7e00ddeb45p+3, 0x1.49a59ae41593ep-7, 0x1.64a2b9de432e6p-11, 0x1.982c10e8f6d25p-7, 0x1.36cc47466fb88p-2, 0x1.a5679e3ac8311p-10, 0x1.0a2d784901a3cp-6, 0x1.584c91234361ep-10, -0x1p+0}},
    {0, "verylong", 1, 1, 0x1.1eb9a1305d414p+13, {0x1.13dc02b751343p-6, 0x1.3c6c72ef652b8p+5, 0x1.18506f364414dp-5, 0x1.07c69c66fe236p-9, 0x1.788317abddaefp-6, 0x1.6477a4c7b649bp-1, 0x1.ef4690728e972p-9, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "verylong", 2, 1, 0x1.888c4ea07e9c3p+16, {0x1.70e758d9675ep-8, 0x1.9e668d3bdd608p+3, -0x1p+0, 0x1.02141d3daed65p-10, 0x1.54b030d438b07p-6, 0x1.84226893d880bp-2, -0x1p+0, 0x1.30f2c15fe6f76p-6, 0x1.774acf5e0ff35p-10, -0x1p+0}},
    {1, "short_inseason", 0, 1, 0x1.78c60d2b071abp+13, {0x1.371f8490a5168p-6, 0x1.7d8fd7057362p-5, 0x1.144bd0ca1bedap-5, 0x1.a75668c34e239p+1, 0x1.0313bbe02d9cep+2, 0x1.1f3e58f113768p-6, 0x1.9eb5c0d3efecdp-9, 0x1.de5a02bf42808p-7, 0x1.542b416e70b4dp-12, 0x1.ffffffffc2e08p-1}},
    {1, "short_inseason", 1, 1, 0x1.9627c3b7d7481p+15, {0x1.932502b59d135p-3, 0x1.746dcf4057e4p+3, 0x1.add7b60be5d51p-2, 0x1.8aa7e7c505f4p+11, 0x1.aa061ce2a40cbp+11, 0x1.2fe9fcff13246p-3, 0x1.3d3c6e117f4fdp-8, -0x1p+0, -0x1p+0, 0x1.fffffffffffcep-1}},
    {1, "short_inseason", 2, 1, 0x1.102f6ac92898fp+22, {0x1.7038f32f47a5p-2, 0x1.789bc8e15d8b7p-3, -0x1p+0, 0x1.fd18e995b1322p+5, 0x1.1348376b206e8p+6, 0x1.80785596ce4d4p-5, -0x1p+0, 0x1.0ba8df29f0b27p-6, 0x1.a2e9c38701f4bp-12, 0x1.ffffffffc9fa2p-1}},
    {1, "short_ingap", 0, 1, 0x1.d12b50578ce1ap+17, {0x1.94535c1f68f2fp-2, 0x1.854504624d4a4p+4, 0x1.5f165c9eb0b8dp-1, 0x1.e4869a8750c37p+11, 0x1.368ce9a286a4bp+11, 0x1.21194ff61bc9ap-3, 0x1.0706dc78f0138p-8, -0x1p+0, -0x1p+0, 0x1.fffffffffffd3p-1}},
    {1, "short_ingap", 1, 1, 0x1.d12b50578ce1ap+17, {0x1.94535c1f68f2fp-2, 0x1.854504624d4a4p+4, 0x1.5f165c9eb0b8dp-1, 0x1.e4869a8750c37p+11, 0x1.368ce9a286a4bp+11, 0x1.21194ff61bc9ap-3, 0x1.0706dc78f0138p-8, -0x1p+0, -0x1p+0, 0x1.fffffffffffd3p-1}},
    {1, "short_ingap", 2, 0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}},
    {1, "mid_inseason", 0, 1, 0x1.55028ff60bb7bp+14, {0x1.b1438e7098aeep-8, 0x1.9ee7bf4ceee48p-1, 0x1.51fb4d7750ad1p-7, 0x1.ce5aa3cad0acp-1, 0x1.28476460eda33p+0, 0x1.19129ddf5f6f8p-7, 0x1.672c5f12aebd3p-10, 0x1.9d8c8fe44bf1cp-7, 0x1.070daad5373ap-12, 0x1.ffffffffc8d08p-1}},
    {1, "mid_inseason", 1, 1, 0x1.3c0f8d5f7923cp+13, {0x1.2c38cf4c87636p-4, 0x1.8544075edd0c7p+3, 0x1.3ec835b724796p-3, 0x1.842bdde5ba6c5p+3, 0x1.ebc5614edf5cbp+3, 0x1.7bb285ef6c323p-3, 0x1.9af027ee7bfa7p-10, -0x1p+0, -0x1p+0, 0x1.ffffffffffe47p-1}},
    {1, "mid_inseason", 2, 1, 0x1.97d4f25cd525ap+16, {0x1.c3dc536372597p-7, 0x1.c2cd4764b3573p+0, -0x1p+0, 0x1.fae4b8b9ceed2p+0, 0x1.44635c3de0465p+1, 0x1.752b11747f9f1p-7, -0x1p+0, 0x1.7af750bcab465p-6, 0x1.2456a270a92c4p-12, 0x1.fffffffffc511p-1}},
    {1, "mid_ingap", 0, 1, 0x1.6ac7ba67176f8p+20, {0x1.63d5511126e16p+7, 0x1.0326f0bc9ac19p+10, 0x1.1a8e94311c1aap+6, 0x1.46bd24cadefe3p+7, 0x1.fd5e7d45ea84cp+6, 0x1.ac03486ab0a87p+5, 0x1.0bb3636f1f796p-9, -0x1p+0, -0x1p+0, 0x1p+0}},
    {1, "mid_ingap", 1, 1, 0x1.6ac7ba67176f8p+20, {0x1.63d5511126e16p+7, 0x1.0326f0bc9ac19p+10, 0x1.1a8e94311c1aap+6, 0x1.46bd24cadefe3p+7, 0x1.fd5e7d45ea84cp+6, 0x1.ac03486ab0a87p+5, 0x1.0bb3636f1f796p-9, -0x1p+0, -0x1p+0, 0x1p+0}},
    {1, "mid_ingap", 2, 0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}},
    {1, "long_inseason", 0, 1, 0x1.49f4f7a420ec6p+13, {0x1.db10fde9ed1a5p-9, 0x1.17dc09e9fd839p+0, 0x1.202979c4e4a7cp-7, 0x1.07a0b5d02a44dp-5, 0x1.0ea1f739d0563p-4, 0x1.f220cb630d19bp-6, 0x1.7a476a8810ad9p-11, 0x1.82361de37d1ap-7, 0x1.f8601fb01c9f5p-14, 0x1.fffffffffeb36p-1}},
    {1, "long_inseason", 1, 1, 0x1.6e20375241198p+11, {0x1.0a81b37b20658p-5, 0x1.8d38f94fe9f65p+3, 0x1.1936022b825dp-4, 0x1.8ec70eb31e507p-3, 0x1.07dfc87967127p-2, 0x1.dd095f5661ce5p-2, 0x1.a4d6e96a7a531p-11, -0x1p+0, -0x1p+0, 0x1.ffffffffffefbp-1}},
    {1, "long_inseason", 2, 1, 0x1.b86274aeb8292p+13, {0x1.123ed18fda548p-8, 0x1.1ba31565031fdp+0, -0x1p+0, 0x1.208899b92dac2p-5, 0x1.3c7d05e5b48f2p-4, 0x1.09d2b18ebabcp-5, -0x1p+0, 0x1.bdf63bd2e47ffp-7, 0x1.f8832a6837db3p-14, 0x1.ffffffffff8a4p-1}},
    {1, "long_ingap", 0, 1, 0x1.72919cf14d3bp+13, {0x1.a8b423b5e4a4bp-3, 0x1.c1301c2573b5bp+3, 0x1.3037792be78bfp-2, 0x1.13ca505928f31p-1, 0x1.02d5c22bc6e84p-1, 0x1.f4556773f5a65p+1, 0x1.cbafd99997cb1p-11, 0x1.56f209b13eacbp-1, 0x1.01937da6444c6p-12, 0x1.fffffffffef7ep-1}},
    {1, "long_ingap", 1, 1, 0x1.72521e10c549bp+13, {0x1.a9364ec44370ap-3, 0x1.c1786447c9a2ap+3, 0x1.30d8907585f9ep-2, 0x1.14a1839b53437p-1, 0x1.03a1667fc169fp-1, 0x1.f5541a1b43e72p+1, 0x1.cbafedfa1a0e8p-11, -0x1p+0, -0x1p+0, 0x1.fffffffffefb2p-1}},
    {1, "long_ingap", 2, 1, 0x1.3309bcf628abap+26, {0x1.0f3478eecc619p+8, 0x1.f80e8039e39fap+13, -0x1p+0, 0x1.2624a3f1500fbp+5, 0x1.ee1ed6263a179p+4, 0x1.0589f688c0e61p+11, -0x1p+0, 0x1.0b801c8e58327p+9, 0x1.ae077e6c1efa8p-9, 0x1.fffffffffffffp-1}},
    {1, "verylong", 0, 1, 0x1.2fa5efee817d3p+16, {0x1.461cd4cfcbfbcp-8, 0x1.77f7a34c2c822p+3, 0x1.49a540aa3fe71p-7, 0x1.64a2f6a4e6726p-11, 0x1.982ba23e28644p-7, 0x1.36cbf6ecbe4cdp-2, 0x1.a5678e7ce55eap-10, 0x1.0a2d2e644275cp-6, 0x1.584c73b70659fp-10, 0x1.ffffffcc245d1p-1}},
    {1, "verylong", 1, 1, 0x1.1eb968a21b6d9p+13, {0x1.13db7d14d593fp-6, 0x1.3c6c3ba1203ap+5, 0x1.18502448317ap-5, 0x1.07c6d3c1aac6dp-9, 0x1.7882c3b8f83d4p-6, 0x1.647767bbfb4d2p-1, 0x1.ef4663fd6783ep-9, -0x1p+0, -0x1p+0, 0x1.fffffffb5e193p-1}},
    {1, "verylong", 2, 1, 0x1.8899a5a7afdddp+16, {0x1.70e69ce0a40a8p-8, 0x1.9e6644ab2ae54p+3, -0x1p+0, 0x1.0214340c24db4p-10, 0x1.54afe2097cf6p-6, 0x1.84220433a4473p-2, -0x1p+0, 0x1.30f266fee9c4fp-6, 0x1.774aaeffa5c96p-10, 0x1.ffffffe401031p-1}},
};

// One event of the table, built and ready for FisherM; the callers change its epochs and run it.
struct EventRun {
    std::unique_ptr<source>   s = std::make_unique<source>();
    std::unique_ptr<lens>     l = std::make_unique<lens>();
    std::unique_ptr<astromet> as = std::make_unique<astromet>();
    std::unique_ptr<covarian> co = std::make_unique<covarian>();
    int ndw = 0, nL = 0, nR = 0;
    explicit EventRun(const FixtureEvent& ev)
    {
        setupStatic(*s, *l);
        l->tE = ev.tE;  l->t0 = ev.t0;  l->u0 = ev.u0;  l->piE = ev.piE;
        ndw = buildLightCurve(*s, *l, *as, nL, nR);
        l->tE = ev.tE;  l->t0 = ev.t0;  l->u0 = ev.u0;  l->piE = ev.piE;   // FisherM perturbs from whatever the structs hold
        s->xi = kXi;
    }
    // Keep the epochs for which keep(i) is true, compacting the buffers in place (stored order kept).
    template <class Pred> void keepOnly(Pred keep)
    {
        int w = 0;
        for (int i = 0; i < ndw; ++i) {
            if (not keep(i)) continue;
            l->timn[w] = l->timn[i]; l->magn[w] = l->magn[i]; l->errm[w] = l->errm[i];
            l->soux[w] = l->soux[i]; l->souy[w] = l->souy[i]; l->erra[w] = l->erra[i];
            l->tele[w] = l->tele[i]; l->filt[w] = l->filt[i]; l->rseas[w] = l->rseas[i]; l->rroll[w] = l->rroll[i];
            ++w;
        }
        ndw = w;
    }
    void addRubin(double tim, int filt) { appendEpoch(*s, *l, *as, tim, 0, ndw, -1, filt); }
    void run()
    {
        FisherM(*s, *l, *as, *co, ndw);
        ErrorCal(*co, *l, *s);
    }
};

bool sameDouble(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

// Every photometric output of two runs is identical to the last bit: the epoch counts, okA, condA and all of Era.
bool samePhotometry(const covarian& a, const covarian& b, const char* what, bool jointOnly)
{
    bool ok = true;
    for (int q = 0; q < NSURV; ++q) {
        if (jointOnly and q != SJOINT) continue;
        bool same = a.okA[q] == b.okA[q] and a.nepochA[q] == b.nepochA[q] and sameDouble(a.condA[q], b.condA[q]);
        for (int k = 0; k < Nx; ++k) same = same and sameDouble(a.Era[q][k], b.Era[q][k]);
        if (not same) {
            std::cerr << "FAIL [" << what << "] partition " << q << " differs from the reference in okA, nepochA, condA or Era\n";
            ok = false;
        }
    }
    return ok;
}

// With the prior on the blend fractions switched off, every kEvents event (all of Rubin's epochs in r) gives the forecast
// of the build before the bands existed, bit for bit: the per-filter blend/baseline arrays, the compact derivative rows,
// the sparse accumulation and the shared flux-pair function change no number for an event with one Rubin band.
bool checkPriorOffReproducesToday()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    int n = 0;
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        for (const auto& ev : kEvents) {
            EventRun r(ev);
            r.co->fbPriorSig = std::numeric_limits<double>::infinity();
            r.run();
            for (const auto& ref : kPriorOffRef) {
                if (ref.espl != (model == SourceModel::ESPL ? 1 : 0) or std::strcmp(ref.name, ev.name) != 0) continue;
                const int q = ref.q;
                bool same = r.co->okA[q] == ref.ok and sameDouble(r.co->condA[q], ref.cond);
                for (int k = 0; k < 10; ++k) same = same and sameDouble(r.co->Era[q][k], ref.era[k]);
                for (int k = 10; k < Nx; ++k) same = same and r.co->Era[q][k] == -1.0;   // no other band, no sigma
                ++n;
                if (not same) {
                    std::cerr << "FAIL [prior-off] " << ev.name << (model == SourceModel::ESPL ? " (espl)" : " (pspl)")
                              << " partition " << q << " is not bit-identical to the pre-band forecast\n";
                    ok = false;
                }
            }
        }
    }
    vbmSetSourceModel(keep);
    if (n != 2 * NSURV * int(sizeof(kEvents) / sizeof(kEvents[0]))) {
        std::cerr << "FAIL [prior-off] compared " << n << " rows, not every event/partition/model\n";
        ok = false;
    }
    std::cout << "# prior-off reproduces the f9e9938 forecast bit for bit: " << n << " rows, " << (ok ? "ok\n" : "FAIL\n");
    return ok;
}

// Five u-band visits far from the peak (A - 1 < 1e-6): the band's baseline is measured, its blend fraction is not (the
// prior holds it at sigma = FB_PRIOR_SIG), and no physical parameter moves.
bool checkBandAtBaseline()
{
    bool ok = true;
    const auto& ev = kEvents[2];   // mid_inseason
    EventRun ref(ev), r(ev);
    ref.run();
    r.s->blend[0] = 0.40;  r.s->magb[0] = 21.0;   // the u band's own blend fraction and baseline
    for (int k = 0; k < 5; ++k) r.addRubin(ev.t0 + 40.0 * ev.tE + 0.5 * k, 0);
    r.run();
    const covarian& c = *r.co;
    if (c.nRubinBand[0] != 5) { std::cerr << "FAIL [band-at-baseline] nRubinBand[u] = " << c.nRubinBand[0] << ", not 5\n"; ok = false; }
    for (int q = 0; q < NSURV; ++q) {
        if (not c.okA[q] and ref.co->okA[q]) { std::cerr << "FAIL [band-at-baseline] partition " << q << " not characterised\n"; ok = false; continue; }
        if (not c.okA[q]) continue;
        for (const int k : {0, 1, 3, 4, 5}) {
            const double a = c.Era[q][k], b = ref.co->Era[q][k];
            if (not (std::fabs(a - b) <= 1e-6 * b)) {
                std::cerr << "FAIL [band-at-baseline] partition " << q << " sigma(" << kPhotNames[k] << ") " << a << " vs " << b
                          << " without the u visits\n";
                ok = false;
            }
        }
    }
    // sigma(f_S,u) is the prior's, sigma(mbs_u) the five visits' own: 0.02 mag / sqrt(5)
    for (const int q : {SJOINT, SRUBIN}) {
        if (not (std::fabs(c.Era[q][10] - FB_PRIOR_SIG) < 1e-3 * FB_PRIOR_SIG)) {
            std::cerr << "FAIL [band-at-baseline] partition " << q << " sigma(fb_u) = " << c.Era[q][10] << ", not the prior's\n";
            ok = false;
        }
        if (not (std::fabs(c.Era[q][11] - kErrMagRubin / std::sqrt(5.0)) < 1e-3 * kErrMagRubin)) {
            std::cerr << "FAIL [band-at-baseline] partition " << q << " sigma(mbs_u) = " << c.Era[q][11] << "\n";
            ok = false;
        }
    }
    if (c.Era[SROMAN][10] != -1.0 or c.Era[SROMAN][11] != -1.0) { std::cerr << "FAIL [band-at-baseline] Roman reports a u-band sigma\n"; ok = false; }
    std::cout << "# band at baseline: sigma(fb_u) " << c.Era[SJOINT][10] << ", sigma(mbs_u) " << c.Era[SJOINT][11]
              << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// Two g-band visits at the peak are too few to fit the band: they are left out of the photometric matrices and
// nothing about the photometric forecast changes. (The astrometric matrices leave the band out too: Task 3.)
bool checkFewBandVisits()
{
    const auto& ev = kEvents[2];
    EventRun ref(ev), r(ev);
    ref.run();
    r.addRubin(ev.t0 - 0.1 * ev.tE, 1);
    r.addRubin(ev.t0 + 0.1 * ev.tE, 1);
    r.run();
    const bool ok = samePhotometry(*r.co, *ref.co, "few-band-visits", false) and r.co->nRubinBand[1] == 0;
    std::cout << "# two g visits at the peak leave the photometric forecast bit for bit unchanged" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// Rubin with only two r and two i visits has no band to fit: its partition is empty and the joint forecast is Roman's
// alone, bit for bit.
bool checkRubinBandsAllSparse()
{
    bool ok = true;
    const auto& ev = kEvents[2];
    EventRun ref(ev), r(ev);
    ref.keepOnly([&](int i) { return ref.l->tele[i] == 1; });
    ref.run();
    r.keepOnly([&](int i) { return r.l->tele[i] == 1; });
    r.addRubin(ev.t0 - 3.0, 2);  r.addRubin(ev.t0 + 3.0, 2);
    r.addRubin(ev.t0 - 2.0, 3);  r.addRubin(ev.t0 + 2.0, 3);
    r.run();
    if (r.co->okA[SRUBIN] != 0 or r.co->nepochA[SRUBIN] != 0) {
        std::cerr << "FAIL [rubin-all-sparse] Rubin partition okA " << r.co->okA[SRUBIN] << ", nepochA " << r.co->nepochA[SRUBIN] << "\n";
        ok = false;
    }
    if (nRubinBandsFitted(*r.co) != 0) { std::cerr << "FAIL [rubin-all-sparse] a band was fitted\n"; ok = false; }
    ok = samePhotometry(*r.co, *ref.co, "rubin-all-sparse", true) and ok;
    std::cout << "# 2 r + 2 i Rubin visits: Rubin out, joint = Roman alone" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// Rubin's visits all in i (at the bump) plus Roman: the joint and Rubin fits hold i's pair and not r's. The r blend
// fraction is then not measured, and must read as the -1 sentinel everywhere, in Era and in the relative error resu[2]
// (a negative value there would look like a measurement to a sum or a histogram).
bool checkROmittedSentinel()
{
    bool ok = true;
    EventRun r(kEvents[2]);
    r.s->blend[3] = 0.50;  r.s->magb[3] = 19.9;
    for (int i = 0; i < r.ndw; ++i) if (r.l->tele[i] == 0) r.l->filt[i] = 3;
    r.run();
    const covarian& c = *r.co;
    if (not (c.okA[SJOINT] and c.okA[SRUBIN])) { std::cerr << "FAIL [r-omitted] okA " << c.okA[SJOINT] << c.okA[SRUBIN] << "\n"; ok = false; }
    for (const int q : {SJOINT, SRUBIN}) {
        if (c.Era[q][2] != -1.0) { std::cerr << "FAIL [r-omitted] partition " << q << " Era[2] = " << c.Era[q][2] << ", not -1\n"; ok = false; }
        if (not (c.Era[q][14] > 0.0)) { std::cerr << "FAIL [r-omitted] partition " << q << " sigma(fb_i) = " << c.Era[q][14] << "\n"; ok = false; }
    }
    if (c.resu[2] != -1.0) { std::cerr << "FAIL [r-omitted] resu[2] = " << c.resu[2] << ", not the -1 sentinel\n"; ok = false; }
    std::cout << "# r absent from the fit: Era[2] and resu[2] are the -1 sentinel, sigma(fb_i) " << c.Era[SJOINT][14]
              << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// A chromatic event: Rubin's visits cycle r, i, z with their own blend fractions (0.3, 0.5, 0.6). The forecast has a
// pair per band, the joint is no worse than either survey for every parameter, F[J] = F[L] + F[R] up to the priors,
// F F^-1 = 1, and the bands' own sigmas exist.
bool checkChromaticEvent()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        EventRun r(kEvents[2]);
        makeChromatic(*r.s, *r.l, r.ndw);
        r.run();
        const covarian& c = *r.co;
        const char* name = (model == SourceModel::ESPL) ? "chromatic (espl)" : "chromatic (pspl)";
        if (not (c.okA[SJOINT] and c.okA[SRUBIN] and c.okA[SROMAN])) {
            std::cerr << "FAIL [" << name << "] okA " << c.okA[SJOINT] << c.okA[SRUBIN] << c.okA[SROMAN] << "\n";
            ok = false;
        }
        if (nRubinBandsFitted(c) != 3) { std::cerr << "FAIL [" << name << "] bands fitted " << nRubinBandsFitted(c) << ", not 3\n"; ok = false; }
        for (const int q : {SJOINT, SRUBIN})
            for (const int k : {2, 6, 14, 15, 16, 17})
                if (not (c.Era[q][k] > 0.0)) { std::cerr << "FAIL [" << name << "] partition " << q << " sigma(" << kPhotNames[k] << ") = " << c.Era[q][k] << "\n"; ok = false; }
        for (const int k : {10, 11, 12, 13, 18, 19})
            if (c.Era[SJOINT][k] != -1.0) { std::cerr << "FAIL [" << name << "] sigma(" << kPhotNames[k] << ") for a band without visits\n"; ok = false; }
        for (const int q : {SRUBIN, SROMAN})
            for (const int k : {0, 1, 2, 3, 4, 5, 9}) {
                if (k == 2 and q == SROMAN) continue;
                const double j = c.Era[SJOINT][k], p = c.Era[q][k];
                if (j > 0.0 and p > 0.0 and j > p * (1.0 + 1e-9)) {
                    std::cerr << "FAIL [" << name << "] sigma(" << kPhotNames[k] << ") joint " << j << " > " << (q == SRUBIN ? "rubin " : "roman ") << p << "\n";
                    ok = false;
                }
            }
        // Each band's flux pair is in that band's own epochs only, so the information between two different bands'
        // pairs is exactly zero in every partition's stored matrix (a leak would put the other band's derivative here).
        for (int q = 0; q < NSURV; ++q) {
            static const int kBandPair[3][2] = {{2, 6}, {14, 15}, {16, 17}};   // r, i, z: fb, mbs
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < a; ++b)
                    for (const int ja : kBandPair[a]) for (const int jb : kBandPair[b]) {
                        const double v = gsl_matrix_get(c.inputA[q].get(), ja, jb), w = gsl_matrix_get(c.inputA[q].get(), jb, ja);
                        if (v != 0.0 or w != 0.0) {
                            std::cerr << "FAIL [" << name << "] partition " << q << " F[" << kPhotNames[ja] << "][" << kPhotNames[jb]
                                      << "] = " << v << ", not exactly 0 between different bands\n";
                            ok = false;
                        }
                    }
        }
        ok = checkAdditivity(c, name, Nx, true) and ok;
        ok = checkJointNoWorse(c, name, kPhotNames, Nx, true) and ok;
        ok = checkInverseRoundTrip(c, name, Nx, true) and ok;
        std::cout << "# " << name << ": sigma(fb_i) " << c.Era[SJOINT][14] << ", sigma(fb_z) " << c.Era[SJOINT][16]
                  << ", sigma(fb0) " << c.Era[SJOINT][2] << (ok ? "  -> ok\n" : "  -> FAIL\n");
    }
    vbmSetSourceModel(keep);
    return ok;
}

// ---------------------------------------------------------------------------------------------
// Astrometric dilution per band. Each Rubin visit's centroid is diluted with its own band's fractions, each Rubin band
// carries its own unknown neighbours' offset (all sharing one reference position), and the visits of a Rubin band with
// fewer than kMinTeleEpochs visits leave the astrometric matrix as they leave the photometric one.
// ---------------------------------------------------------------------------------------------

// The astrometric forecast (Erb, okB, condB and the three noise variants' ErbV, okBV, condBV of every partition) the
// 2ccd875 build printed for every kEvents event, at full precision (hexadecimal floats, exact), under the point source
// (espl 0) and the uniform disc (espl 1): every Rubin visit in r, where the astrometry has always diluted with r's
// fractions. Generated once from that build; never regenerated to make a check pass.
struct AstroRefV { int ok; double cond; double erb[Ny]; };
struct AstroRef { int espl; const char* name; int q; int ok; double cond; double erb[Ny]; AstroRefV v[NAVAR]; };
const AstroRef kAstroRef[] = {
    {0, "short_inseason", 0, 1, 0x1.2df999ee5c1b6p+6, {0x1.23e360f1d3efap-6, 0x1.4713e488df562p-12, 0x1.c2668e2500c75p-12, 0x1.f087c7409d503p+2}, {{1, 0x1.2df999ee5c1b6p+6, {0x1.23e360f1d3efap-6, 0x1.4713e488df562p-12, 0x1.c2668e2500c75p-12, 0x1.f087c7409d503p+2}}, {1, 0x1.5871468ac3224p+2, {0x1.80d85e3ff3a9ap-3, 0x1.984975649aef6p-10, 0x1.b1acbd2ddafcep-10, 0x1.127e44111496ap+5}}, {1, 0x1.c7964c71a0af8p+0, {0x1.3845337fe58b2p-2, 0x1.e272f6e8ba22cp-10, 0x1.df1eb6319e0bfp-10, 0x1.17e6523087af6p+6}}}},
    {0, "short_inseason", 1, 1, 0x1.4b23148eb5accp+5, {0x1.08dd668c7f492p+0, 0x1.1e120a77fcf37p-8, 0x1.3ef6872526fcdp-8, 0x1.721bcba7dae17p+10}, {{1, 0x1.4b23148eb5accp+5, {0x1.08dd668c7f492p+0, 0x1.1e120a77fcf37p-8, 0x1.3ef6872526fcdp-8, 0x1.721bcba7dae17p+10}}, {1, 0x1.4b23148eb5accp+5, {0x1.08dd668c7f492p+0, 0x1.1e120a77fcf37p-8, 0x1.3ef6872526fcdp-8, 0x1.721bcba7dae17p+10}}, {1, 0x1.4b23148eb5accp+5, {0x1.08dd668c7f492p+0, 0x1.1e120a77fcf37p-8, 0x1.3ef6872526fcdp-8, 0x1.721bcba7dae17p+10}}}},
    {0, "short_inseason", 2, 1, 0x1.46fc13388394ep+6, {0x1.245ce42a6d8e4p-6, 0x1.524591dbec603p-12, 0x1.d4458fda71102p-12, 0x1.0214d4bf976aap+3}, {{1, 0x1.46fc13388394ep+6, {0x1.245ce42a6d8e4p-6, 0x1.524591dbec603p-12, 0x1.d4458fda71102p-12, 0x1.0214d4bf976aap+3}}, {1, 0x1.4e1862c3d2e51p+5, {0x1.b218c5e6516f7p-3, 0x1.d301852397161p-9, 0x1.352c23e64b6edp-8, 0x1.4a8c21bb88cecp+6}}, {1, 0x1.d9d7be37c80d6p+3, {0x1.5c4b1631120bp-2, 0x1.e043f58fe0c4dp-8, 0x1.2354a366d9445p-7, 0x1.1d20a787d22f1p+7}}}},
    {0, "short_ingap", 0, 1, 0x1.16d6b4d20cb91p+3, {0x1.1153c59831c26p+0, 0x1.657951dfd3b89p-10, 0x1.91860738f85a2p-10, 0x1.a0d7d2c3b9addp+9}, {{1, 0x1.16d6b4d20cb91p+3, {0x1.1153c59831c26p+0, 0x1.657951dfd3b89p-10, 0x1.91860738f85a2p-10, 0x1.a0d7d2c3b9addp+9}}, {1, 0x1.16d6b4d20cb91p+3, {0x1.1153c59831c26p+0, 0x1.657951dfd3b89p-10, 0x1.91860738f85a2p-10, 0x1.a0d7d2c3b9addp+9}}, {1, 0x1.16d6b4d20cb91p+3, {0x1.1153c59831c26p+0, 0x1.657951dfd3b89p-10, 0x1.91860738f85a2p-10, 0x1.a0d7d2c3b9addp+9}}}},
    {0, "short_ingap", 1, 1, 0x1.16d6b4d20cb91p+3, {0x1.1153c59831c26p+0, 0x1.657951dfd3b89p-10, 0x1.91860738f85a2p-10, 0x1.a0d7d2c3b9addp+9}, {{1, 0x1.16d6b4d20cb91p+3, {0x1.1153c59831c26p+0, 0x1.657951dfd3b89p-10, 0x1.91860738f85a2p-10, 0x1.a0d7d2c3b9addp+9}}, {1, 0x1.16d6b4d20cb91p+3, {0x1.1153c59831c26p+0, 0x1.657951dfd3b89p-10, 0x1.91860738f85a2p-10, 0x1.a0d7d2c3b9addp+9}}, {1, 0x1.16d6b4d20cb91p+3, {0x1.1153c59831c26p+0, 0x1.657951dfd3b89p-10, 0x1.91860738f85a2p-10, 0x1.a0d7d2c3b9addp+9}}}},
    {0, "short_ingap", 2, 0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}, {{0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}, {0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}, {0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}}},
    {0, "mid_inseason", 0, 1, 0x1.03e5866af86b5p+2, {0x1.a5b6ef7cf62c9p-7, 0x1.008c9edcbd7b2p-15, 0x1.fc4f64a1874a3p-16, 0x1.3f80cf7c27aa1p-1}, {{1, 0x1.03e5866af86b5p+2, {0x1.a5b6ef7cf62c9p-7, 0x1.008c9edcbd7b2p-15, 0x1.fc4f64a1874a3p-16, 0x1.3f80cf7c27aa1p-1}}, {1, 0x1.211e821f546d7p+1, {0x1.3c168e8e5fap-3, 0x1.d8d870b091bfbp-14, 0x1.dcfd9aa556d2p-14, 0x1.1d7d21721af9fp+2}}, {1, 0x1.237406233d695p+1, {0x1.4743d75e14051p-2, 0x1.e5833982f6923p-14, 0x1.e6fe9a0836786p-14, 0x1.3f6228ede6426p+2}}}},
    {0, "mid_inseason", 1, 1, 0x1.4ff45bb1f86dbp+1, {0x1.ac0e2444402c7p-2, 0x1.f09f2571d0548p-14, 0x1.ee80cfbf4e474p-14, 0x1.48b2e3957eb33p+2}, {{1, 0x1.4ff45bb1f86dbp+1, {0x1.ac0e2444402c7p-2, 0x1.f09f2571d0548p-14, 0x1.ee80cfbf4e474p-14, 0x1.48b2e3957eb33p+2}}, {1, 0x1.4ff45bb1f86dbp+1, {0x1.ac0e2444402c7p-2, 0x1.f09f2571d0548p-14, 0x1.ee80cfbf4e474p-14, 0x1.48b2e3957eb33p+2}}, {1, 0x1.4ff45bb1f86dbp+1, {0x1.ac0e2444402c7p-2, 0x1.f09f2571d0548p-14, 0x1.ee80cfbf4e474p-14, 0x1.48b2e3957eb33p+2}}}},
    {0, "mid_inseason", 2, 1, 0x1.07f19b25add7fp+2, {0x1.a7c6b8f90b74bp-7, 0x1.0ad61752e8791p-15, 0x1.07d28c97a79f1p-15, 0x1.4298c985131ccp-1}, {{1, 0x1.07f19b25add7fp+2, {0x1.a7c6b8f90b74bp-7, 0x1.0ad61752e8791p-15, 0x1.07d28c97a79f1p-15, 0x1.4298c985131ccp-1}}, {1, 0x1.3d5c54a4df2afp+3, {0x1.00a3c94d685ccp-2, 0x1.c3d82d354e6e4p-10, 0x1.78fd3a912d317p-10, 0x1.7a758c5975e23p+3}}, {1, 0x1.755d2d987bc7ep+2, {0x1.503aeed79f6b8p-1, 0x1.25ba7467df75p-8, 0x1.107f7c13cb435p-8, 0x1.867dd35b98cddp+4}}}},
    {0, "mid_ingap", 0, 1, 0x1.4dbb68986be57p+2, {0x1.24dfcf8216c28p-1, 0x1.3ca29f3a4d4afp-13, 0x1.343cc1fc1c63ap-13, 0x1.788924b369e7cp+2}, {{1, 0x1.4dbb68986be57p+2, {0x1.24dfcf8216c28p-1, 0x1.3ca29f3a4d4afp-13, 0x1.343cc1fc1c63ap-13, 0x1.788924b369e7cp+2}}, {1, 0x1.4dbb68986be57p+2, {0x1.24dfcf8216c28p-1, 0x1.3ca29f3a4d4afp-13, 0x1.343cc1fc1c63ap-13, 0x1.788924b369e7cp+2}}, {1, 0x1.4dbb68986be57p+2, {0x1.24dfcf8216c28p-1, 0x1.3ca29f3a4d4afp-13, 0x1.343cc1fc1c63ap-13, 0x1.788924b369e7cp+2}}}},
    {0, "mid_ingap", 1, 1, 0x1.4dbb68986be57p+2, {0x1.24dfcf8216c28p-1, 0x1.3ca29f3a4d4afp-13, 0x1.343cc1fc1c63ap-13, 0x1.788924b369e7cp+2}, {{1, 0x1.4dbb68986be57p+2, {0x1.24dfcf8216c28p-1, 0x1.3ca29f3a4d4afp-13, 0x1.343cc1fc1c63ap-13, 0x1.788924b369e7cp+2}}, {1, 0x1.4dbb68986be57p+2, {0x1.24dfcf8216c28p-1, 0x1.3ca29f3a4d4afp-13, 0x1.343cc1fc1c63ap-13, 0x1.788924b369e7cp+2}}, {1, 0x1.4dbb68986be57p+2, {0x1.24dfcf8216c28p-1, 0x1.3ca29f3a4d4afp-13, 0x1.343cc1fc1c63ap-13, 0x1.788924b369e7cp+2}}}},
    {0, "mid_ingap", 2, 0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}, {{0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}, {0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}, {0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}}},
    {0, "long_inseason", 0, 1, 0x1.43262f6e0b06p+1, {0x1.32b72d2b619cbp-7, 0x1.9ba48f97e6d41p-21, 0x1.96f336ba7d108p-21, 0x1.095ac2ef4b9dp-5}, {{1, 0x1.43262f6e0b06p+1, {0x1.32b72d2b619cbp-7, 0x1.9ba48f97e6d41p-21, 0x1.96f336ba7d108p-21, 0x1.095ac2ef4b9dp-5}}, {1, 0x1.5d54fc549083p+1, {0x1.045fb15bc76fep-3, 0x1.45010b7ef88ddp-17, 0x1.400d4e1fb8d8ap-17, 0x1.a7bd56ffa4959p-2}}, {1, 0x1.981bf80c23619p+1, {0x1.cbde0007d2716p-3, 0x1.1aeb8e5c3e063p-16, 0x1.143c6fb99095ap-16, 0x1.4d1a9e2fad2e4p-1}}}},
    {0, "long_inseason", 1, 1, 0x1.b0b3d4d55aa65p+1, {0x1.fdf0590e19d4dp-3, 0x1.366168d80d743p-16, 0x1.2dd7bfa579861p-16, 0x1.65d3da8a51752p-1}, {{1, 0x1.b0b3d4d55aa65p+1, {0x1.fdf0590e19d4dp-3, 0x1.366168d80d743p-16, 0x1.2dd7bfa579861p-16, 0x1.65d3da8a51752p-1}}, {1, 0x1.b0b3d4d55aa65p+1, {0x1.fdf0590e19d4dp-3, 0x1.366168d80d743p-16, 0x1.2dd7bfa579861p-16, 0x1.65d3da8a51752p-1}}, {1, 0x1.b0b3d4d55aa65p+1, {0x1.fdf0590e19d4dp-3, 0x1.366168d80d743p-16, 0x1.2dd7bfa579861p-16, 0x1.65d3da8a51752p-1}}}},
    {0, "long_inseason", 2, 1, 0x1.4302575302a97p+1, {0x1.32f04457d1e3bp-7, 0x1.9c023f36dcf7p-21, 0x1.9751d2b38d863p-21, 0x1.09a4f28eb797cp-5}, {{1, 0x1.4302575302a97p+1, {0x1.32f04457d1e3bp-7, 0x1.9c023f36dcf7p-21, 0x1.9751d2b38d863p-21, 0x1.09a4f28eb797cp-5}}, {1, 0x1.442f5fae864dbp+1, {0x1.2fee391031bfbp-3, 0x1.7eb715b2d9bf9p-17, 0x1.7a2db0882db8fp-17, 0x1.0777ef08ce0e3p-1}}, {1, 0x1.3455c1b5936b5p+1, {0x1.0c4a4a798cd53p-1, 0x1.5c7578668601ep-15, 0x1.5919464efcdb5p-15, 0x1.caadaa2c745dcp+0}}}},
    {0, "long_ingap", 0, 1, 0x1.1ed8ecec9e37cp+6, {0x1.0f27cd267585fp-4, 0x1.b2b270977efap-18, 0x1.97129978a5616p-18, 0x1.f6852688a0ef5p-4}, {{1, 0x1.1ed8ecec9e37cp+6, {0x1.0f27cd267585fp-4, 0x1.b2b270977efap-18, 0x1.97129978a5616p-18, 0x1.f6852688a0ef5p-4}}, {1, 0x1.187b8f5466504p+3, {0x1.986e53ac16a6bp-3, 0x1.7ef921beeff0fp-16, 0x1.7c21b6213a6e1p-16, 0x1.a2bc8c0a66648p-2}}, {1, 0x1.0a35697495cd7p+2, {0x1.b772f796d250ep-3, 0x1.9c1e13df642a9p-16, 0x1.9be45df0a574fp-16, 0x1.bfdc11bc829edp-2}}}},
    {0, "long_ingap", 1, 1, 0x1.e3db65e8c47ecp+1, {0x1.bd47bba6238ebp-3, 0x1.a156585819caep-16, 0x1.a166e2f17ab4ap-16, 0x1.c4aa2acd92f53p-2}, {{1, 0x1.e3db65e8c47ecp+1, {0x1.bd47bba6238ebp-3, 0x1.a156585819caep-16, 0x1.a166e2f17ab4ap-16, 0x1.c4aa2acd92f53p-2}}, {1, 0x1.e3db65e8c47ecp+1, {0x1.bd47bba6238ebp-3, 0x1.a156585819caep-16, 0x1.a166e2f17ab4ap-16, 0x1.c4aa2acd92f53p-2}}, {1, 0x1.e3db65e8c47ecp+1, {0x1.bd47bba6238ebp-3, 0x1.a156585819caep-16, 0x1.a166e2f17ab4ap-16, 0x1.c4aa2acd92f53p-2}}}},
    {0, "long_ingap", 2, 1, 0x1.428a096d81c19p+6, {0x1.203551e1ace85p-4, 0x1.ca943899ebcfap-18, 0x1.ac986819de92ep-18, 0x1.065acd2283af7p-3}, {{1, 0x1.428a096d81c19p+6, {0x1.203551e1ace85p-4, 0x1.ca943899ebcfap-18, 0x1.ac986819de92ep-18, 0x1.065acd2283af7p-3}}, {1, 0x1.8fff2403b8912p+10, {0x1.c7c8844767e5dp+1, 0x1.4309d6a87b2cap-12, 0x1.1f67d7ded3edap-12, 0x1.2d976ff8b9e1fp+1}}, {1, 0x1.1910c68fc2b92p+10, {0x1.5d81fe717fb73p+3, 0x1.f4d3f7404b7bdp-11, 0x1.c1274a51b09a4p-11, 0x1.057981c2ad715p+3}}}},
    {0, "verylong", 0, 1, 0x1.ac5766d0bc249p+8, {0x1.58b84d99b5b87p-5, 0x1.67a4652885656p-18, 0x1.35d4a50e2a22ep-18, 0x1.444cbcfec6df4p-7}, {{1, 0x1.ac5766d0bc249p+8, {0x1.58b84d99b5b87p-5, 0x1.67a4652885656p-18, 0x1.35d4a50e2a22ep-18, 0x1.444cbcfec6df4p-7}}, {1, 0x1.a3c7df2aa4f9ep+6, {0x1.12909c956bf65p-2, 0x1.29935516c94c3p-15, 0x1.e7ce50c81df55p-16, 0x1.ded853bbb04f5p-4}}, {1, 0x1.6375e28a2b6e9p+5, {0x1.46e927857b3adp-2, 0x1.708d30d5c9bb7p-15, 0x1.132de4f0742ecp-15, 0x1.a7ebe6f6e62f2p-3}}}},
    {0, "verylong", 1, 1, 0x1.3b6ab9e28261ap+5, {0x1.5551a194d82c1p-2, 0x1.83244bef8787bp-15, 0x1.18082a08e2eadp-15, 0x1.dcd6c600f5dbap-3}, {{1, 0x1.3b6ab9e28261ap+5, {0x1.5551a194d82c1p-2, 0x1.83244bef8787bp-15, 0x1.18082a08e2eadp-15, 0x1.dcd6c600f5dbap-3}}, {1, 0x1.3b6ab9e28261ap+5, {0x1.5551a194d82c1p-2, 0x1.83244bef8787bp-15, 0x1.18082a08e2eadp-15, 0x1.dcd6c600f5dbap-3}}, {1, 0x1.3b6ab9e28261ap+5, {0x1.5551a194d82c1p-2, 0x1.83244bef8787bp-15, 0x1.18082a08e2eadp-15, 0x1.dcd6c600f5dbap-3}}}},
    {0, "verylong", 2, 1, 0x1.b51ee25344535p+8, {0x1.5c7be5af70898p-5, 0x1.6b7b56b619f44p-18, 0x1.393afd90b675fp-18, 0x1.457c9b002a7e6p-7}, {{1, 0x1.b51ee25344535p+8, {0x1.5c7be5af70898p-5, 0x1.6b7b56b619f44p-18, 0x1.393afd90b675fp-18, 0x1.457c9b002a7e6p-7}}, {1, 0x1.bc1eabb9fe1c4p+8, {0x1.43af63e2efbccp-1, 0x1.5224dcf85d61dp-14, 0x1.2358dbf9d3e02p-14, 0x1.3b21f7e39ce7ap-3}}, {1, 0x1.af3ab8038929p+8, {0x1.23f100708307fp+1, 0x1.314a698ade81dp-12, 0x1.06b67ec90db5dp-12, 0x1.1a46bc977a41dp-1}}}},
    {1, "short_inseason", 0, 1, 0x1.2df9fc360fa99p+6, {0x1.23e36e952d769p-6, 0x1.47140670aa31dp-12, 0x1.c266dab3631dcp-12, 0x1.f0881a5a4e0a1p+2}, {{1, 0x1.2df9fc360fa99p+6, {0x1.23e36e952d769p-6, 0x1.47140670aa31dp-12, 0x1.c266dab3631dcp-12, 0x1.f0881a5a4e0a1p+2}}, {1, 0x1.58716ebbc06efp+2, {0x1.80d871341df5ap-3, 0x1.98497a78a2a14p-10, 0x1.b1acd20e9d233p-10, 0x1.127e5e1cfcd8dp+5}}, {1, 0x1.c796299e8d3e5p+0, {0x1.38454ae46b173p-2, 0x1.e272f82a76f0ap-10, 0x1.df1ebab44ab41p-10, 0x1.17e69a4bc15a6p+6}}}},
    {1, "short_inseason", 1, 1, 0x1.4b22f94813506p+5, {0x1.08dd6d4b370bp+0, 0x1.1e120182509d9p-8, 0x1.3ef67b14c8f57p-8, 0x1.721bbd138456fp+10}, {{1, 0x1.4b22f94813506p+5, {0x1.08dd6d4b370bp+0, 0x1.1e120182509d9p-8, 0x1.3ef67b14c8f57p-8, 0x1.721bbd138456fp+10}}, {1, 0x1.4b22f94813506p+5, {0x1.08dd6d4b370bp+0, 0x1.1e120182509d9p-8, 0x1.3ef67b14c8f57p-8, 0x1.721bbd138456fp+10}}, {1, 0x1.4b22f94813506p+5, {0x1.08dd6d4b370bp+0, 0x1.1e120182509d9p-8, 0x1.3ef67b14c8f57p-8, 0x1.721bbd138456fp+10}}}},
    {1, "short_inseason", 2, 1, 0x1.46fc869e56c05p+6, {0x1.245cf1f32ac3fp-6, 0x1.5245b910ba613p-12, 0x1.d445e5e2ec545p-12, 0x1.0215037a86229p+3}, {{1, 0x1.46fc869e56c05p+6, {0x1.245cf1f32ac3fp-6, 0x1.5245b910ba613p-12, 0x1.d445e5e2ec545p-12, 0x1.0215037a86229p+3}}, {1, 0x1.4e19cc789dd8dp+5, {0x1.b218e01217a39p-3, 0x1.d30224f98866p-9, 0x1.352cc3f457f9ap-8, 0x1.4a8cddaf939d5p+6}}, {1, 0x1.d9db2315fbb95p+3, {0x1.5c4b2db211f01p-2, 0x1.e044db12a93bcp-8, 0x1.235580da3c5d2p-7, 0x1.1d21d8beac7d9p+7}}}},
    {1, "short_ingap", 0, 1, 0x1.16d6a2210e857p+3, {0x1.1153caac5c953p+0, 0x1.6579529d5d2f8p-10, 0x1.918602a1c0e58p-10, 0x1.a0d7c0566b756p+9}, {{1, 0x1.16d6a2210e857p+3, {0x1.1153caac5c953p+0, 0x1.6579529d5d2f8p-10, 0x1.918602a1c0e58p-10, 0x1.a0d7c0566b756p+9}}, {1, 0x1.16d6a2210e857p+3, {0x1.1153caac5c953p+0, 0x1.6579529d5d2f8p-10, 0x1.918602a1c0e58p-10, 0x1.a0d7c0566b756p+9}}, {1, 0x1.16d6a2210e857p+3, {0x1.1153caac5c953p+0, 0x1.6579529d5d2f8p-10, 0x1.918602a1c0e58p-10, 0x1.a0d7c0566b756p+9}}}},
    {1, "short_ingap", 1, 1, 0x1.16d6a2210e857p+3, {0x1.1153caac5c953p+0, 0x1.6579529d5d2f8p-10, 0x1.918602a1c0e58p-10, 0x1.a0d7c0566b756p+9}, {{1, 0x1.16d6a2210e857p+3, {0x1.1153caac5c953p+0, 0x1.6579529d5d2f8p-10, 0x1.918602a1c0e58p-10, 0x1.a0d7c0566b756p+9}}, {1, 0x1.16d6a2210e857p+3, {0x1.1153caac5c953p+0, 0x1.6579529d5d2f8p-10, 0x1.918602a1c0e58p-10, 0x1.a0d7c0566b756p+9}}, {1, 0x1.16d6a2210e857p+3, {0x1.1153caac5c953p+0, 0x1.6579529d5d2f8p-10, 0x1.918602a1c0e58p-10, 0x1.a0d7c0566b756p+9}}}},
    {1, "short_ingap", 2, 0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}, {{0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}, {0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}, {0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}}},
    {1, "mid_inseason", 0, 1, 0x1.03e581a2dba23p+2, {0x1.a5b710f2de538p-7, 0x1.008ca1deef1dap-15, 0x1.fc4f5b38cb78fp-16, 0x1.3f80d5d6ce5a8p-1}, {{1, 0x1.03e581a2dba23p+2, {0x1.a5b710f2de538p-7, 0x1.008ca1deef1dap-15, 0x1.fc4f5b38cb78fp-16, 0x1.3f80d5d6ce5a8p-1}}, {1, 0x1.211e7bd5c755dp+1, {0x1.3c16a75860c67p-3, 0x1.d8d871a8134c1p-14, 0x1.dcfd9b672e682p-14, 0x1.1d7d25bbe1d8dp+2}}, {1, 0x1.23740e26f1f51p+1, {0x1.4743ec673cf4p-2, 0x1.e5833bb824a5bp-14, 0x1.e6fe9ba1ad369p-14, 0x1.3f622c3e8ee0cp+2}}}},
    {1, "mid_inseason", 1, 1, 0x1.4ff45f16e3a08p+1, {0x1.ac0e3289a7d52p-2, 0x1.f09f26d8c6094p-14, 0x1.ee80d1ae64608p-14, 0x1.48b2e4ccf654bp+2}, {{1, 0x1.4ff45f16e3a08p+1, {0x1.ac0e3289a7d52p-2, 0x1.f09f26d8c6094p-14, 0x1.ee80d1ae64608p-14, 0x1.48b2e4ccf654bp+2}}, {1, 0x1.4ff45f16e3a08p+1, {0x1.ac0e3289a7d52p-2, 0x1.f09f26d8c6094p-14, 0x1.ee80d1ae64608p-14, 0x1.48b2e4ccf654bp+2}}, {1, 0x1.4ff45f16e3a08p+1, {0x1.ac0e3289a7d52p-2, 0x1.f09f26d8c6094p-14, 0x1.ee80d1ae64608p-14, 0x1.48b2e4ccf654bp+2}}}},
    {1, "mid_inseason", 2, 1, 0x1.07f1968a2638cp+2, {0x1.a7c6dae48cf8cp-7, 0x1.0ad61aaf597b3p-15, 0x1.07d28758cd3cap-15, 0x1.4298d0063177fp-1}, {{1, 0x1.07f1968a2638cp+2, {0x1.a7c6dae48cf8cp-7, 0x1.0ad61aaf597b3p-15, 0x1.07d28758cd3cap-15, 0x1.4298d0063177fp-1}}, {1, 0x1.3d5c7518a5bcp+3, {0x1.00a3e9839adcdp-2, 0x1.c3d846bedb409p-10, 0x1.78fd37f41c876p-10, 0x1.7a75ad693a0c3p+3}}, {1, 0x1.755da20c03405p+2, {0x1.503b304a2b0f2p-1, 0x1.25ba8abe4ffb8p-8, 0x1.107f7acfd9b53p-8, 0x1.867e119aecf8p+4}}}},
    {1, "mid_ingap", 0, 1, 0x1.4dbb67f2bed39p+2, {0x1.24dfd00c84a72p-1, 0x1.3ca29f57a2483p-13, 0x1.343cc1f417c63p-13, 0x1.788924771da48p+2}, {{1, 0x1.4dbb67f2bed39p+2, {0x1.24dfd00c84a72p-1, 0x1.3ca29f57a2483p-13, 0x1.343cc1f417c63p-13, 0x1.788924771da48p+2}}, {1, 0x1.4dbb67f2bed39p+2, {0x1.24dfd00c84a72p-1, 0x1.3ca29f57a2483p-13, 0x1.343cc1f417c63p-13, 0x1.788924771da48p+2}}, {1, 0x1.4dbb67f2bed39p+2, {0x1.24dfd00c84a72p-1, 0x1.3ca29f57a2483p-13, 0x1.343cc1f417c63p-13, 0x1.788924771da48p+2}}}},
    {1, "mid_ingap", 1, 1, 0x1.4dbb67f2bed39p+2, {0x1.24dfd00c84a72p-1, 0x1.3ca29f57a2483p-13, 0x1.343cc1f417c63p-13, 0x1.788924771da48p+2}, {{1, 0x1.4dbb67f2bed39p+2, {0x1.24dfd00c84a72p-1, 0x1.3ca29f57a2483p-13, 0x1.343cc1f417c63p-13, 0x1.788924771da48p+2}}, {1, 0x1.4dbb67f2bed39p+2, {0x1.24dfd00c84a72p-1, 0x1.3ca29f57a2483p-13, 0x1.343cc1f417c63p-13, 0x1.788924771da48p+2}}, {1, 0x1.4dbb67f2bed39p+2, {0x1.24dfd00c84a72p-1, 0x1.3ca29f57a2483p-13, 0x1.343cc1f417c63p-13, 0x1.788924771da48p+2}}}},
    {1, "mid_ingap", 2, 0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}, {{0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}, {0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}, {0, -0x1p+0, {-0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}}}},
    {1, "long_inseason", 0, 1, 0x1.43261868790b4p+1, {0x1.32b7378103d05p-7, 0x1.9ba49009b7acp-21, 0x1.96f338e1f58a9p-21, 0x1.095ac3f0d4276p-5}, {{1, 0x1.43261868790b4p+1, {0x1.32b7378103d05p-7, 0x1.9ba49009b7acp-21, 0x1.96f338e1f58a9p-21, 0x1.095ac3f0d4276p-5}}, {1, 0x1.5d54e7feae498p+1, {0x1.045fbac07ccfdp-3, 0x1.45010bf6201dap-17, 0x1.400d5085405bep-17, 0x1.a7bd5b197860fp-2}}, {1, 0x1.981be97884465p+1, {0x1.cbde0c375f408p-3, 0x1.1aeb8ec51123ap-16, 0x1.143c726a5233dp-16, 0x1.4d1aa4f8a0798p-1}}}},
    {1, "long_inseason", 1, 1, 0x1.b0b3c7732c474p+1, {0x1.fdf063d5eaca1p-3, 0x1.3661692e0d5ecp-16, 0x1.2dd7c29a2703dp-16, 0x1.65d3e2dad3143p-1}, {{1, 0x1.b0b3c7732c474p+1, {0x1.fdf063d5eaca1p-3, 0x1.3661692e0d5ecp-16, 0x1.2dd7c29a2703dp-16, 0x1.65d3e2dad3143p-1}}, {1, 0x1.b0b3c7732c474p+1, {0x1.fdf063d5eaca1p-3, 0x1.3661692e0d5ecp-16, 0x1.2dd7c29a2703dp-16, 0x1.65d3e2dad3143p-1}}, {1, 0x1.b0b3c7732c474p+1, {0x1.fdf063d5eaca1p-3, 0x1.3661692e0d5ecp-16, 0x1.2dd7c29a2703dp-16, 0x1.65d3e2dad3143p-1}}}},
    {1, "long_inseason", 2, 1, 0x1.4302404812f7bp+1, {0x1.32f04eb07cdcbp-7, 0x1.9c023fa89df5dp-21, 0x1.9751d4dabfec7p-21, 0x1.09a4f38daae85p-5}, {{1, 0x1.4302404812f7bp+1, {0x1.32f04eb07cdcbp-7, 0x1.9c023fa89df5dp-21, 0x1.9751d4dabfec7p-21, 0x1.09a4f38daae85p-5}}, {1, 0x1.442f47043afabp+1, {0x1.2fee45811bc71p-3, 0x1.7eb7161a15d59p-17, 0x1.7a2db2ef286d6p-17, 0x1.0777efa40e8aap-1}}, {1, 0x1.3455a4da23573p+1, {0x1.0c4a57a8ee5fp-1, 0x1.5c7578bc46023p-15, 0x1.591948e1f7e91p-15, 0x1.caadab27911f1p+0}}}},
    {1, "long_ingap", 0, 1, 0x1.1ed8ee1d502bfp+6, {0x1.0f27cdeb42049p-4, 0x1.b2b27204e812p-18, 0x1.97129a1331886p-18, 0x1.f68527d08ffa5p-4}, {{1, 0x1.1ed8ee1d502bfp+6, {0x1.0f27cdeb42049p-4, 0x1.b2b27204e812p-18, 0x1.97129a1331886p-18, 0x1.f68527d08ffa5p-4}}, {1, 0x1.187b8f64764f7p+3, {0x1.986e5ba68af2p-3, 0x1.7ef9271cc1bf7p-16, 0x1.7c21b6f436c72p-16, 0x1.a2bc8f49e7984p-2}}, {1, 0x1.0a35679b22a84p+2, {0x1.b77301bfc6519p-3, 0x1.9c1e170345e8dp-16, 0x1.9be45e6ac4533p-16, 0x1.bfdc15c7f7debp-2}}}},
    {1, "long_ingap", 1, 1, 0x1.e3db61cb44d2ep+1, {0x1.bd47c642a3527p-3, 0x1.a1565b1f57799p-16, 0x1.a166e35d27106p-16, 0x1.c4aa2efd47f78p-2}, {{1, 0x1.e3db61cb44d2ep+1, {0x1.bd47c642a3527p-3, 0x1.a1565b1f57799p-16, 0x1.a166e35d27106p-16, 0x1.c4aa2efd47f78p-2}}, {1, 0x1.e3db61cb44d2ep+1, {0x1.bd47c642a3527p-3, 0x1.a1565b1f57799p-16, 0x1.a166e35d27106p-16, 0x1.c4aa2efd47f78p-2}}, {1, 0x1.e3db61cb44d2ep+1, {0x1.bd47c642a3527p-3, 0x1.a1565b1f57799p-16, 0x1.a166e35d27106p-16, 0x1.c4aa2efd47f78p-2}}}},
    {1, "long_ingap", 2, 1, 0x1.428a09afc6206p+6, {0x1.2035521c2f07dp-4, 0x1.ca9438ff58d8ap-18, 0x1.ac98688941babp-18, 0x1.065acdae385aap-3}, {{1, 0x1.428a09afc6206p+6, {0x1.2035521c2f07dp-4, 0x1.ca9438ff58d8ap-18, 0x1.ac98688941babp-18, 0x1.065acdae385aap-3}}, {1, 0x1.8fff243f5c851p+10, {0x1.c7c884761a333p+1, 0x1.4309d6c123715p-12, 0x1.1f67d7f352575p-12, 0x1.2d9770071ca27p+1}}, {1, 0x1.1910c6b439be9p+10, {0x1.5d81fe920b045p+3, 0x1.f4d3f7630022dp-11, 0x1.c1274a6bf9881p-11, 0x1.057981cde86dcp+3}}}},
    {1, "verylong", 0, 1, 0x1.ac5875204a108p+8, {0x1.58b8c6f676531p-5, 0x1.67a5064a09211p-18, 0x1.35d4f4ed711p-18, 0x1.444c43a46bf6p-7}, {{1, 0x1.ac5875204a108p+8, {0x1.58b8c6f676531p-5, 0x1.67a5064a09211p-18, 0x1.35d4f4ed711p-18, 0x1.444c43a46bf6p-7}}, {1, 0x1.a3c899df8b45bp+6, {0x1.1290e7a6325eep-2, 0x1.2993ba3ea1e85p-15, 0x1.e7ce8fb3c4f9dp-16, 0x1.ded7fc2ce987dp-4}}, {1, 0x1.6376b082fabb3p+5, {0x1.46e9a74da815dp-2, 0x1.708db55639deep-15, 0x1.132e061032358p-15, 0x1.a7ebbb6a24bbdp-3}}}},
    {1, "verylong", 1, 1, 0x1.3b6b90f533c29p+5, {0x1.55523d5be851dp-2, 0x1.8324e3c70b54bp-15, 0x1.18084de32a3b4p-15, 0x1.dcd6a571aded2p-3}, {{1, 0x1.3b6b90f533c29p+5, {0x1.55523d5be851dp-2, 0x1.8324e3c70b54bp-15, 0x1.18084de32a3b4p-15, 0x1.dcd6a571aded2p-3}}, {1, 0x1.3b6b90f533c29p+5, {0x1.55523d5be851dp-2, 0x1.8324e3c70b54bp-15, 0x1.18084de32a3b4p-15, 0x1.dcd6a571aded2p-3}}, {1, 0x1.3b6b90f533c29p+5, {0x1.55523d5be851dp-2, 0x1.8324e3c70b54bp-15, 0x1.18084de32a3b4p-15, 0x1.dcd6a571aded2p-3}}}},
    {1, "verylong", 2, 1, 0x1.b51ff93c092ffp+8, {0x1.5c7c618c30dbbp-5, 0x1.6b7bfaddb4acap-18, 0x1.393b4f71244d5p-18, 0x1.457c1ff4209d5p-7}, {{1, 0x1.b51ff93c092ffp+8, {0x1.5c7c618c30dbbp-5, 0x1.6b7bfaddb4acap-18, 0x1.393b4f71244d5p-18, 0x1.457c1ff4209d5p-7}}, {1, 0x1.bc1fec482d523p+8, {0x1.43afe1e619a2dp-1, 0x1.5225860ad6158p-14, 0x1.23592f1ee096dp-14, 0x1.3b2178f56a202p-3}}, {1, 0x1.af3bf5c0bd72bp+8, {0x1.23f1746d4098ep+1, 0x1.314b02ca2bf79p-12, 0x1.06b6cc8cc7e4fp-12, 0x1.1a4650a171c17p-1}}}},
};

// Every astrometric output of two runs is identical to the last bit.
bool sameAstrometry(const covarian& a, const covarian& b, const char* what)
{
    bool ok = true;
    for (int q = 0; q < NSURV; ++q) {
        bool same = a.okB[q] == b.okB[q] and sameDouble(a.condB[q], b.condB[q]);
        for (int k = 0; k < Ny; ++k) same = same and sameDouble(a.Erb[q][k], b.Erb[q][k]);
        for (int v = 0; v < NAVAR; ++v) {
            same = same and a.okBV[v][q] == b.okBV[v][q] and sameDouble(a.condBV[v][q], b.condBV[v][q]);
            for (int k = 0; k < Ny; ++k) same = same and sameDouble(a.ErbV[v][q][k], b.ErbV[v][q][k]);
            for (int jm = 0; jm < Ny * Ny; ++jm) same = same and sameDouble(a.FBV[v][q][jm], b.FBV[v][q][jm]);
        }
        if (not same) {
            std::cerr << "FAIL [" << what << "] partition " << q << " differs from the reference in an astrometric output\n";
            ok = false;
        }
    }
    return ok;
}

// With every Rubin visit in r (every kEvents event) the astrometric forecast is the one of the build before the
// astrometry was diluted per band, bit for bit: one band takes today's path.
bool checkAstroAllRUnchanged()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    int n = 0;
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        for (const auto& ev : kEvents) {
            EventRun r(ev);
            r.run();
            for (const auto& ref : kAstroRef) {
                if (ref.espl != (model == SourceModel::ESPL ? 1 : 0) or std::strcmp(ref.name, ev.name) != 0) continue;
                const int q = ref.q;
                const covarian& c = *r.co;
                bool same = c.okB[q] == ref.ok and sameDouble(c.condB[q], ref.cond);
                for (int k = 0; k < Ny; ++k) same = same and sameDouble(c.Erb[q][k], ref.erb[k]);
                for (int v = 0; v < NAVAR; ++v) {
                    same = same and c.okBV[v][q] == ref.v[v].ok and sameDouble(c.condBV[v][q], ref.v[v].cond);
                    for (int k = 0; k < Ny; ++k) same = same and sameDouble(c.ErbV[v][q][k], ref.v[v].erb[k]);
                }
                ++n;
                if (not same) {
                    std::cerr << "FAIL [astro-all-r] " << ev.name << (model == SourceModel::ESPL ? " (espl)" : " (pspl)")
                              << " partition " << q << " is not bit-identical to the 2ccd875 astrometric forecast\n";
                    ok = false;
                }
            }
        }
    }
    vbmSetSourceModel(keep);
    if (n != 2 * NSURV * int(sizeof(kEvents) / sizeof(kEvents[0]))) {
        std::cerr << "FAIL [astro-all-r] compared " << n << " rows, not every event/partition/model\n";
        ok = false;
    }
    std::cout << "# astrometry with every Rubin visit in r reproduces the 2ccd875 forecast bit for bit: " << n << " rows, "
              << (ok ? "ok\n" : "FAIL\n");
    return ok;
}

// Rubin's visits relabelled cycling through `cycle` (in stored order), the bands' fractions as the caller set them.
void relabelRubin(EventRun& r, const std::vector<int>& cycle)
{
    int c = 0;
    for (int i = 0; i < r.ndw; ++i) if (r.l->tele[i] == 0) r.l->filt[i] = cycle[c++ % cycle.size()];
}

// i and z set to r's fractions and baseline (the lens's share too).
void copyRBandToIZ(source& s)
{
    for (const int b : {3, 4}) { s.blend[b] = s.blend[2];  s.magb[b] = s.magb[2];  s.fLensF[b] = s.fLensF[2]; }
}

// The astrometric part of the standard checks (F[J] = F[L] + F[R], sigma_J <= single) on one finished run.
bool astroStandardChecks(const covarian& c, const char* name)
{
    bool ok = checkAdditivity(c, name, Ny, false);
    ok = checkJointNoWorse(c, name, kAstNames, Ny, false) and ok;
    return ok;
}

// Splitting Rubin's visits over three bands with IDENTICAL fractions changes nothing in what the visits say about the
// centroid, but each band now carries an offset of its own: more nuisance parameters cannot add astrometric
// information, and (decision a: one offset per band, in FisherM) with neighbours present they cost some. With no
// neighbours (f_N = 0 in every band) there is no offset to marginalise and the split changes nothing.
bool checkAstroBandSplitNoGain()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        EventRun a(kEvents[2]), b(kEvents[2]);
        copyRBandToIZ(*b.s);
        relabelRubin(b, {2, 3, 4});
        a.run();  b.run();
        for (const int q : {SRUBIN, SJOINT})
            for (int k = 0; k < Ny; ++k) {
                const double s1 = a.co->Erb[q][k], s3 = b.co->Erb[q][k];
                if (not (s3 >= s1 * (1.0 - 1e-12))) {
                    std::cerr << "FAIL [astro-band-split] partition " << q << " sigma[" << k << "]: three bands " << s3
                              << " < one band " << s1 << "\n";
                    ok = false;
                }
            }
        // One offset per band, not one shared by the bands: three offsets cost Rubin's sigma(thetaE) strictly.
        const double g1 = a.co->Erb[SRUBIN][0], g3 = b.co->Erb[SRUBIN][0];
        if (not (g3 > g1 * (1.0 + 1e-8))) {
            std::cerr << "FAIL [astro-band-split] Rubin sigma(thetaE) with three offsets " << g3
                      << " does not exceed the one with one offset " << g1 << " by 1e-8 relative: the bands share an offset\n";
            ok = false;
        }
        ok = astroStandardChecks(*b.co, "astro-band-split") and ok;
        std::cout << "# astro band split, identical fractions" << (model == SourceModel::ESPL ? " (espl)" : " (pspl)")
                  << ": sigma(thetaE) Rubin " << g1 << " (r) -> " << g3 << " (r,i,z), joint "
                  << a.co->Erb[SJOINT][0] << " -> " << b.co->Erb[SJOINT][0]
                  << ";  sigma(piE) Rubin " << a.co->Erb[SRUBIN][3] << " -> " << b.co->Erb[SRUBIN][3] << "\n";

        // No neighbours at all: f_S = 1 in r, i and z (dark lens): h = 0, no offset, the split changes only the order of sums.
        EventRun c(kEvents[2]), d(kEvents[2]);
        for (EventRun* e : {&c, &d}) {
            for (const int band : {2, 3, 4}) { e->s->blend[band] = 1.0;  e->s->fLensF[band] = 0.0; }
            e->s->fb[0] = 1.0;
        }
        relabelRubin(d, {2, 3, 4});
        c.run();  d.run();
        double worst = 0.0;
        for (const int q : {SRUBIN, SJOINT})
            for (int k = 0; k < Ny; ++k)
                worst = std::max(worst, std::fabs(d.co->Erb[q][k] - c.co->Erb[q][k]) / c.co->Erb[q][k]);
        if (not (worst < 1e-10)) {
            std::cerr << "FAIL [astro-band-split] without neighbours the split moves sigma by " << worst << " relative\n";
            ok = false;
        }
        std::cout << "#   f_N = 0 in every band: split vs all-r, worst relative difference " << worst << "\n";
    }
    vbmSetSourceModel(keep);
    std::cout << "#   no sigma smaller with three offsets, Rubin sigma(thetaE) strictly larger" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// Every Rubin visit relabelled i, with i's fractions equal to r's: one band, the same sums in the same order, every
// astrometric output bit for bit the all-r event's.
bool checkAstroSingleOtherBand()
{
    EventRun a(kEvents[2]), b(kEvents[2]);
    copyRBandToIZ(*b.s);
    relabelRubin(b, {3});
    a.run();  b.run();
    const bool ok = sameAstrometry(*b.co, *a.co, "astro-single-i");
    std::cout << "# every Rubin visit in i (fractions as r's): astrometry bit for bit the all-r event's" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// The neighbours' offset is an unknown of the fit and its true value cancels: each visit's centroid minus the
// neighbours' share (with that visit's band's f_S and f_N) is independent of blendOff. A share built with another
// band's (or the telescope's) fractions would leave blendOff in the derivatives. Also with a luminous lens whose
// share differs by band (the lens's light enters f_N).
bool checkAstroShareCancels()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        for (const bool luminous : {false, true}) {
            const char* name = luminous ? "astro-share-cancels (luminous)" : "astro-share-cancels";
            EventRun a(kEvents[2]), b(kEvents[2]);
            for (EventRun* e : {&a, &b}) {
                makeChromatic(*e->s, *e->l, e->ndw);
                if (luminous) {
                    e->s->fLensF[2] = 0.10;  e->s->fLensF[3] = 0.05;  e->s->fLensF[4] = 0.03;
                    e->s->fLens[0] = e->s->fLensF[2];
                }
            }
            b.s->blendOff[0] = {3.0, -2.0};                //Rubin's neighbours 3, -2 mas from the source
            a.run();  b.run();
            double worst = 0.0;
            const auto cmp = [&](double x, double y, double scale) { worst = std::max(worst, std::fabs(x - y) / scale); };
            for (int q = 0; q < NSURV; ++q) {
                for (int k = 0; k < Ny; ++k) cmp(a.co->Erb[q][k], b.co->Erb[q][k], std::fabs(a.co->Erb[q][k]) + 1e-300);
                for (int v = 0; v < NAVAR; ++v) {
                    for (int k = 0; k < Ny; ++k)
                        cmp(a.co->ErbV[v][q][k], b.co->ErbV[v][q][k], std::fabs(a.co->ErbV[v][q][k]) + 1e-300);
                    for (int j = 0; j < Ny; ++j)
                        for (int m = 0; m < Ny; ++m)
                            cmp(a.co->FBV[v][q][j * Ny + m], b.co->FBV[v][q][j * Ny + m],
                                std::sqrt(std::fabs(a.co->FBV[v][q][j * Ny + j] * a.co->FBV[v][q][m * Ny + m])) + 1e-300);
                }
            }
            if (not (worst < 1e-9)) {
                std::cerr << "FAIL [" << name << "] blendOff moves an astrometric output by " << worst << " relative\n";
                ok = false;
            }
            ok = astroStandardChecks(*b.co, name) and ok;
            std::cout << "# " << name << (model == SourceModel::ESPL ? " (espl)" : " (pspl)")
                      << ": blendOff (3, -2) mas vs 0, worst relative difference of any astrometric output " << worst << "\n";
        }
    }
    vbmSetSourceModel(keep);
    std::cout << "#   the neighbours' share cancels" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// The chromatic event (f_S 0.3 / 0.5 / 0.6 in r / i / z, dark lens) against the same event whose i and z visits are
// diluted like r's: i and z are less diluted, so Rubin's sigma(thetaE) is lower with the true fractions. Expected, not a
// theorem (the offset's derivative changes with the fractions as well).
bool checkAstroChromaticDilution()
{
    bool ok = true;
    EventRun a(kEvents[2]), b(kEvents[2]);
    makeChromatic(*a.s, *a.l, a.ndw);
    makeChromatic(*b.s, *b.l, b.ndw);
    copyRBandToIZ(*b.s);                       // the labels stay r, i, z; i and z are diluted like r
    a.run();  b.run();
    ok = astroStandardChecks(*a.co, "astro-chromatic") and ok;
    const double sTrue = a.co->Erb[SRUBIN][0], sSame = b.co->Erb[SRUBIN][0];
    if (not (sTrue < sSame)) {
        std::cerr << "FAIL [astro-chromatic] sigma(thetaE) Rubin with the true fractions " << sTrue
                  << " is not below the one with r's fractions in every band " << sSame << "\n";
        ok = false;
    }
    std::cout << "# astro chromatic dilution: sigma(thetaE) Rubin " << sTrue << " (true fractions) vs " << sSame
              << " (r's in every band); joint " << a.co->Erb[SJOINT][0] << " vs " << b.co->Erb[SJOINT][0]
              << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// Two g visits at the peak: the band has too few visits for the photometry and leaves the astrometry too, every
// astrometric output bit for bit as without them.
bool checkAstroSparseBandDropped()
{
    const auto& ev = kEvents[2];
    EventRun ref(ev), r(ev);
    ref.run();
    r.addRubin(ev.t0 - 0.1 * ev.tE, 1);
    r.addRubin(ev.t0 + 0.1 * ev.tE, 1);
    r.run();
    const bool ok = sameAstrometry(*r.co, *ref.co, "astro-sparse-band");
    std::cout << "# two g visits at the peak leave the astrometric forecast bit for bit unchanged" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// Rubin with only two r and two i visits has no band in the astrometry either: its partition is empty and the joint
// astrometric forecast is Roman's, bit for bit (the matrices too, in every noise variant).
bool checkAstroRubinAllSparse()
{
    bool ok = true;
    const auto& ev = kEvents[2];
    EventRun r(ev);
    r.keepOnly([&](int i) { return r.l->tele[i] == 1; });
    r.addRubin(ev.t0 - 3.0, 2);  r.addRubin(ev.t0 + 3.0, 2);
    r.addRubin(ev.t0 - 2.0, 3);  r.addRubin(ev.t0 + 2.0, 3);
    r.run();
    const covarian& c = *r.co;
    if (c.okB[SRUBIN] != 0) { std::cerr << "FAIL [astro-all-sparse] okB[Rubin] = " << int(c.okB[SRUBIN]) << ", not 0\n"; ok = false; }
    for (int k = 0; k < Ny; ++k)
        if (not sameDouble(c.Erb[SJOINT][k], c.Erb[SROMAN][k])) {
            std::cerr << "FAIL [astro-all-sparse] Erb[joint][" << k << "] = " << c.Erb[SJOINT][k] << " != Erb[roman] = " << c.Erb[SROMAN][k] << "\n";
            ok = false;
        }
    for (int v = 0; v < NAVAR; ++v)
        for (int jm = 0; jm < Ny * Ny; ++jm)
            if (not sameDouble(c.FBV[v][SJOINT][jm], c.FBV[v][SROMAN][jm])) {
                std::cerr << "FAIL [astro-all-sparse] FBV[" << v << "] joint != Roman at " << jm << "\n";
                ok = false;
            }
    std::cout << "# 2 r + 2 i Rubin visits: Rubin out of the astrometry, joint = Roman alone" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// ---- the closed form for several Rubin bands, on synthetic sums ----
using fisher_detail::AstSums;
constexpr int kD = fisher_detail::kAstD;

// One band's sums: n epochs, derivative vectors a smooth deterministic function of the epoch, offset derivative h(m).
template <class H> AstSums syntheticBand(int n, double shift, H hOf)
{
    AstSums k;
    for (int m = 0; m < n; ++m) {
        std::array<double, kD> dx{}, dy{};
        for (int j = 0; j < Ny; ++j) {
            dx[j] = (1.0 + 0.3 * j) * std::sin(0.37 * m + 1.1 * j + shift) + 0.2 * std::cos(0.11 * m * (j + 1));
            dy[j] = (1.0 - 0.1 * j) * std::cos(0.29 * m + 0.7 * j + 2.0 * shift) - 0.1 * std::sin(0.05 * m * (j + 2));
        }
        const double h = hOf(m);
        dx[Ny] = h;  dy[Ny] = 0.8 * h;                       // y: a different, still positive, offset derivative
        k.add(1.0 + 0.3 * std::cos(0.9 * m + shift), dx, dy);
    }
    return k;
}

// Solve C X = Bt (n x n, n x Ny) by Gauss elimination with partial pivoting; X returned in place of Bt.
void gaussSolve(std::vector<double>& C, std::vector<double>& Bt, int n)
{
    for (int p = 0; p < n; ++p) {
        int piv = p;
        for (int i = p + 1; i < n; ++i) if (std::fabs(C[i * n + p]) > std::fabs(C[piv * n + p])) piv = i;
        for (int j = 0; j < n; ++j) std::swap(C[p * n + j], C[piv * n + j]);
        for (int j = 0; j < Ny; ++j) std::swap(Bt[p * Ny + j], Bt[piv * Ny + j]);
        for (int i = p + 1; i < n; ++i) {
            const double f = C[i * n + p] / C[p * n + p];
            for (int j = p; j < n; ++j) C[i * n + j] -= f * C[p * n + j];
            for (int j = 0; j < Ny; ++j) Bt[i * Ny + j] -= f * Bt[p * Ny + j];
        }
    }
    for (int p = n - 1; p >= 0; --p)
        for (int j = 0; j < Ny; ++j) {
            double v = Bt[p * Ny + j];
            for (int m = p + 1; m < n; ++m) v -= C[p * n + m] * Bt[m * Ny + j];
            Bt[p * Ny + j] = v / C[p * n + p];
        }
}

// The brute force: per coordinate the full information matrix on (the four physical parameters, one offset for every
// band whose offset derivative is not identically zero, ONE reference position), the unknowns beyond the physical ones
// marginalised by a direct solve. Valid when that block is not singular (some band's h not constant).
std::array<double, Ny * Ny> bruteForceRubin(const std::array<AstSums, kRubinBands>& bands)
{
    std::array<double, Ny * Ny> out{};
    for (int cc = 0; cc < 2; ++cc) {
        std::vector<int> off;
        for (int b = 0; b < kRubinBands; ++b) if (bands[b].Sw > 0.0 and bands[b].Swdd[cc][Ny * kD + Ny] > 0.0) off.push_back(b);
        const int n = int(off.size()) + 1;                   // the offsets, then the reference
        std::vector<double> C(n * n, 0.0), Bt(n * Ny, 0.0);
        std::array<double, Ny * Ny> A{};
        for (int b = 0; b < kRubinBands; ++b) {
            if (bands[b].Sw <= 0.0) continue;
            const auto& S = bands[b];
            for (int j = 0; j < Ny; ++j) {
                for (int m = 0; m < Ny; ++m) A[j * Ny + m] += S.Swdd[cc][j * kD + m];
                Bt[(n - 1) * Ny + j] += S.Swd[cc][j];                    // reference row: sum w d_j
            }
            C[(n - 1) * n + (n - 1)] += S.Sw;
        }
        for (int a = 0; a < int(off.size()); ++a) {
            const auto& S = bands[off[a]];
            for (int j = 0; j < Ny; ++j) Bt[a * Ny + j] = S.Swdd[cc][j * kD + Ny];
            C[a * n + a] = S.Swdd[cc][Ny * kD + Ny];
            C[a * n + (n - 1)] = C[(n - 1) * n + a] = S.Swd[cc][Ny];
        }
        std::vector<double> X = Bt;
        gaussSolve(C, X, n);
        for (int j = 0; j < Ny; ++j)
            for (int m = 0; m < Ny; ++m) {
                double v = A[j * Ny + m];
                for (int p = 0; p < n; ++p) v -= Bt[p * Ny + j] * X[p * Ny + m];
                out[j * Ny + m] += v;
            }
    }
    return out;
}

// Today's route for one band (Group at sigma_c = 0, marginalised, then the Schur complement of the offset with the
// raw-h guard), reimplemented as the reference of the one-band case.
std::array<double, Ny * Ny> todayRoute(const AstSums& k)
{
    std::array<double, Ny * Ny> out{};
    for (int cc = 0; cc < 2; ++cc) {
        std::array<double, kD * kD> G{};
        for (int jm = 0; jm < kD * kD; ++jm)
            G[jm] = k.Swdd[cc][jm] - (k.Sw > 0.0 ? k.Swd[cc][jm / kD] * k.Swd[cc][jm % kD] / k.Sw : 0.0);
        const double raw = k.Swdd[cc][Ny * kD + Ny], hh = G[Ny * kD + Ny];
        const bool fit = AST_BLEND_OFFSET_FREE and raw > 0.0 and hh > 1.0e-9 * raw;
        for (int j = 0; j < Ny; ++j)
            for (int m = 0; m < Ny; ++m)
                out[j * Ny + m] += G[j * kD + m] - (fit ? G[j * kD + Ny] * G[m * kD + Ny] / hh : 0.0);
    }
    return out;
}

double worstRelative(const std::array<double, Ny * Ny>& a, const std::array<double, Ny * Ny>& b)
{
    double worst = 0.0;
    for (int j = 0; j < Ny; ++j)
        for (int m = 0; m < Ny; ++m) {
            const double scale = std::sqrt(std::fabs(b[j * Ny + j] * b[m * Ny + m])) + 1e-300;
            worst = std::max(worst, std::fabs(a[j * Ny + m] - b[j * Ny + m]) / scale);
        }
    return worst;
}

bool checkAstroClosedForm()
{
    bool ok = true;
    const auto expectClose = [&](const char* what, const std::array<double, Ny * Ny>& got,
                                 const std::array<double, Ny * Ny>& want) {
        const double w = worstRelative(got, want);
        const bool good = std::isfinite(w) and w < 1e-12;
        std::cout << "#   closed form, " << what << ": worst relative difference " << w << (good ? "  -> ok\n" : "  -> FAIL\n");
        if (not good) { std::cerr << "FAIL [astro-closed-form] " << what << ": " << w << "\n"; ok = false; }
    };
    std::cout << "# astrometric closed form for several Rubin bands, on synthetic sums\n";

    // (i) Every band flat (h constant in each band, different between them): each band's offset takes out its own
    // mean, the reference is degenerate and nothing more is subtracted: the sum over bands of Swdd - Swd Swd^T / Sw.
    {
        std::array<AstSums, kRubinBands> bands{};
        bands[2] = syntheticBand(40, 0.1, [](int) { return 0.45; });
        bands[3] = syntheticBand(25, 0.9, [](int) { return 0.30; });
        bands[4] = syntheticBand(31, 1.7, [](int) { return 0.20; });
        std::array<double, Ny * Ny> want{};
        for (const int b : {2, 3, 4})
            for (int cc = 0; cc < 2; ++cc)
                for (int j = 0; j < Ny; ++j)
                    for (int m = 0; m < Ny; ++m)
                        want[j * Ny + m] += bands[b].Swdd[cc][j * kD + m] - bands[b].Swd[cc][j] * bands[b].Swd[cc][m] / bands[b].Sw;
        std::array<double, Ny * Ny> got{};
        fisher_detail::rubinBandInformation(bands, got);
        expectClose("every band flat: each band keeps its own reference", got, want);
    }
    // (ii) A band with no neighbours' light (h = 0, Swhh = 0) carries no offset and nothing divides by zero; the others
    // are as the brute-force solve with that band's offset column absent.
    {
        std::array<AstSums, kRubinBands> bands{};
        bands[2] = syntheticBand(40, 0.1, [](int m) { return 0.6 / (1.0 + 0.2 * m); });
        bands[3] = syntheticBand(25, 0.9, [](int) { return 0.0; });
        bands[4] = syntheticBand(31, 1.7, [](int m) { return 0.3 / (1.0 + 0.05 * m * m); });
        std::array<double, Ny * Ny> got{};
        fisher_detail::rubinBandInformation(bands, got);
        expectClose("a band with Swhh = 0 carries no offset", got, bruteForceRubin(bands));
    }
    // (iii) Several bands, none flat: the closed form is the brute-force solve over (physical, offsets, one reference).
    {
        std::array<AstSums, kRubinBands> bands{};
        bands[2] = syntheticBand(40, 0.1, [](int m) { return 0.6 / (1.0 + 0.2 * m); });
        bands[3] = syntheticBand(25, 0.9, [](int m) { return 0.5 / (1.0 + 0.1 * m); });
        bands[4] = syntheticBand(31, 1.7, [](int m) { return 0.3 / (1.0 + 0.05 * m * m); });
        bands[5] = syntheticBand(12, 2.3, [](int m) { return 0.2 + 0.02 * m; });
        std::array<double, Ny * Ny> got{};
        fisher_detail::rubinBandInformation(bands, got);
        expectClose("four bands, none flat, against the brute-force solve", got, bruteForceRubin(bands));
    }
    // (iv) One band: the closed form is today's Group/marginalise/physical route.
    {
        std::array<AstSums, kRubinBands> bands{};
        bands[2] = syntheticBand(40, 0.1, [](int m) { return 0.6 / (1.0 + 0.2 * m); });
        std::array<double, Ny * Ny> got{};
        fisher_detail::rubinBandInformation(bands, got);
        expectClose("one band, against today's route", got, todayRoute(bands[2]));
    }
    return ok;
}


// ---------------------------------------------------------------------------------------------
// The lens's own light fraction as a free parameter of the astrometric fit (--lens-light free).
// ---------------------------------------------------------------------------------------------

// One band's sums with the lens's-light derivative l (slot Ny + 1) as well: lOf(m) in x, a different one in y.
template <class H, class L> AstSums syntheticBandL(int n, double shift, H hOf, L lOf)
{
    AstSums k;
    for (int m = 0; m < n; ++m) {
        std::array<double, kD> dx{}, dy{};
        for (int j = 0; j < Ny; ++j) {
            dx[j] = (1.0 + 0.3 * j) * std::sin(0.37 * m + 1.1 * j + shift) + 0.2 * std::cos(0.11 * m * (j + 1));
            dy[j] = (1.0 - 0.1 * j) * std::cos(0.29 * m + 0.7 * j + 2.0 * shift) - 0.1 * std::sin(0.05 * m * (j + 2));
        }
        const double h = hOf(m), l = lOf(m);
        dx[Ny] = h;  dy[Ny] = 0.8 * h;
        dx[Ny + 1] = l;  dy[Ny + 1] = 0.6 * l + 0.05 * std::cos(0.2 * m + shift);
        k.add(1.0 + 0.3 * std::cos(0.9 * m + shift), dx, dy);
    }
    return k;
}

// Dense solve A X = B (n x n, n x nrhs), Gauss elimination with partial pivoting; X returned in B.
void solveDense(std::vector<double> A, std::vector<double>& B, int n, int nrhs)
{
    for (int p = 0; p < n; ++p) {
        int piv = p;
        for (int i = p + 1; i < n; ++i) if (std::fabs(A[i * n + p]) > std::fabs(A[piv * n + p])) piv = i;
        for (int j = 0; j < n; ++j) std::swap(A[p * n + j], A[piv * n + j]);
        for (int j = 0; j < nrhs; ++j) std::swap(B[p * nrhs + j], B[piv * nrhs + j]);
        for (int i = p + 1; i < n; ++i) {
            const double f = A[i * n + p] / A[p * n + p];
            for (int j = p; j < n; ++j) A[i * n + j] -= f * A[p * n + j];
            for (int j = 0; j < nrhs; ++j) B[i * nrhs + j] -= f * B[p * nrhs + j];
        }
    }
    for (int p = n - 1; p >= 0; --p)
        for (int j = 0; j < nrhs; ++j) {
            double v = B[p * nrhs + j];
            for (int m = p + 1; m < n; ++m) v -= A[p * n + m] * B[m * nrhs + j];
            B[p * nrhs + j] = v / A[p * n + p];
        }
}

// X (n x n symmetric, row-major) with the last `nn` unknowns marginalised: the Schur complement on the first n - nn.
std::vector<double> schurOut(const std::vector<double>& X, int n, int nn)
{
    const int nk = n - nn;
    std::vector<double> out(nk * nk);
    for (int i = 0; i < nk; ++i) for (int j = 0; j < nk; ++j) out[i * nk + j] = X[i * n + j];
    if (nn == 0) return out;
    std::vector<double> Ann(nn * nn), Ank(nn * nk);
    for (int i = 0; i < nn; ++i) {
        for (int j = 0; j < nn; ++j) Ann[i * nn + j] = X[(nk + i) * n + nk + j];
        for (int j = 0; j < nk; ++j) Ank[i * nk + j] = X[(nk + i) * n + j];
    }
    std::vector<double> sol = Ank;
    solveDense(Ann, sol, nn, nk);
    for (int i = 0; i < nk; ++i)
        for (int j = 0; j < nk; ++j)
            for (int p = 0; p < nn; ++p) out[i * nk + j] -= X[i * n + nk + p] * sol[p * nk + j];
    return out;
}

// The brute force for the lens's-light fit: per coordinate the full information on (the four physical parameters, the
// kept bands' fractions, one offset for every band whose offset derivative is not zero, ONE reference position), the
// offsets and the reference marginalised by a direct solve; the two coordinates added, the priors added to the
// fractions, and the fractions marginalised. flInfo < 0: the band keeps no fraction.
std::array<double, Ny * Ny> bruteForceLensFree(const std::array<AstSums, kRubinBands>& bands,
                                               const std::array<double, kRubinBands>& flInfo)
{
    std::array<int, kRubinBands> fSlot{};
    int k = 0;
    for (int b = 0; b < kRubinBands; ++b) fSlot[b] = (bands[b].Sw > 0.0 and flInfo[b] >= 0.0) ? k++ : -1;
    const int nk = Ny + k;
    std::vector<double> T(nk * nk, 0.0);
    for (int cc = 0; cc < 2; ++cc) {
        std::array<int, kRubinBands> oSlot{};
        int no = 0;
        for (int b = 0; b < kRubinBands; ++b)
            oSlot[b] = (bands[b].Sw > 0.0 and bands[b].Swdd[cc][Ny * kD + Ny] > 0.0) ? no++ : -1;
        const int n = nk + no + 1;
        std::vector<double> X(n * n, 0.0);
        for (int b = 0; b < kRubinBands; ++b) {
            const AstSums& S = bands[b];
            if (not (S.Sw > 0.0)) continue;
            // the band's components: slots of the full vector (-1: not a parameter) and the sums' own slots
            std::vector<int> full, mine;
            for (int j = 0; j < Ny; ++j) { full.push_back(j);  mine.push_back(j); }
            if (fSlot[b] >= 0) { full.push_back(Ny + fSlot[b]);  mine.push_back(Ny + 1); }
            if (oSlot[b] >= 0) { full.push_back(nk + oSlot[b]);  mine.push_back(Ny); }
            full.push_back(n - 1);  mine.push_back(-1);              // the reference: the constant 1
            for (size_t a = 0; a < full.size(); ++a)
                for (size_t c = 0; c < full.size(); ++c) {
                    const int ma = mine[a], mc = mine[c];
                    double v;
                    if (ma >= 0 and mc >= 0) v = S.Swdd[cc][ma * kD + mc];
                    else if (ma >= 0)        v = S.Swd[cc][ma];
                    else if (mc >= 0)        v = S.Swd[cc][mc];
                    else                     v = S.Sw;
                    X[full[a] * n + full[c]] += v;
                }
        }
        const std::vector<double> M = schurOut(X, n, no + 1);
        for (int i = 0; i < nk * nk; ++i) T[i] += M[i];
    }
    for (int b = 0; b < kRubinBands; ++b) if (fSlot[b] >= 0) T[(Ny + fSlot[b]) * nk + Ny + fSlot[b]] += flInfo[b];
    const std::vector<double> R = schurOut(T, nk, k);
    std::array<double, Ny * Ny> out{};
    for (int j = 0; j < Ny; ++j) for (int m = 0; m < Ny; ++m) out[j * Ny + m] = R[j * Ny + m];
    return out;
}

// l -> l + alpha x_src on the sums of one band: src = Ny (the offset's derivative h) or -1 (the constant 1, the reference).
void shiftLensColumn(AstSums& S, int src, double alpha)
{
    for (int cc = 0; cc < 2; ++cc) {
        const auto old = S.Swdd[cc];
        const auto oldd = S.Swd[cc];
        const auto cross = [&](int a) { return src >= 0 ? old[a * kD + src] : oldd[a]; };   // sum w x_a x_src
        const double oneSrc = src >= 0 ? oldd[src] : S.Sw;                                 // sum w x_src
        const double srcSrc = src >= 0 ? old[src * kD + src] : S.Sw;
        for (int a = 0; a < kD; ++a) if (a != Ny + 1) {
            S.Swdd[cc][a * kD + Ny + 1] = old[a * kD + Ny + 1] + alpha * cross(a);
            S.Swdd[cc][(Ny + 1) * kD + a] = S.Swdd[cc][a * kD + Ny + 1];
        }
        S.Swdd[cc][(Ny + 1) * kD + Ny + 1] = old[(Ny + 1) * kD + Ny + 1] + 2.0 * alpha * cross(Ny + 1) + alpha * alpha * srcSrc;
        S.Swd[cc][Ny + 1] = oldd[Ny + 1] + alpha * oneSrc;
    }
}

bool checkLensLightFreeFunction()
{
    bool ok = true;
    const auto expectClose = [&](const char* what, const std::array<double, Ny * Ny>& got,
                                 const std::array<double, Ny * Ny>& want) {
        const double w = worstRelative(got, want);
        const bool good = std::isfinite(w) and w < 1e-12;
        std::cout << "#   lens-light function, " << what << ": worst relative difference " << w << (good ? "  -> ok\n" : "  -> FAIL\n");
        if (not good) { std::cerr << "FAIL [lens-light-function] " << what << ": " << w << "\n"; ok = false; }
    };
    std::cout << "# the lens's-light reduction on synthetic sums\n";
    const auto lOf = [](double s0) { return [s0](int m) { return 0.4 * std::sin(0.13 * m + s0) + 0.3; }; };
    std::array<AstSums, kRubinBands> four{};
    four[2] = syntheticBandL(40, 0.1, [](int m) { return 0.6 / (1.0 + 0.2 * m); }, lOf(0.3));
    four[3] = syntheticBandL(25, 0.9, [](int m) { return 0.5 / (1.0 + 0.1 * m); }, lOf(1.1));
    four[4] = syntheticBandL(31, 1.7, [](int m) { return 0.3 / (1.0 + 0.05 * m * m); }, lOf(2.5));
    four[5] = syntheticBandL(12, 2.3, [](int m) { return 0.2 + 0.02 * m; }, lOf(0.8));
    const std::array<double, kRubinBands> info = {-1.0, -1.0, 12.0, 50.0, 3.0, 25.0};
    {
        std::array<double, Ny * Ny> got{};
        const bool good = fisher_detail::rubinBandInformationLensFree(four, info, got);
        if (not good) { std::cerr << "FAIL [lens-light-function] four bands did not invert\n"; ok = false; }
        expectClose("four bands against the brute-force solve", got, bruteForceLensFree(four, info));
    }
    {   // a band without a fraction parameter (negative precision) keeps its physical part, drops its l column
        std::array<double, kRubinBands> info2 = info;  info2[3] = -1.0;
        std::array<double, Ny * Ny> got{};
        fisher_detail::rubinBandInformationLensFree(four, info2, got);
        expectClose("a band with no fraction parameter", got, bruteForceLensFree(four, info2));
    }
    {   // one band: the same as the brute-force solve over (physical, fraction, offset, reference)
        std::array<AstSums, kRubinBands> one{};
        one[2] = four[2];
        std::array<double, Ny * Ny> got{};
        fisher_detail::rubinBandInformationLensFree(one, info, got);
        expectClose("one band", got, bruteForceLensFree(one, info));
    }
    // Adding alpha h (the offset's derivative) to every l of a band changes nothing, whatever the other bands do: the band's own
    // offset takes that direction out. Adding alpha 1 (the reference) does so for a single band, whose epochs the reference's
    // column covers whole; with several bands the constant on one band's epochs is not the shared reference's direction.
    {
        std::array<AstSums, kRubinBands> moved = four;
        shiftLensColumn(moved[3], Ny, 0.7);
        shiftLensColumn(moved[5], Ny, -1.9);
        std::array<double, Ny * Ny> a{}, b{};
        fisher_detail::rubinBandInformationLensFree(four, info, a);
        fisher_detail::rubinBandInformationLensFree(moved, info, b);
        expectClose("l + alpha h in two bands leaves FL unchanged", b, a);
    }
    {
        std::array<AstSums, kRubinBands> one{}, moved{};
        one[3] = four[3];  moved[3] = four[3];
        shiftLensColumn(moved[3], -1, -1.9);
        std::array<double, Ny * Ny> a{}, b{};
        fisher_detail::rubinBandInformationLensFree(one, info, a);
        fisher_detail::rubinBandInformationLensFree(moved, info, b);
        expectClose("l + alpha 1 in a single band leaves FL unchanged", b, a);
    }
    return ok;
}

// ---- the events ----
// The cases of the lens's-light checks: the kEvents events (every Rubin visit in r, dark lens), then the chromatic event
// (r, i, z at f_S 0.3 / 0.5 / 0.6) with a dark lens and with a luminous one (f_L 0.10 / 0.05 / 0.03 in r / i / z, 0.05 in F146).
constexpr int kLensCases = int(sizeof(kEvents) / sizeof(kEvents[0])) + 2;

const char* lensCaseName(int c)
{
    const int n = kLensCases - 2;
    return c < n ? kEvents[c].name : (c == n ? "chromatic_dark" : "chromatic_luminous");
}

std::unique_ptr<EventRun> makeLensCase(int c, bool lensFree, double priorFrac)
{
    const int n = kLensCases - 2;
    auto r = std::make_unique<EventRun>(c < n ? kEvents[c] : kEvents[2]);
    if (c >= n) makeChromatic(*r->s, *r->l, r->ndw);
    if (c == n + 1) {
        r->s->fLensF[2] = 0.10;  r->s->fLensF[3] = 0.05;  r->s->fLensF[4] = 0.03;
        r->s->fLens[0] = r->s->fLensF[2];
        r->s->fLensF[6] = 0.05;  r->s->fLens[1] = r->s->fLensF[6];
    }
    r->co->lensLightFree = lensFree;
    r->co->lensLightPriorFrac = priorFrac;
    r->run();
    return r;
}

constexpr double kUniformFrac = kLensLightUniformFrac;

// The chromatic events' astrometric forecast of the 5dd83b6 build, hexadecimal floats: dark and luminous, point source
// (espl 0) and disc (1). Generated once from that build (a temporary mode of this file); never regenerated to make a check pass.
const AstroRef kAstroRefChrom[] = {
    {0, "chrom_dark", 0, 1, 0x1.03e3de41f03b8p+2, {0x1.a5b9c5ee23f66p-7, 0x1.0088bdf6c4ef5p-15, 0x1.fc3ecbd9b9d87p-16, 0x1.3fa24f7f6a12fp-1}, {{1, 0x1.03e3de41f03b8p+2, {0x1.a5b9c5ee23f66p-7, 0x1.0088bdf6c4ef5p-15, 0x1.fc3ecbd9b9d87p-16, 0x1.3fa24f7f6a12fp-1}}, {1, 0x1.1ea99a669a882p+1, {0x1.3d916e4baa649p-3, 0x1.d7fa7c31ec548p-14, 0x1.dbd5016a0fb6ap-14, 0x1.24a1356d2a157p+2}}, {1, 0x1.1b4161c49dee6p+1, {0x1.4d2da1ad0b179p-2, 0x1.e3919a83528dap-14, 0x1.e51d88ea13889p-14, 0x1.49a069c182f53p+2}}}},
    {0, "chrom_dark", 1, 1, 0x1.459677648262p+1, {0x1.b9b2b839278cap-2, 0x1.ee47205023d29p-14, 0x1.ec5b9f631e031p-14, 0x1.53eab663a7f79p+2}, {{1, 0x1.459677648262p+1, {0x1.b9b2b839278cap-2, 0x1.ee47205023d29p-14, 0x1.ec5b9f631e031p-14, 0x1.53eab663a7f79p+2}}, {1, 0x1.459677648262p+1, {0x1.b9b2b839278cap-2, 0x1.ee47205023d29p-14, 0x1.ec5b9f631e031p-14, 0x1.53eab663a7f79p+2}}, {1, 0x1.459677648262p+1, {0x1.b9b2b839278cap-2, 0x1.ee47205023d29p-14, 0x1.ec5b9f631e031p-14, 0x1.53eab663a7f79p+2}}}},
    {0, "chrom_dark", 2, 1, 0x1.07f19b25add7fp+2, {0x1.a7c6b8f90b74bp-7, 0x1.0ad61752e8791p-15, 0x1.07d28c97a79f1p-15, 0x1.4298c985131ccp-1}, {{1, 0x1.07f19b25add7fp+2, {0x1.a7c6b8f90b74bp-7, 0x1.0ad61752e8791p-15, 0x1.07d28c97a79f1p-15, 0x1.4298c985131ccp-1}}, {1, 0x1.3d5c54a4df2afp+3, {0x1.00a3c94d685ccp-2, 0x1.c3d82d354e6e4p-10, 0x1.78fd3a912d317p-10, 0x1.7a758c5975e23p+3}}, {1, 0x1.755d2d987bc7ep+2, {0x1.503aeed79f6b8p-1, 0x1.25ba7467df75p-8, 0x1.107f7c13cb435p-8, 0x1.867dd35b98cddp+4}}}},
    {0, "chrom_lum", 0, 1, 0x1.9f1ca21215c9p+3, {0x1.6bf158522d989p-7, 0x1.e28f8a247997ap-16, 0x1.c1d31ef306f98p-15, 0x1.51d7543005ab5p-2}, {{1, 0x1.9f1ca21215c9p+3, {0x1.6bf158522d989p-7, 0x1.e28f8a247997ap-16, 0x1.c1d31ef306f98p-15, 0x1.51d7543005ab5p-2}}, {1, 0x1.9cc37f1a5c349p+3, {0x1.39580627cbc6fp-4, 0x1.a925e168af2d1p-13, 0x1.80643d72f0a27p-13, 0x1.43ec726665c07p+0}}, {1, 0x1.0baa6222c3a6cp+4, {0x1.684a6a3f5ac92p-4, 0x1.c25fe82b16ae8p-13, 0x1.a30d02163506ep-13, 0x1.48beafa54456bp+0}}}},
    {0, "chrom_lum", 1, 1, 0x1.15d7dad5578e2p+4, {0x1.6f7c71fd511cbp-4, 0x1.c663e7a9ba8c7p-13, 0x1.a87b248b6a7c9p-13, 0x1.499d6aaf25798p+0}, {{1, 0x1.15d7dad5578e2p+4, {0x1.6f7c71fd511cbp-4, 0x1.c663e7a9ba8c7p-13, 0x1.a87b248b6a7c9p-13, 0x1.499d6aaf25798p+0}}, {1, 0x1.15d7dad5578e2p+4, {0x1.6f7c71fd511cbp-4, 0x1.c663e7a9ba8c7p-13, 0x1.a87b248b6a7c9p-13, 0x1.499d6aaf25798p+0}}, {1, 0x1.15d7dad5578e2p+4, {0x1.6f7c71fd511cbp-4, 0x1.c663e7a9ba8c7p-13, 0x1.a87b248b6a7c9p-13, 0x1.499d6aaf25798p+0}}}},
    {0, "chrom_lum", 2, 1, 0x1.dc774e6972ad5p+4, {0x1.738d55282a385p-7, 0x1.f72b03f58185p-16, 0x1.51060d6148b5bp-14, 0x1.fc370e2349e6p-2}, {{1, 0x1.dc774e6972ad5p+4, {0x1.738d55282a385p-7, 0x1.f72b03f58185p-16, 0x1.51060d6148b5bp-14, 0x1.fc370e2349e6p-2}}, {1, 0x1.896f937de8bccp+2, {0x1.ed937d88b34b5p-3, 0x1.5b070c7aa189p-10, 0x1.48aaeb4fda2c4p-10, 0x1.08ca22b789c41p+4}}, {1, 0x1.93bfdfcaf6268p+1, {0x1.3a7c733c03cedp-1, 0x1.f8fad008de82fp-9, 0x1.e742a87a581d5p-9, 0x1.f0b8685f50edep+4}}}},
    {1, "chrom_dark", 0, 1, 0x1.03e3d976e4b2ep+2, {0x1.a5b9e761f5e2cp-7, 0x1.0088c0f85e8a7p-15, 0x1.fc3ec271dad86p-16, 0x1.3fa255dcdab6bp-1}, {{1, 0x1.03e3d976e4b2ep+2, {0x1.a5b9e761f5e2cp-7, 0x1.0088c0f85e8a7p-15, 0x1.fc3ec271dad86p-16, 0x1.3fa255dcdab6bp-1}}, {1, 0x1.1ea993cd20407p+1, {0x1.3d91873e76101p-3, 0x1.d7fa7d14c15e7p-14, 0x1.dbd50218d15d8p-14, 0x1.24a139fb27902p+2}}, {1, 0x1.1b41693b2f0ebp+1, {0x1.4d2db75e0362cp-2, 0x1.e3919c8f8fa3dp-14, 0x1.e51d8a5c1acp-14, 0x1.49a06d4eaf0cep+2}}}},
    {1, "chrom_dark", 1, 1, 0x1.45967a16f8fbcp+1, {0x1.b9b2c6ae07484p-2, 0x1.ee472196373aap-14, 0x1.ec5ba126b1d31p-14, 0x1.53eab7a3c3f52p+2}, {{1, 0x1.45967a16f8fbcp+1, {0x1.b9b2c6ae07484p-2, 0x1.ee472196373aap-14, 0x1.ec5ba126b1d31p-14, 0x1.53eab7a3c3f52p+2}}, {1, 0x1.45967a16f8fbcp+1, {0x1.b9b2c6ae07484p-2, 0x1.ee472196373aap-14, 0x1.ec5ba126b1d31p-14, 0x1.53eab7a3c3f52p+2}}, {1, 0x1.45967a16f8fbcp+1, {0x1.b9b2c6ae07484p-2, 0x1.ee472196373aap-14, 0x1.ec5ba126b1d31p-14, 0x1.53eab7a3c3f52p+2}}}},
    {1, "chrom_dark", 2, 1, 0x1.07f1968a2638cp+2, {0x1.a7c6dae48cf8cp-7, 0x1.0ad61aaf597b3p-15, 0x1.07d28758cd3cap-15, 0x1.4298d0063177fp-1}, {{1, 0x1.07f1968a2638cp+2, {0x1.a7c6dae48cf8cp-7, 0x1.0ad61aaf597b3p-15, 0x1.07d28758cd3cap-15, 0x1.4298d0063177fp-1}}, {1, 0x1.3d5c7518a5bcp+3, {0x1.00a3e9839adcdp-2, 0x1.c3d846bedb409p-10, 0x1.78fd37f41c876p-10, 0x1.7a75ad693a0c3p+3}}, {1, 0x1.755da20c03405p+2, {0x1.503b304a2b0f2p-1, 0x1.25ba8abe4ffb8p-8, 0x1.107f7acfd9b53p-8, 0x1.867e119aecf8p+4}}}},
    {1, "chrom_lum", 0, 1, 0x1.9f1cb550abe5ep+3, {0x1.6bf176a0c1e57p-7, 0x1.e28f85c7ccc7ap-16, 0x1.c1d31d5b18b6bp-15, 0x1.51d755bdb38bap-2}, {{1, 0x1.9f1cb550abe5ep+3, {0x1.6bf176a0c1e57p-7, 0x1.e28f85c7ccc7ap-16, 0x1.c1d31d5b18b6bp-15, 0x1.51d755bdb38bap-2}}, {1, 0x1.9cc390720b9d7p+3, {0x1.39580dd993cdbp-4, 0x1.a925e45819ce1p-13, 0x1.806442e893878p-13, 0x1.43ec73632f86ap+0}}, {1, 0x1.0baa64f4f226p+4, {0x1.684a6c4b35092p-4, 0x1.c25fe8fec40fp-13, 0x1.a30d03be25d14p-13, 0x1.48beb01482c73p+0}}}},
    {1, "chrom_lum", 1, 1, 0x1.15d7dbfc0a62cp+4, {0x1.6f7c72bf71819p-4, 0x1.c663e7bf5f546p-13, 0x1.a87b256808425p-13, 0x1.499d6aff4b5dap+0}, {{1, 0x1.15d7dbfc0a62cp+4, {0x1.6f7c72bf71819p-4, 0x1.c663e7bf5f546p-13, 0x1.a87b256808425p-13, 0x1.499d6aff4b5dap+0}}, {1, 0x1.15d7dbfc0a62cp+4, {0x1.6f7c72bf71819p-4, 0x1.c663e7bf5f546p-13, 0x1.a87b256808425p-13, 0x1.499d6aff4b5dap+0}}, {1, 0x1.15d7dbfc0a62cp+4, {0x1.6f7c72bf71819p-4, 0x1.c663e7bf5f546p-13, 0x1.a87b256808425p-13, 0x1.499d6aff4b5dap+0}}}},
    {1, "chrom_lum", 2, 1, 0x1.dc7763796d189p+4, {0x1.738d744698055p-7, 0x1.f72afe96a8efep-16, 0x1.51060b57c1074p-14, 0x1.fc3710c019f55p-2}, {{1, 0x1.dc7763796d189p+4, {0x1.738d744698055p-7, 0x1.f72afe96a8efep-16, 0x1.51060b57c1074p-14, 0x1.fc3710c019f55p-2}}, {1, 0x1.896faf4abf66bp+2, {0x1.ed93bcd9863a4p-3, 0x1.5b070aa780bfap-10, 0x1.48aae82908cbep-10, 0x1.08ca49c1f67c2p+4}}, {1, 0x1.93c04ef2c07c4p+1, {0x1.3a7cab348be45p-1, 0x1.f8face549e2bep-9, 0x1.e742a5a42cddap-9, 0x1.f0b8cc4e9f1c9p+4}}}},
};

bool matchesRef(const covarian& c, const AstroRef& ref)
{
    const int q = ref.q;
    bool same = c.okB[q] == ref.ok and sameDouble(c.condB[q], ref.cond);
    for (int k = 0; k < Ny; ++k) same = same and sameDouble(c.Erb[q][k], ref.erb[k]);
    for (int v = 0; v < NAVAR; ++v) {
        same = same and c.okBV[v][q] == ref.v[v].ok and sameDouble(c.condBV[v][q], ref.v[v].cond);
        for (int k = 0; k < Ny; ++k) same = same and sameDouble(c.ErbV[v][q][k], ref.v[v].erb[k]);
    }
    return same;
}

// In known mode (the covarian default; the run default is free) every astrometric output is the 5dd83b6 one, bit for bit: the kEvents events against the
// table of the Task 3 check and the chromatic events against theirs. A free-mode run before each known one checks that
// nothing the free mode touches leaks into the next.
bool checkLensLightKnownUnchanged()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    int n = 0;
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        const int espl = model == SourceModel::ESPL ? 1 : 0;
        for (int c = 0; c < kLensCases; ++c) {
            makeLensCase(c, true, kUniformFrac);                       // a free-mode run first
            const auto r = makeLensCase(c, false, kUniformFrac);
            const bool chrom = c >= kLensCases - 2;
            const char* refName = chrom ? (c == kLensCases - 2 ? "chrom_dark" : "chrom_lum") : kEvents[c].name;
            const AstroRef* table = chrom ? kAstroRefChrom : kAstroRef;
            const size_t tn = chrom ? sizeof(kAstroRefChrom) / sizeof(kAstroRefChrom[0]) : sizeof(kAstroRef) / sizeof(kAstroRef[0]);
            for (size_t i = 0; i < tn; ++i) {
                const AstroRef& ref = table[i];
                if (ref.espl != espl or std::strcmp(ref.name, refName) != 0) continue;
                ++n;
                if (not matchesRef(*r->co, ref)) {
                    std::cerr << "FAIL [lens-light-known] " << lensCaseName(c) << (espl ? " (espl)" : " (pspl)") << " partition "
                              << ref.q << " is not bit-identical to the 5dd83b6 astrometric forecast\n";
                    ok = false;
                }
            }
        }
    }
    vbmSetSourceModel(keep);
    if (n != 2 * NSURV * kLensCases) { std::cerr << "FAIL [lens-light-known] compared " << n << " rows\n"; ok = false; }
    std::cout << "# --lens-light known: every astrometric output equals the 5dd83b6 one bit for bit: " << n << " rows, "
              << (ok ? "ok\n" : "FAIL\n");
    return ok;
}

// Every astrometric sigma of one run, in a list: Erb of each partition, then the three variants' ErbV; -1 where not ok.
std::vector<double> astroSigmas(const covarian& c)
{
    std::vector<double> out;
    for (int q = 0; q < NSURV; ++q) {
        for (int k = 0; k < Ny; ++k) out.push_back(c.Erb[q][k]);
        for (int v = 0; v < NAVAR; ++v) for (int k = 0; k < Ny; ++k) out.push_back(c.ErbV[v][q][k]);
    }
    return out;
}

// A free fraction is one more nuisance parameter: no sigma can fall below its known-mode value, for either prior.
bool checkLensLightFreeNoGain()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        for (int c = 0; c < kLensCases; ++c) {
            const auto known = makeLensCase(c, false, kUniformFrac);
            const std::vector<double> sk = astroSigmas(*known->co);
            for (const double frac : {kUniformFrac, 1.0}) {
                const auto fr = makeLensCase(c, true, frac);
                const std::vector<double> sf = astroSigmas(*fr->co);
                double worst = 1e300, thetaRatio = -1.0;
                for (size_t i = 0; i < sk.size(); ++i) {
                    if (sk[i] < 0.0) continue;                                    // not measured in known mode
                    if (not (sf[i] >= 0.0)) {
                        std::cerr << "FAIL [lens-light-nogain] " << lensCaseName(c) << ": a sigma measured in known mode is not in free mode\n";
                        ok = false;  continue;
                    }
                    worst = std::min(worst, sf[i] / sk[i] - 1.0);
                }
                if (known->co->okB[SJOINT]) thetaRatio = fr->co->Erb[SJOINT][0] / known->co->Erb[SJOINT][0];
                if (not (worst >= -1e-12)) {
                    std::cerr << "FAIL [lens-light-nogain] " << lensCaseName(c) << (model == SourceModel::ESPL ? " (espl)" : " (pspl)")
                              << ", prior frac " << frac << ": a free sigma is below the known one by " << -worst << " relative\n";
                    ok = false;
                }
                // Not vacuous: where Rubin and Roman both watch the event, the free fraction costs thetaE a visible amount.
                if (c >= kLensCases - 2 and not (thetaRatio > 1.05)) {
                    std::cerr << "FAIL [lens-light-nogain] " << lensCaseName(c) << ": free mode left sigma(thetaE) at " << thetaRatio
                              << " x known (expected > 1.05)\n";
                    ok = false;
                }
                std::cout << "#   " << lensCaseName(c) << (model == SourceModel::ESPL ? " (espl)" : " (pspl)") << ", prior x" << frac
                          << ": smallest free/known - 1 = " << worst << ", joint sigma(thetaE) free/known = " << thetaRatio << "\n";
            }
        }
    }
    vbmSetSourceModel(keep);
    std::cout << "# free lens light never lowers an astrometric sigma" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// The tighter the prior on the fraction, the smaller every sigma: uniform <= blend <= 100 x blend; a prior of ~0 width
// is the known mode.
bool checkLensLightPriorMonotone()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        for (int c = 0; c < kLensCases; ++c) {
            const auto known = makeLensCase(c, false, kUniformFrac);
            const auto tiny  = makeLensCase(c, true, 1e-9);
            const auto uni   = makeLensCase(c, true, kUniformFrac);
            const auto blend = makeLensCase(c, true, 1.0);
            const auto wide  = makeLensCase(c, true, 100.0);
            const auto sk = astroSigmas(*known->co), st = astroSigmas(*tiny->co), su = astroSigmas(*uni->co),
                       sb = astroSigmas(*blend->co), sw = astroSigmas(*wide->co);
            double worstOrder = 1e300, worstTiny = 0.0;
            for (size_t i = 0; i < sk.size(); ++i) {
                if (sk[i] < 0.0) continue;
                worstOrder = std::min(worstOrder, std::min(sb[i] / su[i] - 1.0, sw[i] / sb[i] - 1.0));
                worstTiny = std::max(worstTiny, std::fabs(st[i] / sk[i] - 1.0));
            }
            if (not (worstOrder >= -1e-12)) {
                std::cerr << "FAIL [lens-light-prior] " << lensCaseName(c) << (model == SourceModel::ESPL ? " (espl)" : " (pspl)")
                          << ": a sigma is not monotone in the prior, by " << -worstOrder << " relative\n";
                ok = false;
            }
            if (not (worstTiny < 1e-8)) {
                std::cerr << "FAIL [lens-light-prior] " << lensCaseName(c) << (model == SourceModel::ESPL ? " (espl)" : " (pspl)")
                          << ": a prior of 1e-9 differs from known mode by " << worstTiny << " relative\n";
                ok = false;
            }
            std::cout << "#   " << lensCaseName(c) << (model == SourceModel::ESPL ? " (espl)" : " (pspl)")
                      << ": smallest ordering margin " << worstOrder << ", prior 1e-9 vs known " << worstTiny << "\n";
        }
    }
    vbmSetSourceModel(keep);
    std::cout << "# the prior on the lens's light: uniform <= blend <= 100 x blend, and ~0 reproduces known mode"
              << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// The analytic derivative of the centroid with respect to the lens's light fraction against a central difference in it:
// what is differenced is the centroid minus the neighbours' share (their share, fn blendOff / den, also moves with f_L and
// belongs to the offset's derivative). Rubin's three bands and Roman's F146, with an offset of the neighbours of (3, -2) mas.
bool checkLensLightDerivative()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        EventRun r(kEvents[2]);
        makeChromatic(*r.s, *r.l, r.ndw);
        source& s = *r.s;
        s.fLensF[2] = 0.10;  s.fLensF[3] = 0.05;  s.fLensF[4] = 0.03;  s.fLens[0] = s.fLensF[2];
        s.fLensF[6] = 0.05;  s.fLens[1] = s.fLensF[6];
        s.blendOff[0] = {3.0, -2.0};  s.blendOff[1] = {3.0, -2.0};
        for (int tele = 0; tele < 2; ++tele) {
            std::vector<double> t;  std::vector<int> filt;
            for (int i = 0; i < r.ndw; ++i) if (r.l->tele[i] == tele) { t.push_back(r.l->timn[i]);  filt.push_back(r.l->filt[i]); }
            const int n = int(t.size());
            const int* fp = tele == 0 ? filt.data() : nullptr;
            // the bands tested: Rubin's r, i, z (the epochs with filt == b) and F146 (all of Roman's)
            for (const int b : (tele == 0 ? std::vector<int>{2, 3, 4} : std::vector<int>{6})) {
                std::vector<ModelPoint> ref(n), up(n), dn(n);
                evaluateModel(s, *r.l, *r.as, tele, t.data(), n, true, ref.data(), fp, true);
                const double d = 1e-3;
                const auto setFL = [&](double v) { s.fLensF[b] = v;  if (b == 2) s.fLens[0] = v;  if (b == 6) s.fLens[1] = v; };
                const double keepFL = s.fLensF[b];
                setFL(keepFL + d);  evaluateModel(s, *r.l, *r.as, tele, t.data(), n, true, up.data(), fp);
                setFL(keepFL - d);  evaluateModel(s, *r.l, *r.as, tele, t.data(), n, true, dn.data(), fp);
                setFL(keepFL);
                double worst = 0.0, scale = 0.0;
                int used = 0;
                for (int c = 0; c < 2; ++c) {
                    std::vector<double> fd, an;
                    for (int m = 0; m < n; ++m) {
                        const int fb = fp ? fp[m] : 6;
                        if (fb != b) continue;
                        const double fs = fp ? s.blend[fb] : s.fb[1];
                        const double den = fs * ref[m].A + 1.0 - fs;
                        const auto share = [&](double fL) { return std::max(0.0, 1.0 - fs - fL) * s.blendOff[tele][c] / den; };
                        const double pu = (c ? up[m].pos2c : up[m].pos1c) - share(keepFL + d);
                        const double pd = (c ? dn[m].pos2c : dn[m].pos1c) - share(keepFL - d);
                        fd.push_back((pu - pd) / (2.0 * d));
                        an.push_back(c ? ref[m].lensDy : ref[m].lensDx);
                    }
                    used = int(an.size());
                    for (double a : an) scale = std::max(scale, std::fabs(a));
                    for (size_t m = 0; m < an.size(); ++m) worst = std::max(worst, std::fabs(fd[m] - an[m]));
                }
                const bool good = used > 0 and scale > 0.0 and worst <= 1e-6 * scale;
                if (not good) {
                    std::cerr << "FAIL [lens-light-derivative] " << (tele ? "Roman F146" : "Rubin band ") << (tele ? 0 : b)
                              << ": analytic l differs from the finite difference by " << worst << " (scale " << scale << ", " << used << " epochs)\n";
                    ok = false;
                }
                std::cout << "#   " << (model == SourceModel::ESPL ? "espl " : "pspl ") << (tele ? "F146" : "Rubin band ")
                          << (tele ? 0 : b) << ": " << used << " epochs, largest |l| " << scale << " mas, worst |analytic - difference| "
                          << worst << (good ? "  -> ok\n" : "  -> FAIL\n");
            }
        }
    }
    vbmSetSourceModel(keep);
    std::cout << "# d(centroid)/d(f_L) against a central difference" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// ---- Roman's reduction with the lens's light fraction free, on synthetic day sums ----

// Synthetic Roman epochs: ndays days of perDay epochs, the days' rolls alternating; per epoch a magnification track, the
// neighbours' share h (same in x and y), random physical derivatives and the lens-light derivative l = -0.6 u / den
// (u = 0.3 + t in x, -0.2 t in y). The lists keep what the brute force needs.
struct RomanEpoch { int day; double w, h; std::array<double, kD> dx, dy; };

std::vector<RomanEpoch> syntheticRoman(int ndays, int perDay)
{
    std::mt19937 gen(11);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::normal_distribution<double> nrm(0.0, 1.0);
    std::vector<RomanEpoch> out;
    for (int day = 0; day < ndays; ++day)
        for (int e = 0; e < perDay; ++e) {
            const double t = -3.0 + 6.0 * uni(gen);
            const double A = 1.0 + 2.0 / (1.0 + t * t), fs = 0.4, fL = 0.2, den = fs * A + 1.0 - fs;
            RomanEpoch ep{day, 0.5 + 1.5 * uni(gen), (1.0 - fs - fL) / den, {}, {}};
            for (int c = 0; c < 2; ++c) {
                auto& d = c ? ep.dy : ep.dx;
                for (int j = 0; j < Ny; ++j) d[j] = nrm(gen);
                d[1 + c] = t;
                d[Ny] = ep.h;
                d[Ny + 1] = -0.6 * (c ? -0.2 * t : 0.3 + t) / den;
            }
            out.push_back(ep);
        }
    return out;
}

// The brute force: parameters (physical, f_L if withL, then per coordinate one reference per roll group, one offset per
// day with prior precision 1/s2 (the frame errors), and the neighbours' offset), all but the physical ones marginalised
// by a direct solve; f_L carries the prior precision flInfo.
std::array<double, Ny * Ny> bruteRoman(const std::vector<RomanEpoch>& eps, int ndays, int nroll, double s2, bool withL, double flInfo)
{
    const int per = nroll + ndays + 1, nl = withL ? 1 : 0, n = Ny + nl + 2 * per;
    std::vector<double> F(n * n, 0.0);
    for (const auto& ep : eps)
        for (int c = 0; c < 2; ++c) {
            std::vector<double> g(n, 0.0);
            const auto& d = c ? ep.dy : ep.dx;
            for (int j = 0; j < Ny; ++j) g[j] = d[j];
            if (withL) g[Ny] = d[Ny + 1];
            const int base = Ny + nl + c * per;
            g[base + (nroll == 1 ? 0 : ep.day % 2)] = 1.0;
            g[base + nroll + ep.day] = 1.0;
            g[base + nroll + ndays] = d[Ny];
            for (int a = 0; a < n; ++a) for (int b = 0; b < n; ++b) F[a * n + b] += ep.w * g[a] * g[b];
        }
    if (s2 > 0.0)
        for (int c = 0; c < 2; ++c)
            for (int day = 0; day < ndays; ++day) { const int i = Ny + nl + c * per + nroll + day;  F[i * n + i] += 1.0 / s2; }
    if (withL) F[Ny * n + Ny] += flInfo;
    const std::vector<double> R = schurOut(F, n, n - Ny);
    std::array<double, Ny * Ny> out{};
    for (int i = 0; i < Ny * Ny; ++i) out[i] = R[i];
    return out;
}

bool checkLensLightRoman()
{
    bool ok = true;
    const int nd = 12;
    const auto eps = syntheticRoman(nd, 8);
    std::map<long, AstSums> days;
    for (const auto& ep : eps) {
        AstSums& k = days[ep.day];
        if (k.season < 0) { k.season = 0;  k.roll = ep.day % 2; }
        k.add(ep.w, ep.dx, ep.dy);
    }
    std::cout << "# Roman's reduction with f_L free, on synthetic day sums (day frames at sigma_c, a free reference per roll, the offset, f_L shared by x and y)\n";
    const auto expectClose = [&](const char* what, const std::array<double, Ny * Ny>& got, const std::array<double, Ny * Ny>& want) {
        const double w = worstRelative(got, want);
        const bool good = std::isfinite(w) and w < 1e-12;
        std::cout << "#   " << what << ": worst relative difference " << w << (good ? "  -> ok\n" : "  -> FAIL\n");
        if (not good) { std::cerr << "FAIL [lens-light-roman] " << what << ": " << w << "\n"; ok = false; }
    };
    for (const double s2 : {0.09, 1.0})
        for (const double prior : {0.6, 0.6 * kLensLightUniformFrac})
            for (const bool oneFrame : {false, true}) {
                std::array<double, Ny * Ny> got{};
                fisher_detail::romanInformation(days, oneFrame, std::sqrt(s2), 1.0 / (prior * prior), got);
                char buf[128];
                std::snprintf(buf, sizeof buf, "sigma_c^2 %.2g, prior sigma %.3f, %s", s2, prior, oneFrame ? "one frame" : "frame per roll");
                expectClose(buf, got, bruteRoman(eps, nd, oneFrame ? 1 : 2, s2, true, 1.0 / (prior * prior)));
            }
    {   // no f_L parameter (flInfo < 0): today's reduction
        std::array<double, Ny * Ny> got{};
        fisher_detail::romanInformation(days, false, 0.3, -1.0, got);
        expectClose("no f_L parameter, sigma_c 0.3", got, bruteRoman(eps, nd, 2, 0.09, false, 0.0));
    }
    {   // sigma_c = 0, one frame: the one-band closed form (offset, one reference, f_L) on the same epochs
        AstSums all;
        for (const auto& ep : eps) all.add(ep.w, ep.dx, ep.dy);
        std::array<AstSums, kRubinBands> one{};
        one[2] = all;
        std::array<double, kRubinBands> info = {-1.0, -1.0, 25.0, -1.0, -1.0, -1.0};
        std::array<double, Ny * Ny> got{};
        fisher_detail::romanInformation(days, true, 0.0, 25.0, got);
        expectClose("sigma_c = 0, one frame, against bruteForceLensFree on the same epochs", got, bruteForceLensFree(one, info));
    }
    return ok;
}

// ---- the branches of rubinBandInformationLensFree the four-band test does not reach ----
bool checkLensLightFreeBranches()
{
    bool ok = true;
    const auto expectClose = [&](const char* what, const std::array<double, Ny * Ny>& got, const std::array<double, Ny * Ny>& want) {
        const double w = worstRelative(got, want);
        const bool good = std::isfinite(w) and w < 1e-12;
        std::cout << "#   lens-light function, " << what << ": worst relative difference " << w << (good ? "  -> ok\n" : "  -> FAIL\n");
        if (not good) { std::cerr << "FAIL [lens-light-branches] " << what << ": " << w << "\n"; ok = false; }
    };
    const auto lOf = [](double s0) { return [s0](int m) { return 0.4 * std::sin(0.13 * m + s0) + 0.3; }; };
    {   // a band with no neighbours' light (h = 0, Swhh = 0): no offset, the others as the brute force
        std::array<AstSums, kRubinBands> bands{};
        bands[2] = syntheticBandL(40, 0.1, [](int m) { return 0.6 / (1.0 + 0.2 * m); }, lOf(0.3));
        bands[3] = syntheticBandL(25, 0.9, [](int) { return 0.0; }, lOf(1.1));
        bands[4] = syntheticBandL(31, 1.7, [](int m) { return 0.3 / (1.0 + 0.05 * m * m); }, lOf(2.5));
        const std::array<double, kRubinBands> info = {-1.0, -1.0, 12.0, 50.0, 3.0, -1.0};
        std::array<double, Ny * Ny> got{};
        fisher_detail::rubinBandInformationLensFree(bands, info, got);
        expectClose("a band with h = 0", got, bruteForceLensFree(bands, info));
    }
    {   // h constant in every band: each band's offset takes its own mean out, the reference is degenerate and not subtracted
        std::array<AstSums, kRubinBands> bands{};
        bands[2] = syntheticBandL(40, 0.1, [](int) { return 0.45; }, lOf(0.3));
        bands[3] = syntheticBandL(25, 0.9, [](int) { return 0.30; }, lOf(1.1));
        bands[4] = syntheticBandL(31, 1.7, [](int) { return 0.20; }, lOf(2.5));
        const std::array<double, kRubinBands> info = {-1.0, -1.0, 12.0, 50.0, 3.0, -1.0};
        // want: per coordinate, the sum over bands of Swdd - Swd Swd^T / Sw on (physical, f_L), the priors, f_L out
        const int nk = Ny + 3;
        std::vector<double> T(nk * nk, 0.0);
        for (int cc = 0; cc < 2; ++cc)
            for (int bi = 0; bi < 3; ++bi) {
                const AstSums& S = bands[2 + bi];
                int comp[Ny + 1], emb[Ny + 1];
                for (int j = 0; j < Ny; ++j) { comp[j] = j;  emb[j] = j; }
                comp[Ny] = Ny + 1;  emb[Ny] = Ny + bi;
                for (int a = 0; a <= Ny; ++a)
                    for (int m = 0; m <= Ny; ++m)
                        T[emb[a] * nk + emb[m]] += S.Swdd[cc][comp[a] * kD + comp[m]] - S.Swd[cc][comp[a]] * S.Swd[cc][comp[m]] / S.Sw;
            }
        for (int bi = 0; bi < 3; ++bi) T[(Ny + bi) * nk + Ny + bi] += info[2 + bi];
        const std::vector<double> R = schurOut(T, nk, 3);
        std::array<double, Ny * Ny> want{}, got{};
        for (int i = 0; i < Ny * Ny; ++i) want[i] = R[i];
        fisher_detail::rubinBandInformationLensFree(bands, info, got);
        expectClose("h constant in every band (reference degenerate)", got, want);
    }
    return ok;
}

// ---- --lens-light-floor ----
// A dark event whose F146 has no blend (f_S = 1): with no floor Roman has no f_L parameter and its sigma is the known one, bit for
// bit; with a floor F it carries one of at least that width, and Roman's thetaE gets worse as F grows;.
bool checkLensLightFloor()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        const auto run = [&](bool free, double floorF) {
            EventRun r(kEvents[2]);
            r.s->fb[1] = 1.0;  r.s->blend[6] = 1.0;
            r.co->lensLightFree = free;  r.co->lensLightFloor = floorF;
            r.run();
            return r;
        };
        const EventRun known = run(false, 0.0), f0 = run(true, 0.0), f1 = run(true, 0.01), f2 = run(true, 0.1), tiny = run(true, 1e-12);
        bool same = true;
        for (int k = 0; k < Ny; ++k) same = same and sameDouble(f0.co->Erb[SROMAN][k], known.co->Erb[SROMAN][k]);
        const double a = f1.co->Erb[SROMAN][0] / known.co->Erb[SROMAN][0], b = f2.co->Erb[SROMAN][0] / known.co->Erb[SROMAN][0];
        const double t = tiny.co->Erb[SROMAN][0] / known.co->Erb[SROMAN][0];
        if (not same) { std::cerr << "FAIL [lens-light-floor] floor 0: Roman's sigma is not the known one for an unblended F146\n";  ok = false; }
        if (not (a > 1.001 and b > a)) {
            std::cerr << "FAIL [lens-light-floor] Roman sigma(thetaE) free/known " << a << " (floor 0.01), " << b << " (floor 0.1): not 1.001 < a < b\n";
            ok = false;
        }
        if (not (std::fabs(t - 1.0) < 1e-6)) { std::cerr << "FAIL [lens-light-floor] a floor of 1e-12 moves Roman's sigma by " << t - 1.0 << "\n";  ok = false; }
        std::cout << "#   " << (model == SourceModel::ESPL ? "espl" : "pspl") << ": Roman sigma(thetaE) free/known, unblended F146: floor 0 identical, 0.01 -> " << a
                  << ", 0.1 -> " << b << ", 1e-12 -> " << t << "\n";
    }
    vbmSetSourceModel(keep);
    std::cout << "# --lens-light-floor gives an unblended F146 an f_L parameter" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// ---- the per-band decisions, pinned ----

// The chromatic events' forecast (default free-mode astrometry under both priors from the aab9831 source, and the photometric
// Era/condA of the joint and Roman matrices regenerated once from the commit that put the blend prior on F146 too; Rubin's
// photometric rows are the aab9831 ones), hexadecimal floats, generated once from that source with a temporary mode of this file; never regenerated to
// make a check pass. They pin each epoch's own band blend (not r's) in the photometry, the prior on every Rubin band's blend and
// the width of each band's f_L prior (its own f_S, not r's).
const AstroRef kAstroRefFree[] = {
    {0, "chrom_dark_uni", 0, 1, 0x1.a02bf29485c0fp+9, {0x1.57eeacbc25234p-6, 0x1.420717364bc38p-11, 0x1.153bfdacad27ap-11, 0x1.421c05e4b9434p-1}, {{1, 0x1.a02bf29485c0fp+9, {0x1.57eeacbc25234p-6, 0x1.420717364bc38p-11, 0x1.153bfdacad27ap-11, 0x1.421c05e4b9434p-1}}, {1, 0x1.05cbb4279e989p+7, {0x1.6f9b03c17cce1p-3, 0x1.fd67455d13ea8p-11, 0x1.a12726b146cd7p-11, 0x1.258af181f39d3p+2}}, {1, 0x1.894469f3d2983p+7, {0x1.55f8fcc1f2e63p-2, 0x1.3913f42cfb2c8p-10, 0x1.00b0a0a1c0f12p-10, 0x1.49b40f77de103p+2}}}},
    {0, "chrom_dark_uni", 1, 1, 0x1.c4185a4bfa8b8p+7, {0x1.bd59d301e8a5dp-2, 0x1.4f87c37165fap-10, 0x1.132246255d62ep-10, 0x1.541bfb3399c14p+2}, {{1, 0x1.c4185a4bfa8b8p+7, {0x1.bd59d301e8a5dp-2, 0x1.4f87c37165fap-10, 0x1.132246255d62ep-10, 0x1.541bfb3399c14p+2}}, {1, 0x1.c4185a4bfa8b8p+7, {0x1.bd59d301e8a5dp-2, 0x1.4f87c37165fap-10, 0x1.132246255d62ep-10, 0x1.541bfb3399c14p+2}}, {1, 0x1.c4185a4bfa8b8p+7, {0x1.bd59d301e8a5dp-2, 0x1.4f87c37165fap-10, 0x1.132246255d62ep-10, 0x1.541bfb3399c14p+2}}}},
    {0, "chrom_dark_uni", 2, 1, 0x1.0d29b90549583p+10, {0x1.7e19249b1cd86p-6, 0x1.79d3371db7a5bp-11, 0x1.46671501a3b9fp-11, 0x1.457635e19d84ap-1}, {{1, 0x1.0d29b90549583p+10, {0x1.7e19249b1cd86p-6, 0x1.79d3371db7a5bp-11, 0x1.46671501a3b9fp-11, 0x1.457635e19d84ap-1}}, {1, 0x1.1ce3bdc147694p+3, {0x1.01e6dfbe1d662p-2, 0x1.fee0cb0249dc4p-10, 0x1.a73c2ee13645dp-10, 0x1.7a78028cef00fp+3}}, {1, 0x1.6ba536d93613p+2, {0x1.508d17c252f85p-1, 0x1.2c042a101be81p-8, 0x1.14ca9e50cd0dap-8, 0x1.8680cc735368ep+4}}}},
    {0, "chrom_dark_bl", 0, 1, 0x1.27330225590cbp+11, {0x1.f84bedb1ca4b3p-6, 0x1.0f25fc9c421f9p-10, 0x1.d2bd2164712b6p-11, 0x1.45b04695b1b09p-1}, {{1, 0x1.27330225590cbp+11, {0x1.f84bedb1ca4b3p-6, 0x1.0f25fc9c421f9p-10, 0x1.d2bd2164712b6p-11, 0x1.45b04695b1b09p-1}}, {1, 0x1.cd62ef41b2a8bp+9, {0x1.948ebc68bd5ecp-3, 0x1.519c03720258cp-9, 0x1.1430764cd9983p-9, 0x1.26475952bd364p+2}}, {1, 0x1.3bb4fe4b1c15cp+10, {0x1.6a99133c7b263p-2, 0x1.8c463bdf13c85p-9, 0x1.444d521b42028p-9, 0x1.4a363d4d55433p+2}}}},
    {0, "chrom_dark_bl", 1, 1, 0x1.2db2306f01ccdp+11, {0x1.c171cf285f1c6p-2, 0x1.11cf0e3026c3dp-8, 0x1.c0367031a85aap-9, 0x1.55d7e1b7489ddp+2}, {{1, 0x1.2db2306f01ccdp+11, {0x1.c171cf285f1c6p-2, 0x1.11cf0e3026c3dp-8, 0x1.c0367031a85aap-9, 0x1.55d7e1b7489ddp+2}}, {1, 0x1.2db2306f01ccdp+11, {0x1.c171cf285f1c6p-2, 0x1.11cf0e3026c3dp-8, 0x1.c0367031a85aap-9, 0x1.55d7e1b7489ddp+2}}, {1, 0x1.2db2306f01ccdp+11, {0x1.c171cf285f1c6p-2, 0x1.11cf0e3026c3dp-8, 0x1.c0367031a85aap-9, 0x1.55d7e1b7489ddp+2}}}},
    {0, "chrom_dark_bl", 2, 1, 0x1.4b3d78bcaffdp+11, {0x1.0f144faaf3086p-5, 0x1.28523033df257p-10, 0x1.ffe4166a21e3fp-11, 0x1.499be55091f45p-1}, {{1, 0x1.4b3d78bcaffdp+11, {0x1.0f144faaf3086p-5, 0x1.28523033df257p-10, 0x1.ffe4166a21e3fp-11, 0x1.499be55091f45p-1}}, {1, 0x1.220304f36a51cp+4, {0x1.0f0f3e03ffa98p-2, 0x1.d2667fff58b01p-9, 0x1.7b64ec692c3dap-9, 0x1.7a9263de1befp+3}}, {1, 0x1.6eccc7cb9d1efp+2, {0x1.540a125904c15p-1, 0x1.69afd1ed9acf9p-8, 0x1.3ffc84d46711ep-8, 0x1.86a149f9b0aefp+4}}}},
    {0, "chrom_lum_uni", 0, 1, 0x1.265de6ca4ec06p+10, {0x1.5e40cbac7d344p-6, 0x1.50eebd112ec5bp-11, 0x1.187679142f301p-11, 0x1.7327a87542a2bp-2}, {{1, 0x1.265de6ca4ec06p+10, {0x1.5e40cbac7d344p-6, 0x1.50eebd112ec5bp-11, 0x1.187679142f301p-11, 0x1.7327a87542a2bp-2}}, {1, 0x1.2165b3fb7520dp+7, {0x1.553369d9122ap-3, 0x1.d7c49ea5723dcp-11, 0x1.863904296a861p-11, 0x1.4693612e3a725p+0}}, {1, 0x1.f926d2dcfb5b4p+7, {0x1.421b1098779d3p-2, 0x1.374b92aa9f329p-10, 0x1.0152391aa5fe9p-10, 0x1.497070ef02514p+0}}}},
    {0, "chrom_lum_uni", 1, 1, 0x1.52338cb2e32a1p+8, {0x1.a0f212f54c2dep-2, 0x1.677fc422bfe64p-10, 0x1.29461f0acd7b9p-10, 0x1.4a381f141ea1dp+0}, {{1, 0x1.52338cb2e32a1p+8, {0x1.a0f212f54c2dep-2, 0x1.677fc422bfe64p-10, 0x1.29461f0acd7b9p-10, 0x1.4a381f141ea1dp+0}}, {1, 0x1.52338cb2e32a1p+8, {0x1.a0f212f54c2dep-2, 0x1.677fc422bfe64p-10, 0x1.29461f0acd7b9p-10, 0x1.4a381f141ea1dp+0}}, {1, 0x1.52338cb2e32a1p+8, {0x1.a0f212f54c2dep-2, 0x1.677fc422bfe64p-10, 0x1.29461f0acd7b9p-10, 0x1.4a381f141ea1dp+0}}}},
    {0, "chrom_lum_uni", 2, 1, 0x1.9f2017f461661p+10, {0x1.95a70c507f795p-6, 0x1.95b606d1fcf12p-11, 0x1.4c3137f2618cdp-11, 0x1.1a2bc5224105p-1}, {{1, 0x1.9f2017f461661p+10, {0x1.95a70c507f795p-6, 0x1.95b606d1fcf12p-11, 0x1.4c3137f2618cdp-11, 0x1.1a2bc5224105p-1}}, {1, 0x1.4b68763819086p+2, {0x1.f0b1cba3394f3p-3, 0x1.9d4a46e48a1a5p-10, 0x1.78f5cc414e70bp-10, 0x1.08e845c45b02ep+4}}, {1, 0x1.89c96a017e77ap+1, {0x1.3ad45cad0b3f8p-1, 0x1.02c8f69edc7cap-8, 0x1.f017cd24d3d53p-9, 0x1.f0c4df08650a9p+4}}}},
    {0, "chrom_lum_bl", 0, 1, 0x1.3a58b999b99c4p+12, {0x1.419d980619693p-5, 0x1.5c39abec19cefp-10, 0x1.211052988673fp-10, 0x1.c25da6c7c6c28p-2}, {{1, 0x1.3a58b999b99c4p+12, {0x1.419d980619693p-5, 0x1.5c39abec19cefp-10, 0x1.211052988673fp-10, 0x1.c25da6c7c6c28p-2}}, {1, 0x1.159c8ef4ebfffp+10, {0x1.7cc01245c59f8p-3, 0x1.43b36d6114ca5p-9, 0x1.0b9d246dbfceep-9, 0x1.48bda8754d243p+0}}, {1, 0x1.7d6a897f3f14bp+10, {0x1.5d565d63eae88p-2, 0x1.7c95f6590c303p-9, 0x1.3a85f687431bfp-9, 0x1.4c0d40449590ep+0}}}},
    {0, "chrom_lum_bl", 1, 1, 0x1.872ed570c22e9p+11, {0x1.b38284810d9fdp-2, 0x1.107cd1c21a795p-8, 0x1.c263754f3794bp-9, 0x1.4fe49eba0af89p+0}, {{1, 0x1.872ed570c22e9p+11, {0x1.b38284810d9fdp-2, 0x1.107cd1c21a795p-8, 0x1.c263754f3794bp-9, 0x1.4fe49eba0af89p+0}}, {1, 0x1.872ed570c22e9p+11, {0x1.b38284810d9fdp-2, 0x1.107cd1c21a795p-8, 0x1.c263754f3794bp-9, 0x1.4fe49eba0af89p+0}}, {1, 0x1.872ed570c22e9p+11, {0x1.b38284810d9fdp-2, 0x1.107cd1c21a795p-8, 0x1.c263754f3794bp-9, 0x1.4fe49eba0af89p+0}}}},
    {0, "chrom_lum_bl", 2, 1, 0x1.90a580548e674p+12, {0x1.6ea7d07360191p-5, 0x1.8ed4a6588abbfp-10, 0x1.44c94fbdc62ffp-10, 0x1.5e7168e6ed505p-1}, {{1, 0x1.90a580548e674p+12, {0x1.6ea7d07360191p-5, 0x1.8ed4a6588abbfp-10, 0x1.44c94fbdc62ffp-10, 0x1.5e7168e6ed505p-1}}, {1, 0x1.930dd5dbf855p+3, {0x1.08ab63187e27bp-2, 0x1.a774c2fce3316p-9, 0x1.65ab47dda75cfp-9, 0x1.0a2e4b8982405p+4}}, {1, 0x1.a863c38f6902dp+1, {0x1.3e91395b7047dp-1, 0x1.3fbb6bb255b0ap-8, 0x1.24161ca802087p-8, 0x1.f14d53ff90b66p+4}}}},
    {1, "chrom_dark_uni", 0, 1, 0x1.a02bcea79c982p+9, {0x1.57eec790eb00ap-6, 0x1.420707315b208p-11, 0x1.153bede4afaf8p-11, 0x1.421c0a0d2b2dcp-1}, {{1, 0x1.a02bcea79c982p+9, {0x1.57eec790eb00ap-6, 0x1.420707315b208p-11, 0x1.153bede4afaf8p-11, 0x1.421c0a0d2b2dcp-1}}, {1, 0x1.05cbb52f76f0ap+7, {0x1.6f9b22a9b83c3p-3, 0x1.fd67461ed6e36p-11, 0x1.a127271354223p-11, 0x1.258af67b63cb5p+2}}, {1, 0x1.89446aa109ef1p+7, {0x1.55f9136f27fe4p-2, 0x1.3913f44f1127fp-10, 0x1.00b0a0b699ebap-10, 0x1.49b413010c301p+2}}}},
    {1, "chrom_dark_uni", 1, 1, 0x1.c4185a7f9a78cp+7, {0x1.bd59e1d6eb754p-2, 0x1.4f87c37198315p-10, 0x1.132246288e2bap-10, 0x1.541bfc76787p+2}, {{1, 0x1.c4185a7f9a78cp+7, {0x1.bd59e1d6eb754p-2, 0x1.4f87c37198315p-10, 0x1.132246288e2bap-10, 0x1.541bfc76787p+2}}, {1, 0x1.c4185a7f9a78cp+7, {0x1.bd59e1d6eb754p-2, 0x1.4f87c37198315p-10, 0x1.132246288e2bap-10, 0x1.541bfc76787p+2}}, {1, 0x1.c4185a7f9a78cp+7, {0x1.bd59e1d6eb754p-2, 0x1.4f87c37198315p-10, 0x1.132246288e2bap-10, 0x1.541bfc76787p+2}}}},
    {1, "chrom_dark_uni", 2, 1, 0x1.0d2993f6aa71cp+10, {0x1.7e193b3636283p-6, 0x1.79d31a54cadc2p-11, 0x1.4666f99e6f447p-11, 0x1.457638c511cc8p-1}, {{1, 0x1.0d2993f6aa71cp+10, {0x1.7e193b3636283p-6, 0x1.79d31a54cadc2p-11, 0x1.4666f99e6f447p-11, 0x1.457638c511cc8p-1}}, {1, 0x1.1ce3de4f3a55ep+3, {0x1.01e7001b7ec4bp-2, 0x1.fee0e228671aep-10, 0x1.a73c2bec685ebp-10, 0x1.7a7823a15c0a2p+3}}, {1, 0x1.6ba5a7a086d68p+2, {0x1.508d593ba9072p-1, 0x1.2c04401ab49ccp-8, 0x1.14ca9cf37e4a2p-8, 0x1.86810ab832d36p+4}}}},
    {1, "chrom_dark_bl", 0, 1, 0x1.2732aea316de5p+11, {0x1.f84bebebba9f2p-6, 0x1.0f25d47865a03p-10, 0x1.d2bcd8f67f8ddp-11, 0x1.45b0457182c34p-1}, {{1, 0x1.2732aea316de5p+11, {0x1.f84bebebba9f2p-6, 0x1.0f25d47865a03p-10, 0x1.d2bcd8f67f8ddp-11, 0x1.45b0457182c34p-1}}, {1, 0x1.cd62e91e2f632p+9, {0x1.948edee75bd9p-3, 0x1.519c00ff51efcp-9, 0x1.14307440f11acp-9, 0x1.26475e75db9bfp+2}}, {1, 0x1.3bb4fd026e5efp+10, {0x1.6a992ba4d99e2p-2, 0x1.8c463ae337c03p-9, 0x1.444d51488c5afp-9, 0x1.4a3640cc7ae41p+2}}}},
    {1, "chrom_dark_bl", 1, 1, 0x1.2db22ff4c8c1bp+11, {0x1.c171de54eebb4p-2, 0x1.11cf0dea238a1p-8, 0x1.c0366fc156d64p-9, 0x1.55d7e3028427fp+2}, {{1, 0x1.2db22ff4c8c1bp+11, {0x1.c171de54eebb4p-2, 0x1.11cf0dea238a1p-8, 0x1.c0366fc156d64p-9, 0x1.55d7e3028427fp+2}}, {1, 0x1.2db22ff4c8c1bp+11, {0x1.c171de54eebb4p-2, 0x1.11cf0dea238a1p-8, 0x1.c0366fc156d64p-9, 0x1.55d7e3028427fp+2}}, {1, 0x1.2db22ff4c8c1bp+11, {0x1.c171de54eebb4p-2, 0x1.11cf0dea238a1p-8, 0x1.c0366fc156d64p-9, 0x1.55d7e3028427fp+2}}}},
    {1, "chrom_dark_bl", 2, 1, 0x1.4b3cfd5475307p+11, {0x1.0f1443778ca22p-5, 0x1.2851f6d2db68bp-10, 0x1.ffe3af5bb6c87p-11, 0x1.499be165526f4p-1}, {{1, 0x1.4b3cfd5475307p+11, {0x1.0f1443778ca22p-5, 0x1.2851f6d2db68bp-10, 0x1.ffe3af5bb6c87p-11, 0x1.499be165526f4p-1}}, {1, 0x1.22030c2c525fap+4, {0x1.0f0f5fc1d5c64p-2, 0x1.d266852f42f4fp-9, 0x1.7b64e776f5c69p-9, 0x1.7a928524022fap+3}}, {1, 0x1.6ecd2142376ecp+2, {0x1.540a541830b03p-1, 0x1.69afe57b6f195p-8, 0x1.3ffc825254c5ap-8, 0x1.86a1887af951ap+4}}}},
    {1, "chrom_lum_uni", 0, 1, 0x1.265dccf71fa45p+10, {0x1.5e40dcb13162fp-6, 0x1.50eeaa29547f6p-11, 0x1.1876681c27a85p-11, 0x1.7327a5baf9db2p-2}, {{1, 0x1.265dccf71fa45p+10, {0x1.5e40dcb13162fp-6, 0x1.50eeaa29547f6p-11, 0x1.1876681c27a85p-11, 0x1.7327a5baf9db2p-2}}, {1, 0x1.2165b1ba521cfp+7, {0x1.553386221eaeap-3, 0x1.d7c49ca383503p-11, 0x1.863902d0cd4c8p-11, 0x1.469361af23e8dp+0}}, {1, 0x1.f926d7203a18ap+7, {0x1.421b250ef3ba6p-2, 0x1.374b93dd19db7p-10, 0x1.01523a3237c28p-10, 0x1.4970713aab8ddp+0}}}},
    {1, "chrom_lum_uni", 1, 1, 0x1.523390313797ap+8, {0x1.a0f220c829b96p-2, 0x1.677fc5d1f246p-10, 0x1.294620a03668cp-10, 0x1.4a381f53bf614p+0}, {{1, 0x1.523390313797ap+8, {0x1.a0f220c829b96p-2, 0x1.677fc5d1f246p-10, 0x1.294620a03668cp-10, 0x1.4a381f53bf614p+0}}, {1, 0x1.523390313797ap+8, {0x1.a0f220c829b96p-2, 0x1.677fc5d1f246p-10, 0x1.294620a03668cp-10, 0x1.4a381f53bf614p+0}}, {1, 0x1.523390313797ap+8, {0x1.a0f220c829b96p-2, 0x1.677fc5d1f246p-10, 0x1.294620a03668cp-10, 0x1.4a381f53bf614p+0}}}},
    {1, "chrom_lum_uni", 2, 1, 0x1.9f1fe27afd422p+10, {0x1.95a716b51ad88p-6, 0x1.95b5e7f78a00cp-11, 0x1.4c311caa8bdfdp-11, 0x1.1a2bc3378dcdcp-1}, {{1, 0x1.9f1fe27afd422p+10, {0x1.95a716b51ad88p-6, 0x1.95b5e7f78a00cp-11, 0x1.4c311caa8bdfdp-11, 0x1.1a2bc3378dcdcp-1}}, {1, 0x1.4b68943404bbbp+2, {0x1.f0b20b88bf475p-3, 0x1.9d4a441774005p-10, 0x1.78f5c90980a7ep-10, 0x1.08e86ce99e6d8p+4}}, {1, 0x1.89c9d1fb26e34p+1, {0x1.3ad494b3a384fp-1, 0x1.02c8f59fc11f2p-8, 0x1.f017ca382e0e3p-9, 0x1.f0c54308a503fp+4}}}},
    {1, "chrom_lum_bl", 0, 1, 0x1.3a582c1388231p+12, {0x1.419d6ef064e2dp-5, 0x1.5c395946201ddp-10, 0x1.21100cb6aa32cp-10, 0x1.c25d76040ff0ep-2}, {{1, 0x1.3a582c1388231p+12, {0x1.419d6ef064e2dp-5, 0x1.5c395946201ddp-10, 0x1.21100cb6aa32cp-10, 0x1.c25d76040ff0ep-2}}, {1, 0x1.159c897ee3beap+10, {0x1.7cc031fee3821p-3, 0x1.43b36a21bef25p-9, 0x1.0b9d21ccaa9e3p-9, 0x1.48bda8f3e6a93p+0}}, {1, 0x1.7d6a869903f71p+10, {0x1.5d567496b0ffep-2, 0x1.7c95f4ddfddffp-9, 0x1.3a85f55426ae7p-9, 0x1.4c0d40874f05ap+0}}}},
    {1, "chrom_lum_bl", 1, 1, 0x1.872ed47568dc9p+11, {0x1.b38293c75a5e3p-2, 0x1.107cd165d9c23p-8, 0x1.c26374bbf8f36p-9, 0x1.4fe49ee7e7f5dp+0}, {{1, 0x1.872ed47568dc9p+11, {0x1.b38293c75a5e3p-2, 0x1.107cd165d9c23p-8, 0x1.c26374bbf8f36p-9, 0x1.4fe49ee7e7f5dp+0}}, {1, 0x1.872ed47568dc9p+11, {0x1.b38293c75a5e3p-2, 0x1.107cd165d9c23p-8, 0x1.c26374bbf8f36p-9, 0x1.4fe49ee7e7f5dp+0}}, {1, 0x1.872ed47568dc9p+11, {0x1.b38293c75a5e3p-2, 0x1.107cd165d9c23p-8, 0x1.c26374bbf8f36p-9, 0x1.4fe49ee7e7f5dp+0}}}},
    {1, "chrom_lum_bl", 2, 1, 0x1.90a49eaebca78p+12, {0x1.6ea78b647f1e9p-5, 0x1.8ed43165100adp-10, 0x1.44c8ee862bbd1p-10, 0x1.5e713c1d8697dp-1}, {{1, 0x1.90a49eaebca78p+12, {0x1.6ea78b647f1e9p-5, 0x1.8ed43165100adp-10, 0x1.44c8ee862bbd1p-10, 0x1.5e713c1d8697dp-1}}, {1, 0x1.930dd7b680522p+3, {0x1.08ab85e8ec533p-2, 0x1.a774bd7e0f144p-9, 0x1.65ab4482a8242p-9, 0x1.0a2e73cd8ccd6p+4}}, {1, 0x1.a863e271074bcp+1, {0x1.3e9171f85548ep-1, 0x1.3fbb69422e266p-8, 0x1.24161aa5de44dp-8, 0x1.f14db8b8d3ab6p+4}}}},
};
struct PhotRef { int espl; const char* name; int q; int ok; int nep; double cond; double era[Nx]; };
const PhotRef kPhotRefChrom[] = {
    {0, "chrom_dark", 0, 1, 1102, 0x1.6f469758f3fbbp+14, {0x1.bfe1a3434d1bap-8, 0x1.ae52db00ddbb2p-1, 0x1.5d8fbce74332ap-7, 0x1.dff4d6a5c0051p-1, 0x1.3388fe1ce3e26p+0, 0x1.1b159faeb95bp-7, 0x1.3294ff061b882p-9, 0x1.a747f0982dff6p-7, 0x1.07b134f262c5bp-12, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.c1afe18539ec4p-7, 0x1.362e73348bc35p-9, 0x1.ff7555dcdfe8ep-7, 0x1.3718fde01084cp-9, -0x1p+0, -0x1p+0}},
    {0, "chrom_dark", 1, 1, 238, 0x1.2cf7bc861ddb4p+13, {0x1.31cacea4251edp-4, 0x1.8c565dc65a45ep+3, 0x1.87a96e1e1f817p-4, 0x1.8c88ebf3a87ap+3, 0x1.f576c7f93e2b7p+3, 0x1.8d116432b2113p-3, 0x1.390dce83371b4p-9, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.4330547744c06p-3, 0x1.4687017b20365p-9, 0x1.86bbbf58c3e15p-3, 0x1.4eeab930e4598p-9, -0x1p+0, -0x1p+0}},
    {0, "chrom_dark", 2, 1, 864, 0x1.97a0a258c66e8p+16, {0x1.c3c4c3b88a847p-7, 0x1.c2b0e465799ep+0, -0x1p+0, 0x1.fac626871c221p+0, 0x1.444fca3a55374p+1, 0x1.751f0757378aep-7, -0x1p+0, 0x1.7addb97cc0c4cp-6, 0x1.2453733b4b083p-12, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}},
    {0, "chrom_lum", 0, 1, 1102, 0x1.6f469758f3fbbp+14, {0x1.bfe1a3434d1bap-8, 0x1.ae52db00ddbb2p-1, 0x1.5d8fbce74332ap-7, 0x1.dff4d6a5c0051p-1, 0x1.3388fe1ce3e26p+0, 0x1.1b159faeb95bp-7, 0x1.3294ff061b882p-9, 0x1.a747f0982dff6p-7, 0x1.07b134f262c5bp-12, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.c1afe18539ec4p-7, 0x1.362e73348bc35p-9, 0x1.ff7555dcdfe8ep-7, 0x1.3718fde01084cp-9, -0x1p+0, -0x1p+0}},
    {0, "chrom_lum", 1, 1, 238, 0x1.2cf7bc861ddb4p+13, {0x1.31cacea4251edp-4, 0x1.8c565dc65a45ep+3, 0x1.87a96e1e1f817p-4, 0x1.8c88ebf3a87ap+3, 0x1.f576c7f93e2b7p+3, 0x1.8d116432b2113p-3, 0x1.390dce83371b4p-9, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.4330547744c06p-3, 0x1.4687017b20365p-9, 0x1.86bbbf58c3e15p-3, 0x1.4eeab930e4598p-9, -0x1p+0, -0x1p+0}},
    {0, "chrom_lum", 2, 1, 864, 0x1.97a0a258c66e8p+16, {0x1.c3c4c3b88a847p-7, 0x1.c2b0e465799ep+0, -0x1p+0, 0x1.fac626871c221p+0, 0x1.444fca3a55374p+1, 0x1.751f0757378aep-7, -0x1p+0, 0x1.7addb97cc0c4cp-6, 0x1.2453733b4b083p-12, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}},
    {1, "chrom_dark", 0, 1, 1102, 0x1.6f47ad7698a7fp+14, {0x1.bfe15d9afc115p-8, 0x1.ae525a428e925p-1, 0x1.5d8fadad2054fp-7, 0x1.dff44e923c83ap-1, 0x1.3388a6f581fccp+0, 0x1.1b156a111bec3p-7, 0x1.3294feb1599bap-9, 0x1.a74763c84c839p-7, 0x1.07b1301ec576fp-12, 0x1.ffffffffc9cb6p-1, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.c1afd3511aaep-7, 0x1.362e72fa44693p-9, 0x1.ff754f5650a01p-7, 0x1.3718fdc988592p-9, -0x1p+0, -0x1p+0}},
    {1, "chrom_dark", 1, 1, 238, 0x1.2cf7a2176893ep+13, {0x1.31ca4cdaa433bp-4, 0x1.8c5636dac04f1p+3, 0x1.87a91f1b0946ap-4, 0x1.8c88d24526a8ap+3, 0x1.f576a7962d1fbp+3, 0x1.8d1121f5bcc3ap-3, 0x1.390dce76ae163p-9, -0x1p+0, -0x1p+0, 0x1.ffffffffffddfp-1, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.4330115d83628p-3, 0x1.468700384c237p-9, 0x1.86bb703825c32p-3, 0x1.4eeab8c1fa49p-9, -0x1p+0, -0x1p+0}},
    {1, "chrom_dark", 2, 1, 864, 0x1.97a25f5b1ac23p+16, {0x1.c3c4d1f7da13fp-7, 0x1.c2b0a75a587cap+0, -0x1p+0, 0x1.fac5eb70079b3p+0, 0x1.444fa43c40b9bp+1, 0x1.751ed3d6afba9p-7, -0x1p+0, 0x1.7add5faf855f2p-6, 0x1.24536d85dbce6p-12, 0x1.fffffffffc4f5p-1, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}},
    {1, "chrom_lum", 0, 1, 1102, 0x1.6f47ad7698a7fp+14, {0x1.bfe15d9afc115p-8, 0x1.ae525a428e925p-1, 0x1.5d8fadad2054fp-7, 0x1.dff44e923c83ap-1, 0x1.3388a6f581fccp+0, 0x1.1b156a111bec3p-7, 0x1.3294feb1599bap-9, 0x1.a74763c84c839p-7, 0x1.07b1301ec576fp-12, 0x1.ffffffffc9cb6p-1, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.c1afd3511aaep-7, 0x1.362e72fa44693p-9, 0x1.ff754f5650a01p-7, 0x1.3718fdc988592p-9, -0x1p+0, -0x1p+0}},
    {1, "chrom_lum", 1, 1, 238, 0x1.2cf7a2176893ep+13, {0x1.31ca4cdaa433bp-4, 0x1.8c5636dac04f1p+3, 0x1.87a91f1b0946ap-4, 0x1.8c88d24526a8ap+3, 0x1.f576a7962d1fbp+3, 0x1.8d1121f5bcc3ap-3, 0x1.390dce76ae163p-9, -0x1p+0, -0x1p+0, 0x1.ffffffffffddfp-1, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, 0x1.4330115d83628p-3, 0x1.468700384c237p-9, 0x1.86bb703825c32p-3, 0x1.4eeab8c1fa49p-9, -0x1p+0, -0x1p+0}},
    {1, "chrom_lum", 2, 1, 864, 0x1.97a25f5b1ac23p+16, {0x1.c3c4d1f7da13fp-7, 0x1.c2b0a75a587cap+0, -0x1p+0, 0x1.fac5eb70079b3p+0, 0x1.444fa43c40b9bp+1, 0x1.751ed3d6afba9p-7, -0x1p+0, 0x1.7add5faf855f2p-6, 0x1.24536d85dbce6p-12, 0x1.fffffffffc4f5p-1, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0, -0x1p+0}},
};

bool checkChromaticPinned()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    int nAst = 0, nPhot = 0;
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        const int espl = model == SourceModel::ESPL ? 1 : 0;
        for (const int lum : {0, 1})
            for (const int blend : {0, 1}) {
                const auto r = makeLensCase(kLensCases - 2 + lum, true, blend ? 1.0 : kLensLightUniformFrac);
                const std::string nm = std::string(lum ? "chrom_lum" : "chrom_dark") + (blend ? "_bl" : "_uni");
                for (const auto& ref : kAstroRefFree) {
                    if (ref.espl != espl or nm != ref.name) continue;
                    ++nAst;
                    if (not matchesRef(*r->co, ref)) {
                        std::cerr << "FAIL [chromatic-pinned] astrometry " << nm << (espl ? " (espl)" : " (pspl)") << " partition " << ref.q
                                  << " is not bit-identical to the aab9831 forecast\n";
                        ok = false;
                    }
                }
                if (blend) continue;                                  // the photometry does not read the lens-light prior
                for (const auto& ref : kPhotRefChrom) {
                    if (ref.espl != espl or std::strcmp(ref.name, lum ? "chrom_lum" : "chrom_dark") != 0) continue;
                    ++nPhot;
                    const covarian& c = *r->co;
                    bool same = c.okA[ref.q] == ref.ok and c.nepochA[ref.q] == ref.nep and sameDouble(c.condA[ref.q], ref.cond);
                    for (int k = 0; k < Nx; ++k) same = same and sameDouble(c.Era[ref.q][k], ref.era[k]);
                    if (not same) {
                        std::cerr << "FAIL [chromatic-pinned] photometry " << ref.name << (espl ? " (espl)" : " (pspl)") << " partition " << ref.q
                                  << " is not bit-identical to the aab9831 forecast\n";
                        ok = false;
                    }
                }
            }
    }
    vbmSetSourceModel(keep);
    if (nAst != 24 or nPhot != 12) { std::cerr << "FAIL [chromatic-pinned] compared " << nAst << " astrometric and " << nPhot << " photometric rows\n"; ok = false; }
    std::cout << "# chromatic events, default free-mode astrometry (both priors) and photometry equal the aab9831 forecast bit for bit: "
              << nAst << " + " << nPhot << " rows, " << (ok ? "ok\n" : "FAIL\n");
    return ok;
}

// Every Rubin visit in i with i's own blend fraction and baseline (and r's left different) is the all-r event whose r pair is
// those numbers (to rounding): each epoch reads its band's pair, not r's.
bool checkRubinBandBlendUsed()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        EventRun a(kEvents[2]), b(kEvents[2]);
        a.s->blend[2] = 0.35;  a.s->magb[2] = 20.2;  a.s->fb[0] = 0.35;  a.s->mbs[0] = 20.2;
        b.s->blend[3] = 0.35;  b.s->magb[3] = 20.2;                       // r keeps 0.50 and its baseline
        for (int i = 0; i < b.ndw; ++i) if (b.l->tele[i] == 0) b.l->filt[i] = 3;
        a.run();  b.run();
        for (const int q : {SJOINT, SRUBIN}) {
            // Equal to rounding (the two fits order their parameters differently, so the inversion's last bits differ, ~1e-13).
            const auto close = [](double x, double y) { return std::fabs(x - y) <= 1e-9 * std::fabs(y); };
            bool same = a.co->okA[q] == b.co->okA[q];
            for (const int k : {0, 1, 3, 4, 5, 9}) same = same and close(a.co->Era[q][k], b.co->Era[q][k]);
            same = same and close(a.co->Era[q][fbIndexOfFilter(2)], b.co->Era[q][fbIndexOfFilter(3)])
                        and close(a.co->Era[q][mbsIndexOfFilter(2)], b.co->Era[q][mbsIndexOfFilter(3)]);
            if (not same or not a.co->okA[q]) {
                std::cerr << "FAIL [band-blend-used] partition " << q << (model == SourceModel::ESPL ? " (espl)" : " (pspl)")
                          << ": all visits in i with i's pair differ from the all-r event with that pair\n";
                ok = false;
            }
        }
    }
    vbmSetSourceModel(keep);
    std::cout << "# all Rubin visits in i (own blend) = the all-r event with that pair (to 1e-9)" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// The prior on the blend fraction (sigma = FB_PRIOR_SIG) adds exactly 1/sigma^2 to the fb diagonal of every Rubin band
// in the fit, r included, in the joint and Rubin-only matrices, and to Roman's F146 in the joint and Roman-only ones; and
// nothing to a band in a matrix that does not hold it.
bool checkBlendPriorEveryBand()
{
    bool ok = true;
    EventRun on(kEvents[2]), off(kEvents[2]);
    for (EventRun* e : {&on, &off}) makeChromatic(*e->s, *e->l, e->ndw);
    off.co->fbPriorSig = std::numeric_limits<double>::infinity();
    on.run();  off.run();
    const double want = 1.0 / (FB_PRIOR_SIG * FB_PRIOR_SIG);
    int nb = 0;
    for (int b = 0; b < kRubinBands; ++b) {
        if (on.co->nRubinBand[b] == 0) continue;
        ++nb;
        const int j = fbIndexOfFilter(b);
        for (const int q : {SJOINT, SRUBIN}) {
            const double d = gsl_matrix_get(on.co->inputA[q].get(), j, j) - gsl_matrix_get(off.co->inputA[q].get(), j, j);
            if (not (std::fabs(d - want) <= 1e-12 * want)) {
                std::cerr << "FAIL [blend-prior] band " << b << " partition " << q << ": prior adds " << d << ", not " << want << "\n";
                ok = false;
            }
        }
    }
    const int jR = fbIndexOfFilter(6);
    for (const int q : {SJOINT, SROMAN}) {
        const double d = gsl_matrix_get(on.co->inputA[q].get(), jR, jR) - gsl_matrix_get(off.co->inputA[q].get(), jR, jR);
        if (not (std::fabs(d - want) <= 1e-12 * want)) {
            std::cerr << "FAIL [blend-prior] F146 partition " << q << ": prior adds " << d << ", not " << want << "\n";
            ok = false;
        }
    }
    // a band in a matrix that does not hold it: Rubin's bands in the Roman matrix, F146 in the Rubin one
    for (int b = 0; b < kRubinBands; ++b) {
        const int j = fbIndexOfFilter(b);
        if (gsl_matrix_get(on.co->inputA[SROMAN].get(), j, j) != gsl_matrix_get(off.co->inputA[SROMAN].get(), j, j)) {
            std::cerr << "FAIL [blend-prior] the prior changed Rubin band " << b << "'s diagonal in the Roman matrix\n";  ok = false;
        }
    }
    if (gsl_matrix_get(on.co->inputA[SRUBIN].get(), jR, jR) != gsl_matrix_get(off.co->inputA[SRUBIN].get(), jR, jR)) {
        std::cerr << "FAIL [blend-prior] the prior changed F146's diagonal in the Rubin matrix\n";  ok = false;
    }
    if (nb != 3) { std::cerr << "FAIL [blend-prior] " << nb << " bands in the fit, not r, i, z\n"; ok = false; }
    std::cout << "# the blend prior adds 1/sigma^2 to r, i, z in J and L, to F146 in J and R, and nothing elsewhere" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

// In free mode the joint matrix is still the Rubin one plus the Roman one in every variant (each telescope's fractions are
// its own nuisance parameters), and the joint sigma is no worse than either.
bool checkLensLightAdditivity()
{
    bool ok = true;
    const SourceModel keep = vbmSourceModel();
    for (const SourceModel model : {SourceModel::PSPL, SourceModel::ESPL}) {
        vbmSetSourceModel(model);
        for (const double frac : {kUniformFrac, 1.0}) {
            const auto r = makeLensCase(kLensCases - 1, true, frac);
            ok = astroStandardChecks(*r->co, "lens-light-additivity") and ok;
            double worst = 0.0;
            for (int v = 0; v < NAVAR; ++v)
                for (int jm = 0; jm < Ny * Ny; ++jm) {
                    const double j = r->co->FBV[v][SJOINT][jm];
                    const double sum = r->co->FBV[v][SRUBIN][jm] + r->co->FBV[v][SROMAN][jm];
                    worst = std::max(worst, std::fabs(j - sum) / (std::fabs(j) + std::fabs(sum) + 1e-300));
                }
            if (not (worst < 1e-12)) {
                std::cerr << "FAIL [lens-light-additivity] FBV[J] != FBV[L] + FBV[R], worst relative " << worst << "\n";
                ok = false;
            }
            std::cout << "#   " << (model == SourceModel::ESPL ? "espl" : "pspl") << ", prior x" << frac
                      << ": FBV[J] vs FBV[L] + FBV[R] in all three variants, worst relative " << worst << "\n";
        }
    }
    vbmSetSourceModel(keep);
    std::cout << "# free lens light: J = L + R" << (ok ? "  -> ok\n" : "  -> FAIL\n");
    return ok;
}

} // namespace

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
    if (argc > 1 && std::string(argv[1]) == "--blend") return runBlend();
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
    std::cout << "\n# --- Rubin's bands ---\n";
    if (!checkPriorOffReproducesToday()) ++failures;
    if (!checkBandAtBaseline())          ++failures;
    if (!checkFewBandVisits())           ++failures;
    if (!checkRubinBandsAllSparse())     ++failures;
    if (!checkChromaticEvent())          ++failures;
    if (!checkROmittedSentinel())        ++failures;
    std::cout << "\n# --- Astrometric dilution per band ---\n";
    if (!checkAstroAllRUnchanged())      ++failures;
    if (!checkAstroBandSplitNoGain())    ++failures;
    if (!checkAstroChromaticDilution())  ++failures;
    if (!checkAstroSingleOtherBand())    ++failures;
    if (!checkAstroShareCancels())       ++failures;
    if (!checkAstroSparseBandDropped())  ++failures;
    if (!checkAstroRubinAllSparse())     ++failures;
    if (!checkAstroClosedForm())         ++failures;
    std::cout << "\n# --- The lens's light fraction as a free parameter of the astrometric fit ---\n";
    if (!checkLensLightKnownUnchanged()) ++failures;
    if (!checkLensLightFreeFunction())   ++failures;
    if (!checkLensLightDerivative())     ++failures;
    if (!checkLensLightFreeNoGain())     ++failures;
    if (!checkLensLightPriorMonotone())  ++failures;
    if (!checkLensLightAdditivity())     ++failures;
    if (!checkChromaticPinned())         ++failures;
    if (!checkRubinBandBlendUsed())      ++failures;
    if (!checkBlendPriorEveryBand())     ++failures;
    if (!checkLensLightRoman())          ++failures;
    if (!checkLensLightFreeBranches())   ++failures;
    if (!checkLensLightFloor())          ++failures;

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
