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

int main()
{
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto co = std::make_unique<covarian>();

    std::cout << "# Fisher-matrix fixture -- synthetic events, no data files required\n"
              << "# Nx=" << Nx << " (photometric)   Ny=" << Ny << " (astrometric)\n"
              << "# sigma columns are ABSOLUTE 1-sigma forecast uncertainties from co.Era[]/co.Erb[]\n"
              << "#\n";

    std::cout << std::left << std::setw(16) << "# event"
              << std::right
              << std::setw(8)  << "tE"
              << std::setw(9)  << "t0"
              << std::setw(7)  << "ndw"
              << std::setw(7)  << "nL"
              << std::setw(7)  << "nR"
              << std::setw(13) << "sig_u0"
              << std::setw(13) << "sig_tE"
              << std::setw(13) << "sig_fb"
              << std::setw(13) << "sig_piE"
              << std::setw(13) << "sig_xi";
    if (Nx > 5) std::cout << std::setw(13) << "sig_t0";
    std::cout << std::setw(13) << "sig_tetE" << "\n";

    for (const auto& ev : kEvents) {
        setupStatic(*s, *l);
        l->tE  = ev.tE;
        l->t0  = ev.t0;
        l->u0  = ev.u0;
        l->piE = ev.piE;

        int nL = 0, nR = 0;
        const int ndw = buildLightCurve(*s, *l, *as, nL, nR);

        // Rebuild the true-parameter state before handing off: buildLightCurve left
        // lightcurve() evaluated at the final epoch, and FisherM perturbs from whatever the
        // struct currently holds.
        l->tE  = ev.tE;
        l->t0  = ev.t0;
        l->u0  = ev.u0;
        l->piE = ev.piE;
        s->xi  = kXi;
        s->fb[0] = kFbRubin;
        s->fb[1] = kFbRoman;

        // Guard against an epoch landing exactly on t0. At such an epoch the source-motion
        // terms mus1*(timh - t0) and mus2*(timh - t0) vanish identically, so perturbing mus1
        // (or mus2) leaves BOTH modelled coordinates bit-identical and FisherM's astrometric
        // branch aborts on
        //     CHECK(!(s.pos1c == l.soux[i] && s.pos2c == l.souy[i]))
        // The t0 values above are deliberately off-grid to avoid this. If you change the
        // cadence or the events, this guard tells you what went wrong instead of leaving you
        // with an opaque abort. See OPEN_ITEMS.md -- the underlying fragility is in the
        // production CHECK, not in this fixture.
        for (int i = 0; i < ndw; ++i) {
            if (l->timn[i] == l->t0) {
                std::cerr << "fixture error: event '" << ev.name << "' has an epoch exactly at "
                          << "t0=" << l->t0 << "; nudge t0 off the cadence grid.\n";
                return 1;
            }
        }

        FisherM(*s, *l, *as, *co, ndw);
        ErrorCal(*co, *l, *s);

        std::cout << std::left << std::setw(16) << ev.name
                  << std::right << std::fixed
                  << std::setw(8) << std::setprecision(1) << ev.tE
                  << std::setw(9) << std::setprecision(1) << ev.t0
                  << std::setw(7) << ndw
                  << std::setw(7) << nL
                  << std::setw(7) << nR
                  << std::scientific << std::setprecision(4)
                  << std::setw(13) << co->Era[0]
                  << std::setw(13) << co->Era[1]
                  << std::setw(13) << co->Era[2]
                  << std::setw(13) << co->Era[3]
                  << std::setw(13) << co->Era[4];
        if (Nx > 5) std::cout << std::setw(13) << co->Era[5];
        std::cout << std::setw(13) << co->Erb[0] << "\n";
    }

    std::cout << "#\n# done\n";
    return 0;
}
