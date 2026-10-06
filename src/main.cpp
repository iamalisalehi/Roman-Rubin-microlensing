// The simulator: sightline scan -> star draws -> light curves -> detection -> Fisher forecast -> output.
#include "common.h"
#include "types.h"
#include "util/random.h"
#include "run/config.h"
#include "run/outputs.h"
#include "run/sample_dump.h"
#include "run/histograms.h"
#include "galaxy/catalogue.h"
#include "galaxy/extinction.h"
#include "galaxy/density.h"
#include "galaxy/kinematics.h"
#include "events/source.h"
#include "events/lens.h"
#include "events/lightcurve.h"
#include "surveys/visits.h"
#include "surveys/footprints.h"
#include "surveys/schedule.h"
#include "surveys/noise.h"
#include "fisher/fisher.h"
#include "fisher/linalg.h"
#include <sstream>   // provenance block
#include <utility>   // std::pair, for the Roman field list

time_t _timeNow;
unsigned int _randVal;
unsigned int _dummyVal;
FILE * _randStream;

///==============================================================//
///                                                              //                                                    /
///                  Main program                                //
///                                                              //
///==============================================================//

int main(int argc, char** argv) {
    // NOTE: srand(time(0)) used to be called here. Nothing in this project ever
    // calls rand() -- the RNG is the seeded mt19937_64 in Bulge.h -- so it did
    // nothing except make the run look clock-seeded, which is the opposite of
    // the reproducibility the provenance block below is for.

    RunConfig cfg;
    bool strideGiven = false;
    for (int a = 1; a < argc; ++a) {
        const std::string arg = argv[a];
        auto need = [&](const char* what) -> const char* {
            if (a + 1 >= argc) {
                std::cerr << "ERROR: " << what << " needs a value\n";
                std::exit(2);
            }
            return argv[++a];
        };
        if      (arg == "--stride") { cfg.stride = std::atoi(need("--stride")); strideGiven = true; }
        else if (arg == "--stride-roman") cfg.strideRoman = std::atoi(need("--stride-roman"));
        else if (arg == "--dchi-det") cfg.dchiDet = std::atof(need("--dchi-det"));
        else if (arg == "--events")  cfg.iconTarget  = std::atoi(need("--events"));
        else if (arg == "--lenses")  cfg.nlensTarget = std::atoi(need("--lenses"));
        else if (arg == "--nerr")    cfg.nerrTarget  = std::atof(need("--nerr"));
        else if (arg == "--maxdraws") cfg.maxDraws    = std::atof(need("--maxdraws"));
        else if (arg == "--stub")    cfg.stubPatch   = true;
        else if (arg == "--start-index") cfg.startIndex = std::atol(need("--start-index"));
        else if (arg == "--end-index")   cfg.endIndex   = std::atol(need("--end-index"));
        else if (arg == "--seed")        cfg.seedBase   = std::strtoull(need("--seed"), nullptr, 10);
        else if (arg == "--dry-run") cfg.dryRun      = true;
        else if (arg == "--no-satellite-parallax") cfg.noSatPar = true;
        else if (arg == "--pair-satellite") cfg.pairSat = true;
        else if (arg == "--dump-samples") cfg.dumpSpec = need("--dump-samples");
        else if (arg == "--population") {
            const std::string want = need("--population");
            const LensPopulation* found = nullptr;
            for (const auto& p : POPULATIONS)
                if (want == p.name) { found = &p; break; }
            if (!found) {
                std::cerr << "ERROR: unknown population '" << want << "'. Known:\n";
                for (const auto& p : POPULATIONS)
                    std::cerr << "   " << std::left << std::setw(16) << p.name
                              << p.mlMin << " - " << p.mlMax << " Msun   -> "
                              << "test" << p.tag << ".dat   (" << p.note << ")\n";
                std::exit(2);
            }
            gPop = found;
        }
        else if (arg == "--help")  { printUsage(argv[0]); return 0; }
        else {
            std::cerr << "ERROR: unknown option '" << arg << "'\n";
            printUsage(argv[0]);
            return 2;
        }
    }
    // The stub patch is 0.1x0.1 deg. The default stride of 10 (0.20 deg) would step
    // clean over it and leave a single sightline, so --stub without an explicit
    // --stride falls back to the native dd grid -- which is what the pre-Step-4 loop
    // used, and the only way --stub reproduces those numbers.
    if (cfg.stubPatch and not strideGiven) cfg.stride = 1;

    if (cfg.pairSat and cfg.noSatPar) {
        std::cerr << "ERROR: --pair-satellite compares Roman at L2 against Roman at Earth, but "
                  << "--no-satellite-parallax\n       has already put it at Earth. There "
                  << "would be nothing to compare.\n";
        return 1;
    }
    if (not (cfg.dchiDet > 0.0) or not std::isfinite(cfg.dchiDet)) {
        std::cerr << "ERROR: --dchi-det (" << cfg.dchiDet << ") must be finite and positive. "
                  << "It is a delta-chi2 detection bar; a non-positive bar would declare every "
                  << "event detected.\n";
        return 1;
    }
    if (cfg.stride < 1) {
        std::cerr << "ERROR: --stride must be >= 1\n";
        return 2;
    }
    if (cfg.iconTarget < 1 or cfg.nlensTarget < 0 or cfg.nerrTarget < 0.0) {
        std::cerr << "ERROR: --events must be >= 1, --lenses and --nerr >= 0\n";
        return 2;
    }
    // A cap below the event budget would stop every sightline early, which is not a cap
    // but a silent redefinition of the budget.
    if (cfg.endIndex >= 0 and cfg.endIndex <= cfg.startIndex) {
        std::cerr << "ERROR: --end-index (" << cfg.endIndex << ") must exceed --start-index ("
                  << cfg.startIndex << ").\n";
        return 2;
    }
    if (cfg.startIndex < 0) {
        std::cerr << "ERROR: --start-index (" << cfg.startIndex << ") cannot be negative.\n";
        return 1;
    }
    if (cfg.maxDraws < double(cfg.iconTarget)) {
        std::cerr << "ERROR: --maxdraws (" << cfg.maxDraws << ") is below --events ("
                  << cfg.iconTarget << "); no sightline could reach its budget\n";
        return 2;
    }

    // The GBTDS detector layout (Deviation 69), read first: the footprint grid is checked
    // against the real detector size, and the scan region is built from the field outline.
    const GbtdsLayout gl = readGbtdsLayout();
    std::cout << "**** GBTDS layout read: " << GBTDS_NSCA << " detectors x " << GBTDS_NLAYOUT
              << " rolls, detector side " << gl.scaSide << " deg, field reach " << gl.rField
              << " deg ****\n";
    const double gridStep = cfg.stride * dd;

    // Step E1. strideRoman = 0 means "not given": the footprint grid is the coarse grid, kSub
    // = 1 -- unless that grid is coarser than one Roman detector (Deviation 69), in which case
    // it is refined to the largest divisor of --stride that is not, so that which sightlines
    // fall on a detector and which in a chip gap is sampled at all.
    if (cfg.strideRoman == 0) {
        cfg.strideRoman = cfg.stride;
        if (not cfg.stubPatch and cfg.stride * dd > gl.scaSide) {
            for (int k = cfg.stride; k >= 1; --k)
                if (cfg.stride % k == 0 and k * dd <= gl.scaSide) { cfg.strideRoman = k; break; }
            std::cout << "NOTE: --stride " << cfg.stride << " (" << cfg.stride * dd << " deg) is "
                      << "coarser than one Roman detector (" << gl.scaSide << " deg); using "
                      << "--stride-roman " << cfg.strideRoman << " inside the footprint.\n";
        }
    }
    if (cfg.strideRoman < 1 or cfg.strideRoman > cfg.stride) {
        std::cerr << "ERROR: --stride-roman (" << cfg.strideRoman << ") must be between 1 and "
                  << "--stride (" << cfg.stride << "). It refines the grid inside Roman's "
                  << "footprint; it cannot coarsen it.\n";
        return 2;
    }
    // Must divide exactly, or the fine cells do not tile the coarse ones and the sky-area
    // weights stop summing to the scanned area -- which is the one thing this whole
    // stratification has to preserve.
    if (cfg.stride % cfg.strideRoman != 0) {
        std::cerr << "ERROR: --stride-roman (" << cfg.strideRoman << ") must divide --stride ("
                  << cfg.stride << ") exactly, so that each coarse cell is a whole number of "
                  << "fine cells and the area weights sum to the scanned area.\n";
        return 2;
    }
    const int    kSub     = cfg.stride / cfg.strideRoman; // fine cells per coarse cell, per axis
    const double fineStep = cfg.strideRoman * dd;

    // The footprint grid must not step over whole detectors: with a step wider than one, the
    // grid's point-sampled Roman area is noise and a field can be missed outright. The per-field
    // guard below then re-checks every (field, roll) placement against the detectors.
    if (fineStep > gl.scaSide and not cfg.stubPatch) {
        std::cerr << "ERROR: --stride" << (kSub > 1 ? "-roman " : " ") << cfg.strideRoman
                  << " gives a footprint grid step of " << fineStep << " deg, wider than one "
                  << "Roman detector (" << gl.scaSide << " deg). Use --stride-roman <= "
                  << int(gl.scaSide / dd + 1e-9) << ".\n";
        return 2;
    }

    // --------------------- Allocate objects ------------------------
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto cm = std::make_unique<CMD>();
//    auto ga = std::make_unique<galactic>();
    auto ex = std::make_unique<extin>();
    auto ls = std::make_unique<lsst>();
    auto ro = std::make_unique<roman>();
    auto co = std::make_unique<covarian>();
    // Step H3's second forecast: the same event with the satellite offset zeroed. Allocated
    // once beside `co` rather than per event -- covarian owns several vectors, and building
    // one per detection would cost more than the Fisher call it serves.
    auto coNS = std::make_unique<covarian>();

    // Step H1: Roman's observer position. satScale multiplies L2_OFFSET_AU inside
    // lightcurve(), so 0 puts Roman back at the centre of the Earth -- the pre-H1 behaviour,
    // and the "off" run of Step H3's experiment.
    as->satScale = cfg.noSatPar ? 0.0 : 1.0;
    
    std::vector<EventRecord> records;
    records.reserve(1000); // rough upper bound on icon per field

    // --------------------- Read BulgeBaseline.dat ------------------
    std::ifstream fil(PATH_BULGE_BASELINE);
    if (!fil) {
        std::cerr << "Cannot read BulgeBaseline.dat\n";
        return 1;
    }

    int ID;
    double TV, airm, seeingVal, skyB, texp;
//    int tmpFilter;
    std::string header;

    std::getline(fil, header);   // Skip the header line
    for (int i = 0; i < Nl; ++i) {
        fil >> ID >> ls->RA[i] >> ls->DEC[i] >> ls->l[i] >> ls->b[i]
            >> ls->tim[i] >> ls->filter[i] >> airm >> seeingVal >> skyB
            >> TV >> ls->sig5[i] >> texp >> ls->dist[i] >> ls->rot[i];
//if (i == 0) cout << ID << endl;
        // A failed extraction is a silent no-op that leaves this row zero-initialised,
        // and zeros pass every CHECK below: (l,b)=(0,0) is inside the bulge region,
        // tim=0 is inside [0,Tobs], filter=0 is a valid u-band index. So the stream
        // state is the only thing that can catch a short or malformed baseline.
        if (!fil) {
            std::cerr << "BulgeBaseline.dat: read failed at row " << i << " of " << Nl
                      << ". File has fewer rows than Nl, or contains a stray header "
                      << "(append-mode duplicate). Regenerate with readbaselineBulge.py.\n";
            return 1;
        }
        CHECK(airm >= 0.0);
        CHECK(ls->filter[i] >= 0);
        CHECK(ls->filter[i] < 6);
        CHECK(ls->tim[i] >= 0.0);
        CHECK(ls->tim[i] <= Tobs);
        CHECK(texp >= 0.0);
   }
    fil.close();
    std::cout << "**** File BulgeBaseline.dat was read ****\n";

    // --------------------- Read sigmaA_LSST.txt -------------------
    // Atrometric error?
    fil.open(PATH_SIGMA_A_LSST);
    if (!fil) { std::cerr << "Cannot read sigmaA_LSST.txt\n"; return 1; }

    for (int i = 0; i < Na; ++i) {
        fil >> ls->mag[i] >> ls->err[i];

        CHECK(ls->mag[i] >= 16.0);
        CHECK(ls->mag[i] <= 25.0);
        CHECK(ls->err[i] >= 0.2);
        CHECK(ls->err[i] <= 5.0);
    }
    fil.close();
    std::cout << "**** File sigmaA_LSST.txt was read ****\n";

    // --------------------- Read sigma_Roman.txt ---------------------
    // Photometric error
    fil.open(PATH_SIGMA_ROMAN);
    if (!fil) { std::cerr << "Cannot read sigma_roman.txt\n"; return 1; }

    for (int i = 0; i < NaRoman; ++i) {
        fil >> ro->mag[i] >> ro->err[i];
        if (!fil) { std::cerr << "sigma_roman.txt: read failed at row " << i << "\n"; return 1; }
        CHECK(ro->mag[i] > 12.0);
        CHECK(ro->err[i] > 0.0);
        if (i > 0) CHECK(ro->mag[i] > ro->mag[i - 1]);
    }
    fil.close();
    {
        // Deviation 72: anchor the curve's 5-sigma point to ROMAN_DEPTH5_AB (see config/parameters.h).
        const double e5 = 1.0857 / 5.0;
        double m5 = -1.0;
        for (int i = 1; i < NaRoman; ++i)
            if (ro->err[i - 1] < e5 and ro->err[i] >= e5) {
                const double f = (std::log(e5) - std::log(ro->err[i - 1]))
                               / (std::log(ro->err[i]) - std::log(ro->err[i - 1]));
                m5 = ro->mag[i - 1] + f * (ro->mag[i] - ro->mag[i - 1]);
                break;
            }
        if (m5 < 0.0) { std::cerr << "sigma_roman.txt never reaches 5 sigma\n"; return 1; }
        const double shift = ROMAN_DEPTH5_AB - m5;
        for (int i = 0; i < NaRoman; ++i) {
            const double phot2 = std::max(ro->err[i] * ro->err[i] - ROMAN_PHOT_FLOOR * ROMAN_PHOT_FLOOR, 0.0);
            ro->mag[i] += shift;
            ro->err[i]  = std::sqrt(phot2 + ROMAN_PHOT_FLOOR * ROMAN_PHOT_FLOOR);
        }
        std::cout << "**** File sigma_roman.txt was read: Penny+2019 5-sigma point " << m5
                  << " AB shifted by " << shift << " mag to " << ROMAN_DEPTH5_AB << " ****\n";
    }
 
    // --------------------- Read RomanBaseline.dat -------------------
    // Written by Baseline/generateRomanBaseline.py (adopted GBTDS layout, Deviation 69):
    // ID RA DEC l b time sig5 field layout.
    fil.open(PATH_ROMAN_BASELINE);
    if (!fil) {
        std::cerr << "Cannot read RomanBaseline.dat\n";
        return 1;
    }
    std::getline(fil, header); // skip header line
    for (int i = 0; i < NlRoman; ++i) {
        fil >> ID >> ro->RA[i] >> ro->DEC[i] >> ro->l[i] >> ro->b[i] >> ro->tim[i] >> ro->sig5[i]
            >> ro->field[i] >> ro->layout[i];
        // Same silent zero-fill failure mode as the Rubin read above.
        if (!fil) {
            std::cerr << "RomanBaseline.dat: read failed at row " << i << " of " << NlRoman
                      << ". Regenerate with generateRomanBaseline.py and update NlRoman.\n";
            return 1;
        }

        CHECK(ro->tim[i] >= 0.0);
        CHECK(ro->tim[i] <= Tobs);
        CHECK(ro->layout[i] >= 0 and ro->layout[i] < GBTDS_NLAYOUT);
        CHECK(ro->field[i] >= 0 and ro->field[i] < 6);
        // Deliberately no filter CHECK here — Roman is single-band (F146, index 6) for now.
    }
    fil.close();
    std::cout << "**** File RomanBaseline.dat was read ****\n";

    // --------------------- Roman season geometry (Step D1) ---------------
    // Recovered from the epoch times just read, not restated from the generator.
    // Needed per event to place t0 relative to Roman's observing windows -- the
    // independent variable of the gap-filling result.
    const RomanSchedule sched = buildRomanSchedule(*ro);

    // The clustering is only meaningful if the two spacing populations it separates
    // are genuinely separated: every within-season spacing below the threshold, every
    // between-season gap above it. If a future schedule ever samples a season more
    // sparsely than SEASON_GAP_MIN_DAYS, or packs seasons closer together than it, the
    // seasons come out wrong while dt_edge and t0zone still look perfectly reasonable.
    // Refuse to run rather than emit gap geometry that is quietly fiction.
    if (sched.seasons.size() < 2
        or sched.maxInSeasonSpacing >= SEASON_GAP_MIN_DAYS
        or sched.minSeasonGap       <= SEASON_GAP_MIN_DAYS
        or sched.minSeasonLength    <= 0.0) {
        std::cerr << "FATAL: cannot separate Roman observing seasons from inter-season gaps.\n"
                  << "       SEASON_GAP_MIN_DAYS = " << SEASON_GAP_MIN_DAYS << " d, but this "
                  << "schedule has\n"
                  << "       max in-season spacing " << sched.maxInSeasonSpacing << " d and "
                  << "min inter-season gap " << sched.minSeasonGap << " d,\n"
                  << "       shortest season " << sched.minSeasonLength << " d\n"
                  << "       (" << sched.seasons.size() << " season(s) found over days "
                  << sched.missionStart << " - " << sched.missionEnd << ").\n"
                  << "       Retune SEASON_GAP_MIN_DAYS in config/parameters.h against the cadence in\n"
                  << "       Baseline/generateRomanBaseline.py before trusting dt_edge/t0zone.\n";
        return 2;
    }
    std::cout << "**** Roman schedule: " << sched.seasons.size() << " seasons over days "
              << sched.missionStart << " - " << sched.missionEnd << " ****\n";

    // --------------------- Read extinction ------------------------
    readExtinction(*ex, PATH_EXT_TABLES);
    readLensML(PATH_LENS_ML);   // luminous lenses (Deviation 74)
    readLsstCamMap(PATH_LSSTCAM_FOV);   // Rubin's footprint (Deviation 80)
    {
        // The galactic->ICRS conversion the coverage test relies on, checked against OpSim's own
        // RA/Dec for every Rubin visit (whose l, b readbaselineBulge.py derived with astropy).
        double worst = 0.0;
        for (int i = 0; i < Nl; ++i) {
            double ra, dec;
            galToIcrs(ls->l[i], ls->b[i], ra, dec);
            double dra = std::fabs(ra - ls->RA[i]); if (dra > 180.0) dra = 360.0 - dra;
            worst = std::max(worst, std::hypot(dra * std::cos(dec * M_PI / 180.0), dec - ls->DEC[i]));
        }
        std::cout << "galactic->ICRS check over " << Nl << " Rubin visits: worst " << worst * 3600.0 << " arcsec\n";
        if (worst > 1e-4) { std::cerr << "ERROR: galToIcrs disagrees with OpSim by " << worst << " deg\n"; return 1; }
    }

    // --------------------- Call read_cmd --------------------------
    read_cmd(*cm);
    std::cout << "******* read_cmd was done ************" << std::endl;

    int    dclsEvent = DET_NONE;
    // Bin indices for the seven detection-efficiency axes. gg (tE) was the only one ever
    // computed; the other six were commented out here and at their call site, which is why
    // every efficiency column but tE has been a column of zeros (Deviation 46).
    int    save = 0, flagL, gg = -1, ss = 0, qq = 0, ww = 0, vv = 0, zz = 0, pp = 0;
    int    nri = -1, nde = -1, icon;
    int    nlens;// hh; // nde1, nri1,
    std::array<int, NDETCLASS> nDetClass{}; // per-field detection-taxonomy counts (DetClass, Bulge.h)
    // Run-wide totals. The joint-only class -- events neither telescope can find alone but the
    // combination can -- is the strongest evidence for the joint fit and is expected to be rare,
    // so it needs statistics pooled over every field, not per-field counts that are individually
    // too small to quote. Broken down by tE bin as well, since the whole science case is that
    // the gain is tE-dependent.
    std::array<long, NDETCLASS> NDetClassTot{};
    long nDchiMismatch = 0;   // Deviation 79: events where dchiL != dchiL_L + dchiL_R (should stay 0)
    std::vector<std::array<long, NDETCLASS>> NDetClassTE(GG + 1);
    long nSimTot = 0;
    // ndw MUST start at 0: it doubles as the bound of the per-event buffer clear
    // below, which runs before ndw is reset and therefore reads the PREVIOUS
    // event's value. On the very first event there is no previous value.
    int    gi,       ndw = 0, sq, ndd;
    int    giR,      sqR, nddR; // Roman-side cursor/count, parallel to gi/sq/ndd
    int    flag_det; // nml = 0;
    int    ndw_L, ndw_R;           // per-instrument epoch counts (ndw stays the joint/shared total)
    // Step R1. Per-event, per-survey tallies of epochs at which the two images were both
    // detectable and far enough apart. Reset with ndw_L/ndw_R below -- a counter that leaks
    // across events is the same bug the `ndw` note above guards against.
    long   nres5_L, nres20_L, nresPSF_L, nres5_R, nres20_R, nresPSF_R;
    double dsepMax_L, dsepMax_R;   // largest separation while both were detectable [mas]
    int    flag_det_L, flag_det_R; // per-instrument run-test result (flag_det stays the joint one)
    int    detL, detR, detJ;       // per-instrument / joint detection booleans; FFG[0] = detL or detR or detJ
    int    flagf,  fi; // datf1, datf2;
    double errs,   errg, minc, cade; // fel, , mind
    double fdetRubin, testL, testR; // Step B2: per-survey pre-selection (fdet retired)
    bool   rubinDetectable, romanDetectable, acceptRubin, acceptRoman;
    double mincR, cadeR, errgR; // Roman-side cadence/error tracking, parallel to minc/cade/errg
    double errsR, magnioR; // Roman's per-exposure astrometric error (errRomanA) / noisy magnitude
    double magnio, test, deltaA; // dist,
    double Astar0, As1,  As0;
    double initial;
    double trajm, trajp; //  ddf;
    double chi1,  chi2,   chi3, chi1a, chi2a, chi3a, sil, sil2;
    // Per-instrument duplicates. chi1/chi2/chi3/chi1a/chi2a/chi3a above are the JOINT
    // accumulators (fed by both branches); _L/_R are Rubin-only/Roman-only respectively.
    double chi1_L, chi2_L, chi3_L, chi1a_L, chi2a_L, chi3a_L;
    double chi1_R, chi2_R, chi3_R, chi1a_R, chi2a_R, chi3a_R, silR, sil2R;
    double dchiL, dchiP,  dchiA;
    double dchiL_L, dchiP_L, dchiA_L, dchiL_R, dchiP_R, dchiA_R;
    double flag0, flag1,  flag2;
    double flag0_L, flag1_L, flag2_L, flag0_R, flag1_R, flag2_R;
    double vs1,   vs2,    def1p,  def2p, vsave, dt,    Mpeak;
    double ErtE,  ErpiE,  ErtetE, Erml,  Erdl,  Ermul, Ermus, Eru0, Erfb, nsim;
    double nErAvg; // events the precision means are actually averaged over
    double mbase, fblend, Gamma,  Neven, EFF,   EffiD, EffiL, nerr=0.0;
//    double shib,  Efi;
   
    std::array<int, 3> FFG;
    std::array<double, M> magni, magni0;
    std::array<double, 2> tE, RE, piE, tetE, Vt, u0, Ml, opd, Dl, Ds, vl, vs;
    std::array<double, 2> numd, Struc, murel, vsn, DelT, mbs, fb, fwhm, Map, nbl, Ext;
   
    //fil=fopen("./files/MONTLMC/files/EfLMC2B.dat","r");
    //for(int i=0; i<=GG; ++i){
    //fscanf(fil,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",//18
    //&l.NdtE[i],&l.NstE[i], &l.NdMl[i],&l.NsMl[i], &l.Ndpi[i],&l.Nspi[i],   &l.Ndu0[i], &l.Nsu0[i], &l.Ndmb[i],&l.Nsmb[i], 
    //&l.Ndfb[i],&l.Nsfb[i], &l.Ndmu[i],&l.Nsmu[i], &l.Nhalo[1],&l.Nhalo[0], &l.Nself[1],&l.Nself[0]);}
    //fclose(fil); 
    //cout<<"**** File extinctionf.txt was read ****"<<endl;    
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
// Create/clear files if IMnum == 1
    // Cleared only when the scan STARTS. On a continuation this file, like every other
    // output, holds the earlier chunks' rows and must not be cleared. (cfg.startIndex is
    // read directly here because `resuming` is declared with the other opens below.)

    // File names
//    std::string filnam0 = "./files/MONTLMC/files/BHLSSTMONTS.dat"; // now visible outside the if
    // Named from the population's tag, so a black-hole run cannot append to the bulge run's
    // table. The default population's tag is "5", which is what these files were called
    // before --population existed.
    const std::string tag(gPop->tag);
    std::string fnLDt   = std::string(PATH_OUT_DIR) + "LpLMC"  + tag +  ".dat";
    std::string fnEff   = std::string(PATH_OUT_DIR) + "EfLMC"  + tag +  ".dat";
    std::string fnEffB  = std::string(PATH_OUT_DIR) + "EfLMC"  + tag + "B.dat";
    std::string fnGam   = std::string(PATH_OUT_DIR) + "MapLMC" + tag +  ".dat";
    std::string testf   = "./test"                       + tag +  ".dat";
    std::string fnPair  = "./h3_pair.dat";   //Step H3, written only under --pair-satellite

    // Open files.
    //
    // Every output below ACCUMULATES across the scan: fil2/fil2b write one block per
    // aggregated sightline, fil4/fil5 write per event, fil3 writes one map row per
    // sightline, and the per-event table is appended to per event. So a run that is a
    // CONTINUATION (--start-index > 0) must append to them; truncating would throw away
    // everything the interrupted run wrote.
    //
    // This was not always so, and it destroyed a production table. Before this, all of
    // these except fil3 were opened in the default mode -- ios::out, which truncates --
    // regardless of --start-index. The 2026-09-06 v3 run was paused at scan index 774 and
    // resumed; the resume silently wiped the first 774 sightlines out of test5.dat
    // (1,936,653 rows) and out of EfLMC5/EfLMC5B, and the loss was invisible until the run
    // finished, because the file simply started filling again from the resume point.
    // MapLMC5.dat survived only because fil3 was already ios::app. See DEVIATIONS.md 34.
    const bool resuming = (cfg.startIndex > 0);
    // Deviation 78: a --dry-run must not touch the previous run's outputs, so it opens them all in
    // append mode (and writes nothing); a fresh run truncates EVERY output, MapLMC and LpLMC included
    // (they used to always append, so a re-run in the same directory silently doubled them).
    const bool keepOld = resuming or cfg.dryRun;
    const std::ios::openmode accumulate =
        std::ios::out | (keepOld ? std::ios::app : std::ios::trunc);

    // LpLMC<tag>.dat is an APPEND-MODE OUTPUT -- the per-characterised-event dump inside the
    // sightline loop opens it with ios::app for every row. It was ALSO opened here as an
    // input, and its absence made the startup check below fatal, so a population whose tag
    // had never been run before could not start at all: `--population bh` would abort with
    // "Cannot open one or more files!" before drawing a single star. (The same trap cost a
    // scratch stub run during the Step W3 flush test, where the file had to be copied in by
    // hand.) Create it if it is missing and let the appends do the rest; on a resume it
    // already exists and is left untouched, which is the append-only behaviour recorded in
    // OPEN_ITEMS.md and not changed here.
    { std::ofstream ensureLp(fnLDt, accumulate); }   // create; truncated on a fresh run (Dev. 78)
//    std::ifstream fil1(filnam0);
    std::ofstream fil2(fnEff,   accumulate);
    std::ofstream fil2b(fnEffB, accumulate);
    std::ofstream fil3(fnGam, accumulate);   // was always ios::app (Deviation 78)

    // The per-event table (Step D1). Truncate and write the column header once, in a
    // scope of its own, and leave `filg_in` itself CLOSED.
    //
    // That is not stylistic. The per-event write below calls filg_in.open(); calling
    // open() on an ALREADY-OPEN ofstream sets failbit and does nothing, so the write is
    // silently discarded. Constructing filg_in open (as this used to) therefore threw
    // away the first event of every run -- and only the first, because the matching
    // close() cleared the way for the second open() to succeed. Verified against a
    // three-iteration reproduction, not inferred. The in-memory `records` vector was
    // never affected, so no aggregate statistic was wrong; the flat table just lost a row.
    //
    // The open/append/close per event is deliberate and stays: it flushes each row to
    // disk as it is produced, so a 15-hour run that is interrupted keeps everything it
    // had computed. At ~1.2 s of physics per event the syscalls are not measurable.
    //
    // The header is written when, and only when, the table has no content yet. That single
    // rule covers both uses of --start-index, which the flag itself cannot distinguish:
    //
    //   * CONTINUING an interrupted run, in its directory. The table already holds the
    //     header and every earlier chunk's rows, so it must be appended to and the header
    //     must NOT be rewritten. Opening it here in the default mode -- which truncates --
    //     is what destroyed 1,936,653 rows of the v3 run.
    //   * STARTING a fresh run at an offset, in an empty directory, which is how every A/B
    //     and profiling comparison in this project is set up (see h7_ab.sh, perf_h7.sh).
    //     Here there is nothing to preserve and the header does need writing.
    //
    // An earlier version of this fix made --start-index > 0 mean "continuation" and refused
    // to run against an empty table. That is wrong: it breaks the second case, which is the
    // more common one. Emptiness, not the flag, is the thing worth testing.
    std::streamoff tableBytes = -1;
    {
        std::ifstream probe(testf, std::ios::ate | std::ios::binary);
        tableBytes = probe ? static_cast<std::streamoff>(probe.tellg()) : std::streamoff(-1);
    }
    if (cfg.dryRun) {
        // Deviation 78: leave the previous run's table alone.
    } else if (tableBytes <= 0) {
        std::ofstream head(testf);
        if (!head) {
            std::cerr << "Cannot open " << testf << std::endl;
            return 1;
        }
        head << eventTableHeader() << "\n";
        if (resuming) {
            std::cout << "  NOTE: --start-index " << cfg.startIndex << " with an empty "
                      << testf << ": treating this as a fresh run\n        that begins at "
                      << "an offset, not as a continuation. The table will hold only the "
                      << "sightlines\n        from this index onward." << std::endl;
        }
    } else if (!resuming) {
        // The one remaining way to lose data silently: re-running the scan from the start
        // in a directory that already holds a table. Say so, loudly, before it is gone.
        std::cout << "  NOTE: " << testf << " already held " << tableBytes << " bytes and "
                  << "is being TRUNCATED, because this run\n        starts the scan "
                  << "(--start-index 0). If it was meant to continue one, stop now."
                  << std::endl;
        std::ofstream head(testf);   // truncates
        if (!head) {
            std::cerr << "Cannot open " << testf << std::endl;
            return 1;
        }
        head << eventTableHeader() << "\n";
    } else {
        std::cout << "  Continuing an existing " << testf << " (" << tableBytes
                  << " bytes); rows will be appended." << std::endl;
    }
    // Step H3's side file. Same rule as the event table: the header is written only when
    // there is no content yet, so a continuation appends instead of truncating.
    if (cfg.pairSat) {
        std::ifstream probe(fnPair, std::ios::ate | std::ios::binary);
        const std::streamoff have =
            probe ? static_cast<std::streamoff>(probe.tellg()) : std::streamoff(-1);
        if (have <= 0) {
            std::ofstream h(fnPair);
            h << "# lon lat tE u0 piE tetE du_sat okA_sat okA_nosat okB_sat okB_nosat "
              << "sigtE_sat sigtE_nosat sigpiE_sat sigpiE_nosat "
              << "sigpiER_sat sigpiER_nosat sigtetE_sat sigtetE_nosat "
              << "sigpiEb_sat sigpiEb_nosat relMl_sat relMl_nosat "
              << "condA_sat condA_nosat condB_sat condB_nosat "
            // Ml, Dl, Ds and Vt are appended (not inserted) so the 30-column files written
            // before 2026-09-17 still parse by position. They are here because the pooled
            // event-rate weight needs sqrt(Ml)*Vt*Z(Ds) per event (Deviation 41) and none of
            // it is recoverable from the columns above: theta_E and pi_E give Ml, and
            // theta_E/tE gives mu_rel, but pi_rel = 1/Dl - 1/Ds is one equation in two
            // unknowns, so a paired file without these cannot be weighted at all.
              << "nepL_pk nepR_pk w_area Ml Dl Ds Vt\n";
        }
    }

    std::ofstream filg_in; //opened in append mode per event -- see above

    // Step S1 state. `dumpBuf` holds the CURRENT draw's recorded epochs and is cleared at
    // the top of every draw; `dumpSeq` numbers the events actually written out.
    SampleSpec            dumpSpec;
    std::vector<DumpEpoch> dumpBuf;
    long                  dumpSeq = 0;
    if (!cfg.dumpSpec.empty()) {
        if (!parseSampleSpec(cfg.dumpSpec, dumpSpec)) return 2;
        std::error_code ec;
        std::filesystem::create_directories(dumpSpec.dir, ec);
        if (ec) {
            std::cerr << "ERROR: cannot create sample directory '" << dumpSpec.dir
                      << "': " << ec.message() << "\n";
            return 2;
        }
        std::cout << "  Sample dump ON -> " << dumpSpec.dir << "/  (";
        for (std::size_t i = 0; i < dumpSpec.classes.size(); ++i)
            std::cout << (i ? ", " : "") << dumpSpec.classes[i].name << " x"
                      << dumpSpec.classes[i].quota;
        std::cout << ")" << std::endl;
    }

    // Check all
    if (!fil2 || !fil2b || !fil3) {
        std::cerr << "Cannot open one or more files!" << std::endl;
        return 1;
    }

    save = 0;

    // ----------------------------------------------------------------------
    // Roman field placements and the coverage guard.
    //
    // RomanBaseline.dat repeats a small set of distinct pointings once per visit: the six
    // GBTDS fields at each of the two rolls, i.e. 12 (centre, layout) placements. Collect them,
    // then verify the sightline grid puts at least one sightline ON A DETECTOR of each. A grid
    // that misses a placement produces a run with no Roman epochs there, which does not crash
    // and does not warn: it quietly reports Rubin-only results in joint-labelled columns.
    // ----------------------------------------------------------------------
    std::vector<FieldPlacement> romanFields;
    for (int i = 0; i < NlRoman; ++i) {
        bool seen = false;
        for (const auto& f : romanFields)
            if (std::fabs(f.l - ro->l[i]) < 1e-6 and std::fabs(f.b - ro->b[i]) < 1e-6
                and f.layout == ro->layout[i]) {
                seen = true; break;
            }
        if (!seen) romanFields.push_back({ro->l[i], ro->b[i], ro->layout[i]});
    }

    // Scan bounds (Deviation 69). The region is every point within scanReach of a field centre
    // (see SCAN_RUBIN_REACH in config/parameters.h); the grid's bounding box is that, with its origin on a
    // multiple of the coarse step so sightlines sit at round coordinates. The stub is a 0.1 x
    // 0.1 deg test patch inside field 3, which both rolls image (spring centre l 0.500, autumn
    // 0.350; it was l 0.5-0.6, b -1.0..-0.9 before, which the adopted fields do not reach).
    const double scanReach = SCAN_RUBIN_REACH + gl.rField;
    double fLmin = 1e9, fLmax = -1e9, fBmin = 1e9, fBmax = -1e9;
    for (const auto& f : romanFields) {
        fLmin = std::min(fLmin, f.l); fLmax = std::max(fLmax, f.l);
        fBmin = std::min(fBmin, f.b); fBmax = std::max(fBmax, f.b);
    }
    const double lonMin = cfg.stubPatch ?  0.40  : gridStep * std::floor((fLmin - scanReach) / gridStep);
    const double lonMax = cfg.stubPatch ?  0.50  : fLmax + scanReach;
    const double latMin = cfg.stubPatch ? -1.45  : gridStep * std::floor((fBmin - scanReach) / gridStep);
    const double latMax = cfg.stubPatch ? -1.35  : fBmax + scanReach;

    // Inside the scan region: within scanReach of any field centre (the stub keeps its square).
    auto inScan = [&](double lon, double lat) {
        if (cfg.stubPatch) return true;
        for (const auto& f : romanFields)
            if (std::hypot(lon - f.l, lat - f.b) <= scanReach) return true;
        return false;
    };

    // Every Rubin pointing in the visit list must be one that can image the scan region, i.e.
    // centred within scanReach + FoV of a field centre. A pointing farther out means the list
    // was extracted for a different region (readbaselineBulge.py uses the same rule).
    if (not cfg.stubPatch) {
        int nFar = 0;
        for (int i = 0; i < Nl; ++i) {
            bool near = false;
            for (const auto& f : romanFields)
                if (std::hypot(ls->l[i] - f.l, ls->b[i] - f.b) <= scanReach + RUBIN_MAX_RADIUS + 1e-3) { near = true; break; }
            if (!near) nFar += 1;
        }
        if (nFar > 0) {
            std::cerr << "ERROR: " << nFar << " of " << Nl << " Rubin visits in BulgeBaseline.dat "
                      << "are centred farther than " << scanReach + RUBIN_MAX_RADIUS << " deg from every "
                      << "Roman field, so the visit list was built for another region. "
                      << "Regenerate it with Baseline/readbaselineBulge.py.\n";
            return 1;
        }
    }

    // ----------------------------------------------------------------------
    // The sightline list, with the sky area each sightline stands for (Step E1).
    //
    // WHY THIS IS NOT JUST A NESTED LOOP ANY MORE. The scan region is 68 deg^2; Roman's six
    // GBTDS fields cover ~1.7 deg^2 of it, about 2.6%. A uniform grid therefore spends 97% of
    // its wall clock on sightlines Roman never visits, where the joint Fisher matrix IS
    // Rubin's matrix and nothing whatever can be learned about the combination of the two
    // surveys. The 2026-08-30 production run bought 74,812 joint detections and only 1,950 of
    // them -- 2.6%, exactly the area fraction -- inside the footprint. Every result that is
    // currently sample-limited (F3 panel (a), the F4 precision fractions, the long-tE piE
    // null) is limited by that same 1,950.
    //
    // So the grid is stratified. Two strata:
    //   R  sightlines within FoVRoman of a GBTDS field centre, on the FINE grid (fineStep)
    //   O  everything else, on the COARSE grid (gridStep), one sightline per coarse cell
    // and each carries `area`, the deg^2 of sky it represents. An R sightline carries one
    // fine cell; an O sightline carries however many fine cells of its coarse block fall
    // outside the footprint -- NOT the whole coarse cell, or a block straddling the footprint
    // edge would count its overlapping part twice, once in each stratum.
    //
    // THE ONE INVARIANT: sum(area) over the list equals the scanned area. It is asserted
    // below rather than trusted, because every absolute yield in deg^-2 downstream is that
    // sum in disguise, and an area bookkeeping error does not look like an error -- it looks
    // like a survey that found more events than it did.
    //
    // WHAT STRATIFICATION DOES AND DOES NOT BIAS. Nothing computed *at* a sightline changes:
    // detection efficiency, per-event precision, and every ratio conditional on the sample
    // are untouched, because which sightlines were visited is not an input to any of them.
    // What changes is any quantity POOLED ACROSS sightlines -- a survey-wide yield, a
    // histogram over all events, the "all joint detections" panel of F4. Those must weight
    // each event by its sightline's `w_area` or they will describe a sky that is 2.6% Roman
    // by area and (say) 40% Roman by sample. The weight is written into every row of the
    // per-event table for exactly this reason.
    //
    // With --stride-roman absent, kSub == 1: the fine grid IS the coarse grid, every block is
    // a single cell that is its own representative, every area is gridStep^2, and the list is
    // the same points in the same order as the old nested loop -- so the RNG stream, and
    // therefore the run, is bit-identical to before this step.
    // ----------------------------------------------------------------------
    // The fine stratum (Deviation 69): every fine cell that overlaps the detector outline's
    // bounding rectangle of any placement. A grid point stands for the cell extending from it
    // in +l and +b, so a cell overlaps when its point lies within one fine step below the
    // rectangle; the same margin is kept on the other side. Chip-gap cells are therefore in the
    // fine stratum too: which sightlines actually see Roman is decided by the detector test in
    // matchVisibleEpochs, not here.
    auto insideFootprint = [&](double lon, double lat) {
        for (const auto& f : romanFields) {
            const double dl = lon - f.l, db = lat - f.b;
            const int k = f.layout;
            if (dl >= gl.dlMin[k] - fineStep and dl <= gl.dlMax[k] + fineStep
                and db >= gl.dbMin[k] - fineStep and db <= gl.dbMax[k] + fineStep) return true;
        }
        return false;
    };

    const int nLonGrid = int(std::floor((lonMax - lonMin) / gridStep + 1e-9)) + 1;
    const int nLatGrid = int(std::floor((latMax - latMin) / gridStep + 1e-9)) + 1;
    // Each coarse cell is subdivided into exactly kSub x kSub fine cells, so the fine grid has
    // nLonGrid*kSub columns -- not (nLonGrid-1)*kSub+1. The difference is the edge convention:
    // a grid POINT stands for the cell extending from it, which is what makes the areas tile
    // exactly and the sum below come out on the nose.
    const int    nLonFine = nLonGrid * kSub;
    const int    nLatFine = nLatGrid * kSub;
    const double cellArea = (gridStep * gridStep) / double(kSub * kSub); // deg^2, one fine cell

    struct Sightline { double lon, lat, area; bool inFootprint; int col; int romanClass; };

    // Roman coverage class of a sky point (Deviation 69): bit 0 = on a detector in the spring
    // roll, bit 1 = in the autumn roll. 0 none, 1 spring only, 2 autumn only, 3 both.
    auto romanClassAt = [&](double lon, double lat) {
        int c = 0;
        for (const auto& p : romanFields)
            if (inDetector(gl, p.layout, lon - p.l, lat - p.b)) c |= (1 << p.layout);
        return c;
    };

    // Pass 1: for each coarse block, how many of its fine cells survive the corner cut and lie
    // OUTSIDE the footprint, and which of them represents that area.
    const size_t nBlocks = size_t(nLonGrid) * size_t(nLatGrid);
    std::vector<int>  blockOutCount(nBlocks, 0);
    std::vector<long> blockRep(nBlocks, -1);      // fine cell index (iF*nLatFine + jF)
    long nFineKept = 0;
    for (int iF = 0; iF < nLonFine; ++iF) {
        const double lon = lonMin + iF * fineStep;
        for (int jF = 0; jF < nLatFine; ++jF) {
            const double lat = latMin + jF * fineStep;
            if (!inScan(lon, lat)) continue;       // same region test as the scan below
            nFineKept += 1;
            if (insideFootprint(lon, lat)) continue;
            const size_t b = size_t(iF / kSub) * size_t(nLatGrid) + size_t(jF / kSub);
            if (blockRep[b] < 0) blockRep[b] = long(iF) * nLatFine + jF;
            blockOutCount[b] += 1;
        }
    }

    // Pass 2: the list itself, in the same iLon-major order the old loop walked.
    std::vector<Sightline> scan;
    scan.reserve(size_t(nFineKept));
    std::vector<int> fieldHits(romanFields.size(), 0);
    long nSightlines = 0, nSightlinesRoman = 0;
    double areaFootprint = 0.0, areaOutside = 0.0;
    // Point-sampled sky area on a detector, per roll: the grid's estimate of what Roman images
    // in a spring / an autumn season, to compare with the exact detector area.
    std::array<double, GBTDS_NLAYOUT> areaOnDetector{};
    for (int iF = 0; iF < nLonFine; ++iF) {
        const double lon = lonMin + iF * fineStep;
        for (int jF = 0; jF < nLatFine; ++jF) {
            const double lat = latMin + jF * fineStep;
            if (!inScan(lon, lat)) continue;
            const bool anyField = insideFootprint(lon, lat);
            if (anyField) {
                std::array<bool, GBTDS_NLAYOUT> onDet{};
                for (size_t f = 0; f < romanFields.size(); ++f) {
                    const auto& p = romanFields[f];
                    if (inDetector(gl, p.layout, lon - p.l, lat - p.b)) {
                        fieldHits[f] += 1; onDet[p.layout] = true;
                    }
                }
                for (int k = 0; k < GBTDS_NLAYOUT; ++k) if (onDet[k]) areaOnDetector[k] += cellArea;
                scan.push_back({lon, lat, cellArea, true, iF, int(onDet[0]) | (int(onDet[1]) << 1)});
                areaFootprint    += cellArea;
                nSightlinesRoman += 1;
            } else {
                const size_t b = size_t(iF / kSub) * size_t(nLatGrid) + size_t(jF / kSub);
                if (blockRep[b] != long(iF) * nLatFine + jF) continue;  // not this block's rep
                const double a = blockOutCount[b] * cellArea;
                scan.push_back({lon, lat, a, false, iF, 0});
                areaOutside += a;
            }
            nSightlines += 1;
        }
    }

    // ----------------------------------------------------------------------
    // Post-stratifying the footprint's area weights (Deviation 69).
    //
    // A footprint sightline is a POINT: it either lands on a detector in a given roll or in a
    // chip gap. The detectors are 0.125 deg across with gaps of 0.008-0.026 deg, so a grid of
    // 0.1 deg aliases against them: at --stride-roman 5 the grid put 1.46 deg^2 on a detector per
    // roll where the exact figure is 1.68 -- every Roman yield per deg^2 13% low -- and even a
    // 0.02 deg grid lands 2% high. Rather than chase that with CPU, each footprint sightline's
    // area is rescaled so that each coverage CLASS (none / spring only / autumn only / both)
    // carries its exact sky area inside the footprint stratum:
    //     area_i <- cellArea * exact(c_i) / grid(c_i)
    // exact(c) from sub-sampling every fine cell on a ~0.002 deg raster. The classes partition
    // the stratum, so the stratum's total area -- and the invariant below -- is unchanged. The
    // sightlines of a class still sample its sky at the grid points; only how much sky each
    // stands for changes.
    // ----------------------------------------------------------------------
    std::array<double, 4> classGrid{}, classExact{};
    std::array<long, 4>   classN{};
    for (const auto& sl : scan)
        if (sl.inFootprint) { classGrid[sl.romanClass] += cellArea; classN[sl.romanClass] += 1; }
    {
        const int nSub = std::max(1, int(std::ceil(fineStep / 0.002 - 1e-9)));
        const double h = fineStep / nSub;
        for (const auto& sl : scan) {
            if (!sl.inFootprint) continue;
            for (int a = 0; a < nSub; ++a)
                for (int c = 0; c < nSub; ++c)
                    classExact[romanClassAt(sl.lon + (a + 0.5) * h, sl.lat + (c + 0.5) * h)] += h * h;
        }
    }
    for (int c = 0; c < 4; ++c) {
        if (classExact[c] > 0.0 and classN[c] == 0) {
            std::cerr << "ERROR: Roman coverage class " << c << " has " << classExact[c]
                      << " deg^2 of sky in the footprint but no sightline samples it. Use a finer "
                      << "--stride-roman.\n";
            return 1;
        }
    }
    areaFootprint = 0.0;
    for (auto& sl : scan)
        if (sl.inFootprint) {
            sl.area = cellArea * classExact[sl.romanClass] / classGrid[sl.romanClass];
            areaFootprint += sl.area;
        }
    std::cout << "Roman coverage classes in the footprint stratum (none / spring / autumn / both):"
              << "\n  sightlines  " << classN[0] << " / " << classN[1] << " / " << classN[2]
              << " / " << classN[3]
              << "\n  grid area   " << classGrid[0] << " / " << classGrid[1] << " / "
              << classGrid[2] << " / " << classGrid[3] << " deg^2"
              << "\n  exact area  " << classExact[0] << " / " << classExact[1] << " / "
              << classExact[2] << " / " << classExact[3] << " deg^2 (weights rescaled to these)"
              << std::endl;

    // The invariant. Every fine cell inside the scan region is represented exactly once,
    // either by itself (footprint) or by its block's representative (outside).
    const double areaScanned = areaFootprint + areaOutside;
    {
        const double areaExpect = nFineKept * cellArea;
        if (std::fabs(areaScanned - areaExpect) > 1e-9 * std::max(1.0, areaExpect)) {
            std::cerr << "ERROR: sightline area weights sum to " << areaScanned
                      << " deg^2 but the scanned grid is " << areaExpect << " deg^2. "
                      << "Every absolute yield downstream is this sum; refusing to run.\n";
            return 1;
        }
        std::cout << "Sightline grid: " << scan.size() << " sightlines ("
                  << nSightlinesRoman << " inside Roman's footprint at " << fineStep
                  << " deg, " << (scan.size() - size_t(nSightlinesRoman)) << " outside at "
                  << gridStep << " deg), covering " << areaScanned << " deg^2 ("
                  << areaFootprint << " footprint + " << areaOutside << " outside)." << std::endl;
        // The grid samples the detector mosaic at points, so its on-detector area is an
        // estimate; the exact value is 6 fields x 18 detectors. Their difference is the
        // sampling error of every per-deg^2 Roman yield from this grid.
        for (int k = 0; k < GBTDS_NLAYOUT; ++k) {
            double exact = 0.0;
            for (const ScaRect& r : gl.sca[k]) exact += (r.l1 - r.l0) * (r.b1 - r.b0);
            int nPlaced = 0;
            for (const auto& f : romanFields) nPlaced += (f.layout == k);
            exact *= nPlaced;
            std::cout << "Roman detector area, " << (k == 0 ? "spring" : "autumn") << " roll: "
                      << areaOnDetector[k] << " deg^2 on the grid vs " << exact
                      << " deg^2 exact (" << nPlaced << " fields)." << std::endl;
        }
    }

    const size_t nFieldsCovered = std::count_if(fieldHits.begin(), fieldHits.end(),
                                                [](int h){ return h > 0; });
    // Fatal for a full-region scan, advisory for --stub: the stub patch is 0.1x0.1 deg
    // by design and cannot possibly reach all twelve placements.
    if (nFieldsCovered < romanFields.size() and not cfg.stubPatch) {
        std::cerr << "ERROR: the sightline grid (stride " << cfg.stride << ", step "
                  << gridStep << " deg) misses " << (romanFields.size() - nFieldsCovered)
                  << " of " << romanFields.size() << " Roman field placements:\n";
        for (size_t f = 0; f < romanFields.size(); ++f)
            if (fieldHits[f] == 0)
                std::cerr << "    field at (l,b) = (" << romanFields[f].l << ", "
                          << romanFields[f].b << "), " << (romanFields[f].layout ? "autumn" : "spring")
                          << " roll, has no sightline on a detector\n";
        std::cerr << "Those fields would contribute no Roman epochs, and the run would "
                     "report joint columns built from Rubin data alone. Use a smaller "
                     "--stride.\n";
        return 1;
    }
    if (nFieldsCovered < romanFields.size() and cfg.stubPatch) {
        std::cout << "NOTE: --stub reaches " << nFieldsCovered << " of "
                  << romanFields.size() << " Roman field placements. Expected for a patch this "
                  << "small; joint statistics from it describe those fields only.\n";
    }

    // ----------------------------------------------------------------------
    // Provenance. Written before any science output, so an interrupted run still
    // records what it was. Reproducibility here is not optional: an earlier run
    // (commit e47390a) produced detection statistics that turned out to describe a
    // truncated visit stream rather than the survey, and nothing in its output
    // said which model had produced it.
    //
    // The sky-area weight is the piece that is easy to lose. Neven is a surface
    // density in deg^-2, and nothing in this program sums sightlines into a
    // survey-wide yield -- that happens downstream, where each sightline must be
    // weighted by the sky area it stands for. At stride N that area is
    // (N*dd)^2, NOT dd^2, so a strided run aggregated as if it were unstrided is
    // wrong by a factor of N^2 (100x at the default stride of 10).
    //
    // Under Step E1's stratification that weight is no longer one number for the
    // whole run: footprint sightlines stand for a fine cell and outside ones for
    // (up to) a coarse cell. area_per_sightline below is therefore only meaningful
    // when stratified=0, and the authoritative weight is the per-row `w_area`
    // column of the event table (and the last column of the map file). The header
    // says so, so that a downstream script cannot quietly use the wrong one.
    // ----------------------------------------------------------------------
#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif
    const double areaPerSightline = gridStep * gridStep; // deg^2; unstratified runs only

    {
        std::ostringstream prov;
        prov << "# Roman+Rubin microlensing forecast -- run provenance\n"
             << "# git_commit          " << GIT_COMMIT << "\n"
             << "# built               " << __DATE__ << " " << __TIME__ << "\n"
             // Which lens population produced this table. The analysis layer reads it: the
             // pooled event-rate weight carries a sqrt(Ml) factor that is only correct for
             // the mass function actually sampled, so a figure must know which one that was.
             << "# population          " << gPop->name << "   # " << gPop->note << "\n"
             << "# population_tag      " << gPop->tag << "\n"
             << "# lens_mass_range     " << gPop->mlMin << " " << gPop->mlMax << "   # Msun\n"
             << "# lens_mass_grid      " << (gPop->logGrid ? "log" : "linear") << "\n"
             << "# stride              " << cfg.stride << "\n"
             << "# grid_step_deg       " << gridStep << "\n"
             << "# stratified          " << (kSub > 1 ? 1 : 0)
             << "   # 1 = per-sightline areas differ; use the w_area column, not the scalar below\n"
             << "# area_per_sightline  " << areaPerSightline
             << "   # deg^2 -- weight for Neven; VALID ONLY IF stratified = 0\n"
             << "# stride_roman        " << cfg.strideRoman << "\n"
             << "# fine_step_deg       " << fineStep << "\n"
             << "# k_subdivision       " << kSub << "   # fine cells per coarse cell, per axis\n"
             << "# area_scanned_deg2   " << areaScanned << "\n"
             << "# area_footprint_deg2 " << areaFootprint << "\n"
             << "# area_outside_deg2   " << areaOutside << "\n"
             << "# lon_range_deg       " << lonMin << " " << lonMax << "\n"
             << "# lat_range_deg       " << latMin << " " << latMax << "\n"
             << "# scan_region         within " << scanReach << " deg of a Roman field centre"
             << "   # 2*FoV + field reach " << gl.rField << "\n"
             << "# n_sightlines        " << nSightlines << "\n"
             << "# n_sightlines_roman  " << nSightlinesRoman << "\n"
             << "# roman_fields        " << nFieldsCovered << " of " << romanFields.size()
             << " placements covered\n"
             << "# roman_layout        Baseline/gbtds_layout (gbtds_{spring,autumn}_2026.4.3, "
             << GBTDS_NSCA << " SCAs, side " << gl.scaSide << " deg)\n"
             << "# roman_area_grid     " << areaOnDetector[0] << " " << areaOnDetector[1]
             << "   # deg^2 on a detector, spring autumn, as point-sampled by the grid\n"
             << "# roman_class_exact   " << classExact[0] << " " << classExact[1] << " "
             << classExact[2] << " " << classExact[3]
             << "   # deg^2 none/spring/autumn/both; footprint w_area post-stratified to these\n"
             << "# roman_class_grid    " << classGrid[0] << " " << classGrid[1] << " "
             << classGrid[2] << " " << classGrid[3] << "\n"
             << "# events_target       " << cfg.iconTarget << "   # icon\n"
             << "# lenses_target       " << cfg.nlensTarget << "   # nlens\n"
             << "# nerr_target         " << cfg.nerrTarget << "\n"
             << "# maxdraws            " << cfg.maxDraws << "   # per-sightline draw cap\n"
             << "# seed                " << cfg.seedBase
             << "   # base; each sightline re-seeded from (seed, index) (Deviation 77)\n"
             << "# end_index           " << cfg.endIndex << "   # -1 = to the end\n"
             << "# start_index         " << cfg.startIndex
             << "   # sightlines skipped; >0 means this run RESUMES an earlier one\n"
             << "# region              " << (cfg.stubPatch ? "stub patch" : "full") << "\n"
             << "# Tobs_days           " << Tobs << "\n"
             << "# roman_seasons       " << sched.seasons.size() << "\n"
             << "# roman_mission_days  " << sched.missionStart << " " << sched.missionEnd << "\n"
             << "# season_gap_thresh   " << SEASON_GAP_MIN_DAYS
             << "   # d; max in-season spacing " << sched.maxInSeasonSpacing
             << ", min inter-season gap " << sched.minSeasonGap
             << ", shortest season " << sched.minSeasonLength << "\n"
             << "# Nl                  " << Nl << "\n"
             << "# NlRoman             " << NlRoman << "\n"
             << "# FoV_rubin_deg       " << FoV << "\n"
             << "# roman_noise         ast: errRomanA(m_AB - " << F146_AB_MINUS_VEGA
             << ") [Vega anchors]; phot: Penny+2019 curve anchored to 5-sigma "
             << ROMAN_DEPTH5_AB << " AB (66 s); depth " << thre[6] << ", saturation "
             << satu[6] << " AB   # Deviation 72\n"
             << "# extinction          files/ext/ext_tables.dat: " << ex->nTables << " x "
             << ex->nDist << ", k " << ex->k << " --" << ex->built << "\n"

             << "# rng_seed            " << cfg.seedBase << "\n"
             // Step H1. 1 = Roman at Sun-Earth L2 (physical); 0 = Roman at the centre of the
             // Earth, which is what every run before H1 did. Any piE forecast from a run with
             // 0 here contains only the annual Earth-orbit parallax.
             << "# satellite_parallax  " << (cfg.noSatPar ? 0 : 1)
             << "   # 0 = Roman forced to Earth's position\n"
             << "# L2_offset_AU        " << (cfg.noSatPar ? 0.0 : L2_OFFSET_AU) << "\n"
             << "# dump_samples        " << (cfg.dumpSpec.empty() ? std::string("none")
                                                                    : cfg.dumpSpec) << "\n"
             << "# pair_satellite      " << (cfg.pairSat ? 1 : 0)
             << "   # Step H3: every detection characterised at L2 AND at Earth\n"
             << "# dchi_det            " << cfg.dchiDet
             << "   # Step H7 fixed detection bar; every yield is conditioned on it\n";
        if (!cfg.dryRun) {   // Deviation 78: a dry run must not overwrite the last run's provenance
            std::ofstream fprov(std::string(PATH_OUT_DIR) + "run_provenance.txt");
            if (!fprov) {
                std::cerr << "Cannot write run_provenance.txt\n";
                return 1;
            }
            fprov << prov.str();
        }
        std::cout << prov.str() << std::flush;
    }

    if (cfg.dryRun) {
        // Per-stratum sightline counts, and what they cost. The footprint sightlines are the
        // ones worth buying and the ones that are expensive; printing both together is what
        // makes --stride-roman a decision rather than a guess.
        const long nOutside = long(scan.size()) - nSightlinesRoman;
        std::cout << "\n---- dry run: sightline strata ----\n"
                  << "  footprint  " << nSightlinesRoman << " sightlines at " << fineStep
                  << " deg, " << areaFootprint << " deg^2 ("
                  << (areaFootprint > 0.0 ? cellArea : 0.0) << " deg^2 each)\n"
                  << "  outside    " << nOutside << " sightlines at " << gridStep
                  << " deg, " << areaOutside << " deg^2\n"
                  << "  total      " << scan.size() << " sightlines, " << areaScanned
                  << " deg^2\n"
                  << "  footprint share: " << (100.0 * double(nSightlinesRoman) / double(scan.size()))
                  << "% of sightlines, " << (100.0 * areaFootprint / areaScanned)
                  << "% of area\n"
                  << "No stars drawn. Remove --dry-run to run the simulation.\n";
        return 0;
    }

///HHHHHHHHHHHHHHHHHHHHH Monte Carlo Simulation HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    // Index-driven rather than accumulate-and-compare. The old loop did
    // `for (lon = 0.5; lon <= 0.6; lon += dd)` and lost its last row AND last column
    // to floating-point drift -- 0.02 is not representable in binary, so after five
    // additions the accumulator sits a few ulp above the bound and the comparison
    // fails. That scan advertised 36 sightlines and ran 25. Computing each position
    // as lonMin + i*gridStep keeps the count exact and makes it predictable ahead of
    // the run, which the provenance block and the coverage guard both rely on.

    // Sightline bookkeeping. The scan region reaches 2*FoV + the field reach past the Roman
    // field centres (Deviation 69), so it deliberately extends to sky only a Roman-overlapping
    // Rubin pointing's far edge can catch; part of the grid may have no coverage at all.
    // These counters are what make that visible: a density in deg^-2 computed downstream
    // must know how much of the scanned area yielded nothing, and why.
    int nSkipNoCoverage = 0; //no Rubin AND no Roman epochs -- never entered the star loop
    int nSkipBarren     = 0; //had epochs, but produced no characterised event to aggregate
    int nCapped         = 0; //stopped by --maxdraws with its budget unmet
    int nAggregated     = 0; //reached the per-sightline aggregation block
    // Areas rather than counts, because under stratification the sightlines no longer
    // stand for equal pieces of sky and a count is not an area (Step E1).
    double areaAggregated = 0.0, areaNoCoverage = 0.0, areaBarren = 0.0;
    int lastCol = -1;
    long iScan = -1;
    for (const auto& sightline : scan) {
        iScan += 1;
        s->lon = sightline.lon;
        s->lat = sightline.lat;
        // The sky area THIS sightline stands for. Written into every event row it produces;
        // any statistic pooled across sightlines has to weight by it (see the list build).
        const double wArea = sightline.area;
        // nri/nde keep their old meaning -- longitude column index, and index within that
        // column -- which is what the map file and the event table record. Assigned rather
        // than incremented: the old loop bumped nri once per iLon whether or not any
        // sightline in that column survived the corner cut, so an incrementing counter here
        // would silently renumber the columns the moment a fully-cut column existed. Under
        // stratification the column index is a fine-grid one, so it steps by kSub between
        // coarse columns.
        nri = sightline.col;
        if (sightline.col != lastCol) { nde = -1; lastCol = sightline.col; }
        nde += 1;

        // --start-index. Resuming an interrupted run. The nri/nde bookkeeping above runs for
        // skipped sightlines too, deliberately: nde counts position within a longitude column
        // and is written into the map file, so skipping it would renumber the first partial
        // column and a resumed run would not concatenate onto the interrupted one. Only the
        // simulation is skipped, not the numbering.
        //
        // The index to resume AT is the number of sightlines entered, which is the number of
        // "NEW STEP" lines in the interrupted run's log -- this print sits above the
        // no-coverage and barren skips, so it counts every sightline the scan reached. The map
        // file counts only sightlines that produced something to aggregate and is therefore a
        // smaller number; resuming at it would redo the difference and duplicate rows.
        //
        // This is safe precisely because `scan` is built deterministically from --stride,
        // --stride-roman and --stub and does not depend on the RNG. Note that the random
        // sequence IS advanced by the sightlines a full run would have simulated, so a
        // resumed run is not bit-identical to an uninterrupted one -- it is a valid
        // continuation with a different draw sequence, not a replay.
        if (iScan < cfg.startIndex) continue;
        if (cfg.endIndex >= 0 and iScan >= cfg.endIndex) break;   // Deviation 77
        // Deviation 77: this sightline's own random stream (see sightlineSeed in Bulge.h). The note
        // above, that a resumed run is "not a replay", no longer applies: it IS one.
        rng.seed(sightlineSeed(cfg.seedBase, iScan));
        {
            cout << ">>>>>>>>>>> NEW STEP " << nde << " <<<<<<<<\t nri:  " << nri << endl;
            cout << "longtitude: " << s->lon << "\t latitude: " << s->lat << endl;

            cade = 0.0;

            // Deviation 80: on LSSTCam's active silicon for this visit's pointing and rotation
            // (was: within a 1.75-deg circle). A cheap (l, b) distance cut first.
            double slRA, slDec;
            galToIcrs(s->lon, s->lat, slRA, slDec);
            auto rubinCovers = [&](int i) {
                const double dl = s->lon - ls->l[i], db = s->lat - ls->b[i];
                if (dl * dl + db * db > (RUBIN_MAX_RADIUS + 0.1) * (RUBIN_MAX_RADIUS + 0.1)) return false;
                return onLsstCam(slRA, slDec, ls->RA[i], ls->DEC[i], ls->rot[i]);
            };
            auto romanCovers = [&](int i) {
                return inDetector(gl, ro->layout[i], s->lon - ro->l[i], s->lat - ro->b[i]);
            };
            ndd  = matchVisibleEpochs("LSST",  rubinCovers, ls->tim, Nl,      ls->ct, minc);
            nddR = matchVisibleEpochs("Roman", romanCovers, ro->tim, NlRoman, ro->ct, mincR);

            // Rubin's depth per band at this sightline (Deviation 73): the median of its matched
            // visits' own 5-sigma depths, for the pre-selection below. A band with no visit here
            // gets -inf, so it cannot count toward "detectable in >= 2 bands".
            std::array<double, 6> rubinDepthMed;
            {
                std::array<std::vector<double>, 6> d5;
                for (int k = 0; k < ndd; ++k) {
                    const int v = int(ls->ct[k]);
                    d5[ls->filter[v]].push_back(ls->sig5[v]);
                }
                for (int b = 0; b < 6; ++b) {
                    if (d5[b].empty()) { rubinDepthMed[b] = -std::numeric_limits<double>::infinity(); continue; }
                    std::nth_element(d5[b].begin(), d5[b].begin() + d5[b].size() / 2, d5[b].end());
                    rubinDepthMed[b] = d5[b][d5[b].size() / 2];
                }
            }
            cout << "ndd (LSST): "  << ndd  << "\t minc (LSST): "  << minc  << endl;
            cout << "ndd (Roman): " << nddR << "\t minc (Roman): " << mincR << endl;

            // Empty sky. Neither survey visits this sightline, so no light curve can ever
            // have a datum, no event can be detected, and nerr can never advance -- the
            // do/while below would spin forever. Skipping is not an approximation: a
            // sightline with no epochs contributes exactly zero events to the yield.
            //
            // It has to be a `continue` rather than a run that finds nothing, because the
            // aggregation block at the end of this loop asserts CHECK(numd[1] != 0.0) and
            // CHECK(nerr != 0.0) -- reaching it with an empty field aborts the whole run.
            //
            // The skipped area is NOT written to the map file, so anything converting the
            // per-sightline Neven density into a total count must use the aggregated
            // sightline count reported at the end of the run, not the grid size.
            if (ndd == 0 and nddR == 0) {
                nSkipNoCoverage += 1;
                areaNoCoverage  += wArea;
                cout << "  SKIP: no Rubin and no Roman epochs at this sightline" << endl;
                continue;
            }
 
            icon  = 0;
            nlens = 0;
            nDetClass.fill(0);
            nsim  = 0.0;
            nerr  = 0.0;
            for (int i = 0; i < Num; ++i) { s->nssim[i] = 0.0;  s->nsdet[i] = 0.0; }
            for (int i = 0; i <= GG; ++i) { l->nstE[i]  = 0.0;  l->ndtE[i]  = 0.0; }

            s->TET = (360.0 - s->lon) / RAa;///radian s.lon/RA;//
            s->FI  = s->lat / RAa;
            
            Disk_model(*s, 1);
            int sightlineIdx = nearestSightline(*ex, s->lon, s->lat);

            records.clear();
            do { //Start of visible star
                nsim += 1.0;
                func_source(*s, *cm, *ex, sightlineIdx);
                func_lens(*l, *s, *ex, sightlineIdx);
//                std::cerr << "nsim=" << nsim << "  Ds=" << s->Ds << "  mass=" << s->mass
//                          << "  nums=" << s->nums << "  Ml=" << l->Ml << "  u0=" << l->u0 << "\n";
                optical_depth(*s);

                // tE-histogram bin for THIS event, computed for every draw rather than only
                // for detected ones: nstE is the denominator of the detection efficiency, so
                // it has to count everything simulated. Previously gg was only evaluated
                // inside the detection branch and neither counter was ever incremented, so
                // ndtE stayed identically zero and EFF, Gamma and Neven with it (the run then
                // aborted on CHECK(EFF > 0.0) as soon as a field managed to complete).
                gg = FunctE(*l);
                l->nstE[gg] += 1.0;
                nSimTot += 1;

                // The other six efficiency axes, on the same principle as tE: the DENOMINATOR
                // has to count everything simulated, so the bin is computed here, before the
                // detection test, and the numerator is incremented in the detection branch
                // using these same indices. Every input is already set for this draw --
                // func_source filled s->Map and s->blend, func_lens filled u0, pirel and
                // murel -- so this is six array lookups and no new physics.
                //
                // Unlike the tE pair there is no per-sightline lowercase pair for these; the
                // N* arrays accumulate over the whole run, which is what the EfLMC writer
                // reports and why nothing resets them per sightline.
                ss = FuncMl(*l);                    // lens mass
                qq = FuncPi(*l);                    // log10 relative parallax
                ww = Funcu0(*l);                    // impact parameter
                vv = FuncMu(*l);                    // relative proper motion
                zz = FuncMb(*l, s->Map[2]);         // source baseline magnitude, r band
                pp = FuncFb(*l, s->blend[2]);       // blend fraction, r band
                l->NsMl[ss] += 1.0;
                l->Nspi[qq] += 1.0;
                l->Nsu0[ww] += 1.0;
                l->Nsmu[vv] += 1.0;
                l->Nsmb[zz] += 1.0;
                l->Nsfb[pp] += 1.0;

                s->nssim[s->nums] += 1.0;
                flagf   = 0;
                dumpBuf.clear(); //Step S1: this draw's epoch buffer. See the note on `ndw`.
                dclsEvent = DET_NONE; //DetClass for this draw; stays NONE if no light curve
                initial = 0.0;
                // (Step B2: the old single `test = RandR(0.0,1.0)` draw consumed here by
                // `test <= s->blend[2]` is gone — testL/testR are now drawn fresh right
                // before the per-survey pre-selection check, below.)

                // Clear only the prefix the PREVIOUS event dirtied -- `ndw` is not reset
                // until a few lines below, so it still holds that count here. Clearing all
                // `coun` slots (Nl + NlRoman = 306,092, times seven arrays = 2.1M writes)
                // to reset the ~2,000 an event actually uses was ~150x of pure waste per
                // draw, and got 16x more expensive when NlRoman became the real visit count.
                //
                // Safe because every slot in [0, ndw) is fully written before it is read --
                // both fill branches write all seven arrays at index ndw before incrementing
                // it -- and FisherM reads only [0, ndw). Slots past the previous ndw are
                // therefore untouched since construction, i.e. already zero.
                //
                // Kept as a prefix clear rather than deleted outright so the clean-slate
                // invariant survives: if a future edit ever advances ndw without filling
                // every array, that shows up as a zero instead of as the previous event's
                // photometry silently entering this event's Fisher matrix.
                //
                // NOTE: the untouched tail of tele[] is 0 (its constructed value), not the
                // -1 the old full clear wrote. Unobservable today -- tele is only read at
                // [0, ndw) -- but relevant if anything ever scans the whole array.
                for (int i = 0; i < ndw; ++i) {
                    l->timn[i] = 0.0;  l->magn[i] = 0.0; l->soux[i] = 0.0;  l->souy[i] = 0.0;
                    l->errm[i] = 0.0;  l->erra[i] = 0.0; l->tele[i] = -1;
                }

                ndw     = 0;   flag_det = 0;
                ndw_L   = 0;   ndw_R    = 0;
                nres5_L = 0; nres20_L = 0; nresPSF_L = 0; dsepMax_L = -1.0;
                nres5_R = 0; nres20_R = 0; nresPSF_R = 0; dsepMax_R = -1.0;
                flag_det_L = 0; flag_det_R = 0;
                flag0   = 0.0; flag1    = 0.0; flag2 = 0.0;
                flag0_L = 0.0; flag1_L  = 0.0; flag2_L = 0.0;
                flag0_R = 0.0; flag1_R  = 0.0; flag2_R = 0.0;
                chi1    = 0.0; chi2     = 0.0; chi3  = 0.0;
                chi1_L  = 0.0; chi2_L   = 0.0; chi3_L = 0.0;
                chi1_R  = 0.0; chi2_R   = 0.0; chi3_R = 0.0;
                chi1a   = 0.0; chi2a    = 0.0; chi3a = 0.0;
                chi1a_L = 0.0; chi2a_L  = 0.0; chi3a_L = 0.0;
                chi1a_R = 0.0; chi2a_R  = 0.0; chi3a_R = 0.0;
                def1p   = 0.0; s->def1c = 0.0; vsave = 0.0;
                def2p   = 0.0; s->def2c = 0.0;
                s->errM = 0.0; s->errA  = 0.0;
                dchiL   = 0.0; dchiP    = 0.0; dchiA = 0.0;
                dchiL_L = 0.0; dchiP_L  = 0.0; dchiA_L = 0.0;
                dchiL_R = 0.0; dchiP_R  = 0.0; dchiA_R = 0.0;

                dt   = 60.0;///days
                fdetRubin = 0.0;
                romanDetectable = false;

                for (int i = 0; i < M; ++i) {
                    Mpeak = s->magb[i] - 2.5 * std::log10(l->A0 * s->blend[i] + 1.0 - s->blend[i]);
//                        cout << "i=" << i << "  Mab=" << s->Mab[i] << "  Map=" << s->Map[i]
//                             << "  blend=" << s->blend[i] << "  Mpeak=" << Mpeak << endl;
                    if (i < 6) { // LSST ugrizy
                        if (Mpeak <= rubinDepthMed[i] and s->magb[i] > rubinDepthMed[i] - RUBIN_SATU_BELOW_M5)
                            fdetRubin += 1.0;
                    } else {     // i == 6, Roman F146 — single band, no ">=2 filters" bar applies
                        if (Mpeak <= thre[i] and s->magb[i] > satu[i])    romanDetectable = true;
                    }
                }
                rubinDetectable = (fdetRubin > 1.0); // at least 2 of the 6 LSST bands

                testL = RandR(0.0, 1.0);
                testR = RandR(0.0, 1.0);
                // Independent draws per survey — reusing one draw for both would correlate
                // the Rubin-accept and Roman-accept decisions for no physical reason. Each
                // draw is weighted by that survey's OWN blend fraction (Step B2), replacing
                // the old single test <= s->blend[2] (LSST r-band only, for both surveys).
                acceptRubin = rubinDetectable and (testL <= s->blend[2]);
                acceptRoman = romanDetectable and (testR <= s->blend[6]);

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                if (acceptRubin or acceptRoman) { //detectable by Rubin (>=2 LSST bands) or by Roman (F146)
                    cout << "************** DETECTABLE!!!!!! ********" << endl;
                    s->nsdet[s->nums] += 1.0;
                    flagf = 1;
                    // (The legacy magC/datC/BHLSSTMONTS demo dump that was gated here could never
                    // fire -- save < 0 with save = 0 -- and was removed in Deviation 78, with the
                    // random draw that fed it.)
    
                    gi = 0;
                    gi = 0; giR = 0;
                    for (double tim = float(0.0 * year - 100.0 - initial);  tim < float(10.0 * year + 100.0 + initial); tim = tim + dt) {
                        // Rubin's geometry (observer on Earth). The quantities derived below and
                        // shared across both branches -- Astar0/As0/As1, vs1/vs2, def1p/def2p,
                        // trajm/trajp, magni[] -- are all in this frame; the Roman branch
                        // recomputes the ones it needs in its own frame (Step H1).
                        lightcurve(*s, *l, *as, tim, 0);
                        Astar0   = double(s->ut0 * s->ut0 + 2.0) / std::sqrt(s->ut0 * s->ut0 * (s->ut0 * s->ut0 + 4.0)); //MAgnification equation
                        s->Astar = double(s->ut  * s->ut  + 2.0) / std::sqrt(s->ut  * s->ut  * (s->ut  * s->ut  + 4.0)); //MAgnification equation
                        As0      = double(Astar0   * s->blend[2] + 1.0 - s->blend[2]);
                        As1      = double(s->Astar * s->blend[2] + 1.0 - s->blend[2]); //LSST r-band
                        vs1      = double(s->mus1 + (s->def1c - def1p) / dt); //[mas/days]
                        vs2      = double(s->mus2 + (s->def2c - def2p) / dt); //[mas/days]
                        def1p    = s->def1c; //pervious
                        def2p    = s->def2c;

                        // Computed once per timestep (not just inside the LSST branch) since both
                        // instruments' astrometric chi-square terms need it, and it only depends on
                        // s->pos1b/pos2b/pos1c/pos2c, which lightcurve() already refreshed above.
                        trajm = std::sqrt(s->pos1b * s->pos1b + s->pos2b * s->pos2b); //stright + parallax
                        trajp = std::sqrt(s->pos1c * s->pos1c + s->pos2c * s->pos2c); //stright + parallax+lensing

    
                        for (int i = 0; i < M; ++i) {
                            magni0[i] = s->magb[i] - 2.5 * std::log10(Astar0   * s->blend[i] + 1.0 - s->blend[i]);
                            magni[i]  = s->magb[i] - 2.5 * std::log10(s->Astar * s->blend[i] + 1.0 - s->blend[i]);
                        }
                        // Rubin's representative-band model magnitude (RUBIN_REF_BANDS, config/parameters.h) --
                        // replaces the old hardcoded magni[2] (r-band) at the two use sites below.
                        // Reduces to exactly magni[2] when RUBIN_REF_BANDS = {2} (the default), since
                        // s->mbs[0]/s->fb[0] were built from the same combination in func_source.
                        double magniRubinRef = s->mbs[0] - 2.5 * std::log10(s->Astar * s->fb[0] + 1.0 - s->fb[0]);
                        sq = int(ls->ct[gi]);
                        sq  = int(ls->ct[gi]);
                        sqR = (nddR > 0 and giR < nddR) ? int(ro->ct[giR]) : -1;

                        // ---------------- LSST (ugrizy) ----------------
    
                        if (tim >= 0.0 and tim <= Tobs and tim >= ls->tim[int(ls->ct[0])] and tim <= ls->tim[int(ls->ct[ndd - 1])] and
                            gi < ndd and sq >= 0 and sq <= static_cast<int>(Nl) and tim >= ls->tim[sq]) {
    
                            fi = int(ls->filter[sq]);
    
                            // Deviation 73: this visit's own depth and saturation, not the SRD
                            // minimum -- the same depth that sets errg below.
                            const double m5v   = double(ls->sig5[sq]);
                            const double satuv = m5v - RUBIN_SATU_BELOW_M5;
                            if (magni[fi] >= satuv and magni[fi] <= m5v) {
                                errg = errlsstM(magni[fi], int(fi), m5v); //[mag]
                                errs = errlsstA(*ls, magniRubinRef); ///[mas]

                                // Step R1. Could Rubin have told the two images apart at THIS
                                // epoch? Inside the magnitude gate on purpose: the paper's
                                // criterion counts recorded data points, and an epoch the
                                // survey never recorded is not one.
                                {
                                    const ImagePair ip = imagePair(s->ut, l->tetE, s->magb[fi],
                                                                   s->blend[fi], m5v, satuv);
                                    if (ip.bothDetectable) {
                                        if (ip.sep >= RESOLVE_D_FAINT  * errs)          nres5_L   += 1;
                                        if (ip.sep >= RESOLVE_D_BRIGHT * errs)          nres20_L  += 1;
                                        if (ip.sep >= FWHM[fi] * ARCSEC_TO_MAS)         nresPSF_L += 1;
                                        if (ip.sep >  dsepMax_L)                        dsepMax_L  = ip.sep;
                                    }
                                }

                                deltaA = std::fabs(std::pow(10.0, -0.4 * errg) - 1.0) * (s->blend[fi] * s->Astar + 1.0 - s->blend[fi]);
                                magnio = magni[fi] + RandN(errg, NOISE_TRUNC_NSIGMA);
    
                                chi1 += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg)); //real
                                chi2 += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg)); //real no parallax
                                chi3 += std::fabs((magnio - s->magb[fi]) * (magnio - s->magb[fi]) / (errg * errg)); //baseline
                                chi1_L += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg));
                                chi2_L += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg));
                                chi3_L += std::fabs((magnio - s->magb[fi]) * (magnio - s->magb[fi]) / (errg * errg));

                                // erra is the per-coordinate sigma (Deviation 71; was drawn and
                                // divided as if each coordinate had variance 2 errs^2).
                                sil  = RandN(errs, NOISE_TRUNC_NSIGMA);
                                sil2 = RandN(errs, NOISE_TRUNC_NSIGMA);

                                chi1a += std::fabs((trajp + sil - trajp) * (trajp + sil - trajp) / (errs * errs)); //real trajectory
                                chi2a += std::fabs((trajp + sil - trajm) * (trajp + sil - trajm) / (errs * errs)); //real without lensing
                                chi3a += std::fabs( sil2  * sil2 / (errs * errs));
                                chi1a_L += std::fabs((trajp + sil - trajp) * (trajp + sil - trajp) / (errs * errs));
                                chi2a_L += std::fabs((trajp + sil - trajm) * (trajp + sil - trajm) / (errs * errs));
                                chi3a_L += std::fabs( sil2  * sil2 / (errs * errs));

                                CHECK(ndw < coun); // array-bounds guard — see note on `coun` sizing
                                l->timn[ndw] = tim;
                                l->magn[ndw] = magniRubinRef; // Rubin representative-band model value (RUBIN_REF_BANDS)
                                l->errm[ndw] = errg;
                                l->soux[ndw] = s->pos1c;
                                l->souy[ndw] = s->pos2c;
                                l->erra[ndw] = errs;
                                l->tele[ndw] = 0; // 0 = LSST
                                l->rseas[ndw] = -1; l->rroll[ndw] = -1;

                                // Step S1. Everything here was computed above for the
                                // detection test; nothing new is drawn. magnio in
                                // particular is the noisy datum chi1/chi2/chi3 just used.
                                if (dumpSpec.on)
                                    dumpBuf.push_back(DumpEpoch{
                                        tim, 0, int(fi),
                                        magnio, magni[fi], magni0[fi], errg,
                                        s->ut, s->ut0, s->Astar, Astar0,
                                        s->def1c, s->def2c, s->def1a, s->def2a,
                                        s->pos1b, s->pos2b, s->pos1c, s->pos2c,
                                        l->pos1,  l->pos2,  errs});
    
                                flag2 = 0.0;
                                if (std::fabs(magnio - s->magb[fi]) > std::fabs(OUTLIER_FLAG_NSIGMA * errg))    flag2 = 1.0;
                                if (ndw_L > 2 and float(flag0 + flag1 + flag2) > OUTLIER_RUN_THRESHOLD)   flag_det = 1;
                                flag0 = flag1;
                                flag1 = flag2;
                                flag2_L = flag2; // identical test; kept as an explicit Rubin-labeled copy
                                if (ndw_L > 2 and float(flag0_L + flag1_L + flag2_L) > OUTLIER_RUN_THRESHOLD)   flag_det_L = 1;
                                flag0_L = flag1_L;
                                flag1_L = flag2_L;
    
    
                                CHECK(sq >= 0);
                                CHECK(sq <= int(Nl - 1));
                                CHECK(fi >= 0);
                                CHECK(fi <= 5);
                                CHECK(errs >= 0.0);
                                CHECK(errg >= 0.0);
                                CHECK(deltaA >= 0.0);
                                CHECK(gi <= ndd);
    
                                s->errM += deltaA;
                                s->errA += errs;
                                vsave += std::sqrt(vs1 * vs1 + vs2 * vs2);

                                ndw += 1;
                                ndw_L += 1;
                            }//magnitude limit
                            gi += 1;
                        }
                        // ---------------- Roman (F146) — new ----------------
                        // Roman now feeds its own chi1_R/chi2_R/chi3_R (+ astrometric _R,
                        // placeholder error — see errsR below) and the joint chi1/chi2/chi3,
                        // alongside the Rubin-only _L versions built in the LSST branch above.
                        // See Step B1 of JOINT_FIT_REFACTOR_PLAN.md. s->errM/s->errA remain
                        // Rubin-only by design (ORIENTATION.md: they're specifically the
                        // LSST-only running-error accumulators, unrelated to detection).
                        if (nddR > 0 and tim >= 0.0 and tim <= Tobs and
                            tim >= ro->tim[int(ro->ct[0])] and tim <= ro->tim[int(ro->ct[nddR - 1])] and
                            giR < nddR and sqR >= 0 and sqR <= static_cast<int>(NlRoman) and tim >= ro->tim[sqR]) {

                            constexpr int fiR = 6; // F146 — the only Roman band modeled so far

                            // ---- Step H1: switch to Roman's observer position ----
                            // Everything above was computed with the observer on Earth, which is
                            // right for Rubin and wrong for Roman by the L2 offset. Roman sees a
                            // slightly different trajectory, so the magnification, the F146
                            // magnitude and the astrometric positions all have to be rebuilt here
                            // before any of them is recorded. `magni`/`magni0`/`trajm`/`trajp` are
                            // per-timestep scratch that the Rubin branch above has already
                            // consumed and that the next timestep overwrites, so they are
                            // rewritten in place rather than shadowed.
                            //
                            // No RNG is consumed by this call, so the random stream is untouched
                            // and L2_OFFSET_AU = 0 reproduces the pre-H1 run exactly.
                            lightcurve(*s, *l, *as, tim, 1);
                            {
                                const double Astar0R = double(s->ut0 * s->ut0 + 2.0)
                                                     / std::sqrt(s->ut0 * s->ut0 * (s->ut0 * s->ut0 + 4.0));
                                s->Astar = double(s->ut * s->ut + 2.0)
                                         / std::sqrt(s->ut * s->ut * (s->ut * s->ut + 4.0));
                                magni0[fiR] = s->magb[fiR] - 2.5 * std::log10(Astar0R  * s->blend[fiR] + 1.0 - s->blend[fiR]);
                                magni[fiR]  = s->magb[fiR] - 2.5 * std::log10(s->Astar * s->blend[fiR] + 1.0 - s->blend[fiR]);
                                trajm = std::sqrt(s->pos1b * s->pos1b + s->pos2b * s->pos2b);
                                trajp = std::sqrt(s->pos1c * s->pos1c + s->pos2c * s->pos2c);
                            }

                            if (magni[fiR] >= satu[fiR] and magni[fiR] <= thre[fiR]) {
                                errgR = errRomanM(*ro, magni[fiR]); //[mag] (Deviation 72)

                                magnioR = magni[fiR] + RandN(errgR, NOISE_TRUNC_NSIGMA);
                                chi1 += std::fabs((magnioR -   magni[fiR]) * (magnioR -   magni[fiR]) / (errgR * errgR));
                                chi2 += std::fabs((magnioR -  magni0[fiR]) * (magnioR -  magni0[fiR]) / (errgR * errgR));
                                chi3 += std::fabs((magnioR - s->magb[fiR]) * (magnioR - s->magb[fiR]) / (errgR * errgR));
                                chi1_R += std::fabs((magnioR -   magni[fiR]) * (magnioR -   magni[fiR]) / (errgR * errgR));
                                chi2_R += std::fabs((magnioR -  magni0[fiR]) * (magnioR -  magni0[fiR]) / (errgR * errgR));
                                chi3_R += std::fabs((magnioR - s->magb[fiR]) * (magnioR - s->magb[fiR]) / (errgR * errgR));

                                // Step H4: Roman's own per-exposure astrometric error, replacing the
                                // errlsstA() placeholder (Rubin's curve at Roman's magnitude, which had
                                // no reason to be right). Constants and sources in config/parameters.h; the model
                                // is per EXPOSURE, which is what one row of RomanBaseline.dat is.
                                errsR = errRomanA(magni[fiR]); //[mas]

                                // Step R1, Roman side. s->ut is Roman's OWN impact parameter
                                // here: lightcurve(..., 1) rebuilt the trajectory in the L2
                                // frame above, so this is not the Rubin value reused.
                                {
                                    const ImagePair ip = imagePair(s->ut, l->tetE, s->magb[fiR],
                                                                   s->blend[fiR], thre[fiR], satu[fiR]);
                                    if (ip.bothDetectable) {
                                        if (ip.sep >= RESOLVE_D_FAINT  * errsR)         nres5_R   += 1;
                                        if (ip.sep >= RESOLVE_D_BRIGHT * errsR)         nres20_R  += 1;
                                        if (ip.sep >= FWHM[fiR] * ARCSEC_TO_MAS)        nresPSF_R += 1;
                                        if (ip.sep >  dsepMax_R)                        dsepMax_R  = ip.sep;
                                    }
                                }

                                silR  = RandN(errsR, NOISE_TRUNC_NSIGMA);
                                sil2R = RandN(errsR, NOISE_TRUNC_NSIGMA);
                                chi1a += std::fabs((trajp + silR - trajp) * (trajp + silR - trajp) / (errsR * errsR));
                                chi2a += std::fabs((trajp + silR - trajm) * (trajp + silR - trajm) / (errsR * errsR));
                                chi3a += std::fabs( sil2R * sil2R / (errsR * errsR));
                                chi1a_R += std::fabs((trajp + silR - trajp) * (trajp + silR - trajp) / (errsR * errsR));
                                chi2a_R += std::fabs((trajp + silR - trajm) * (trajp + silR - trajm) / (errsR * errsR));
                                chi3a_R += std::fabs( sil2R * sil2R / (errsR * errsR));

                                flag2_R = 0.0;
                                if (std::fabs(magnioR - s->magb[fiR]) > std::fabs(OUTLIER_FLAG_NSIGMA * errgR))    flag2_R = 1.0;
                                if (ndw_R > 2 and float(flag0_R + flag1_R + flag2_R) > OUTLIER_RUN_THRESHOLD)   flag_det_R = 1;
                                flag0_R = flag1_R;
                                flag1_R = flag2_R;

                                CHECK(ndw < coun); // array-bounds guard — see note on `coun` sizing
                                l->timn[ndw] = tim;
                                l->magn[ndw] = magni[fiR]; // F146 magnitude — NOT magni[2]; FisherM's
                                                            // tt==1 branch compares against s.mbs[1]/s.fb[1],
                                                            // which are F146-based, so the reference point
                                                            // recorded here must be F146 too.
                                l->errm[ndw] = errgR;
                                l->soux[ndw] = s->pos1c;
                                l->souy[ndw] = s->pos2c;
                                // Step H4. This used to store `errs` -- the RUBIN astrometric error,
                                // left over from whichever Rubin epoch last set it, and in general from
                                // a different timestep entirely. So the astrometric Fisher matrix was
                                // being weighted by a stale value from the other telescope, not even by
                                // the errlsstA(magni[fiR]) the comment above it described: errsR was
                                // computed for the chi-squared terms and then thrown away. Now Roman's
                                // own per-exposure error is both used and stored.
                                l->erra[ndw] = errsR;
                                l->tele[ndw] = 1; // 1 = Roman/F146
                                // Season and roll of this exposure, for the astrometric noise
                                // variants' day blocks and frame groups (Deviation 71).
                                l->rseas[ndw] = sched.seasonOf(ro->tim[sqR]);
                                l->rroll[ndw] = ro->layout[sqR];
                                CHECK(l->rseas[ndw] >= 0);

                                // Step S1, Roman side. s->ut, s->def* and s->pos* are the
                                // L2-frame values lightcurve(..., 1) rebuilt above, not the
                                // Rubin ones from earlier in this timestep. Astar0 is
                                // recomputed from s->ut0 because the Roman branch's own
                                // Astar0R went out of scope before the magnitude gate.
                                if (dumpSpec.on)
                                    dumpBuf.push_back(DumpEpoch{
                                        tim, 1, fiR,
                                        magnioR, magni[fiR], magni0[fiR], errgR,
                                        s->ut, s->ut0, s->Astar, magnifOf(s->ut0),
                                        s->def1c, s->def2c, s->def1a, s->def2a,
                                        s->pos1b, s->pos2b, s->pos1c, s->pos2c,
                                        l->pos1,  l->pos2,  errsR});

                                CHECK(sqR >= 0);
                                CHECK(sqR <= int(NlRoman - 1));
                                CHECK(errgR >= 0.0);
                                CHECK(giR <= nddR);

                                ndw += 1;
                                ndw_R += 1;
                            }//magnitude limit
                            giR += 1;
                        }

                        // ---------------- Adaptive step size ----------------
                        // dt must respect whichever instrument's *next* epoch is sooner.
                        // Sizing dt to LSST's cadence alone (the original behavior) would
                        // silently step right over Roman's much denser epochs.
                        if ( tim < -50.0 or tim > (10.0 * year + 50.0)) {
                            dt = 60.0; //days
                        }
    
                        else {
                            bool lsstInWindow  = (tim >= ls->tim[int(ls->ct[0])] and tim <= ls->tim[int(ls->ct[ndd - 1])]);
                            bool romanInWindow = (nddR > 0 and tim >= ro->tim[int(ro->ct[0])] and tim <= ro->tim[int(ro->ct[nddR - 1])]);

                            if (lsstInWindow) {
                                if (gi > 0 and gi < ndd) cade = float(ls->tim[int(ls->ct[gi])] - ls->tim[int(ls->ct[gi - 1])]); //days
                                else cade = minc;
                            } else {
                                cade = TIME_STEP_OUTSIDE_LSST_DAYS; //days — matches the original "outside LSST window" fallback
                            }

                            if (romanInWindow) {
                                if (giR > 0 and giR < nddR) cadeR = float(ro->tim[int(ro->ct[giR])] - ro->tim[int(ro->ct[giR - 1])]); //days
                                else cadeR = mincR;
                            } else {
                                cadeR = cade; // Roman not active right now — don't let it constrain dt
                            }

                            dt = double(std::min(cade, cadeR));
                        }
    
                    }//end of loop time
                }// end of visible star

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH 
                for (int i = 0; i < nq; ++i) co->resu[i] = -1.0;
                // Same sentinel discipline for the per-survey results: an event that never
                // reaches FisherM must not inherit the previous event's sigmas.
                for (int q = 0; q < NSURV; ++q) {
                    co->okA[q] = 0; co->okB[q] = 0; co->nepochA[q] = 0;
                    co->flagi = 0;   // Deviation 78: was stale on uncharacterised events
                    co->condA[q] = -1.0; co->condB[q] = -1.0;
                    co->relMl[q] = -1.0;
                    for (int k = 0; k < Nx; ++k) co->Era[q][k] = -1.0;
                    for (int k = 0; k < Ny; ++k) co->Erb[q][k] = -1.0;
                    for (int v = 0; v < NAVAR; ++v) {
                        co->okBV[v][q] = 0; co->condBV[v][q] = -1.0; co->relMlV[v][q] = -1.0;
                        for (int k = 0; k < Ny; ++k) co->ErbV[v][q][k] = -1.0;
                    }
                }
    
                FFG[0] = 0;
                FFG[1] = 0;
                FFG[2] = 0;
                detL = 0; detR = 0; detJ = 0;

                //cout << "flagf: " << flagf << "ndw: " << ndw << endl;

                if (flagf == 0 or ndw <= 2) {
                    errg    = errlsstM(s->magb[2], 2, double(RUBIN_R_DEPTH5_FALLBACK)); //r-band
                    s->errA = errlsstA(*ls, s->magb[2]); //r-band
                    s->errM = std::fabs(std::pow(10.0, - 0.4 * errg) - 1.0); //r-band
                    dchiL = 0.0;
                    dchiP = 0.0;
                    dchiA = 0.0;
                    dchiL_L = 0.0; dchiP_L = 0.0; dchiA_L = 0.0;
                    dchiL_R = 0.0; dchiP_R = 0.0; dchiA_R = 0.0;
                    vsave = s->mus;
                }

                // Step H3's second forecast for this event. Declared one scope out from the
                // detection block below, because the row is written further down, where the
                // satellite observable and the coverage counts have been computed.
                // Step H3's no-satellite forecast for this same event. Declared here, before
                // the detection block, so the row write further down can see them whether or
                // not the Fisher step ran. -1.0 is the not-measured sentinel used everywhere
                // else in this code: a sigma that never inverted is not a large sigma.
                double sigtE_ns   = -1.0, sigpiE_ns  = -1.0, sigpiER_ns = -1.0;
                double sigtetE_ns = -1.0, sigpiEb_ns = -1.0, relMl_ns   = -1.0;
                double condA_ns   = -1.0, condB_ns   = -1.0;
                int    okNS       = 0,    okNSb      = 0;

                if (flagf > 0 and ndw > 2) { //if star is visible
                    cout << "************** DETECTABLE!!!!!! ********" << endl;
                    icon +=1;
                    vsave   = double(vsave   / (ndw + 0.000065645));
                    s->errM = double(s->errM / (ndw + 0.000065645));
                    s->errA = double(s->errA / (ndw + 0.000065645));
                    // Step H7. The lensing statistic is SIGNED: chi3 is the flat-baseline
                    // residual and chi1 the lensing-model residual, so chi3 - chi1 > 0 means
                    // the lensing model fits better and is the only direction that can count
                    // as evidence of lensing. The previous fabs() meant a lensing model that
                    // fitted WORSE than the baseline by more than the bar would have been
                    // reported as a detection. That never bit in practice, because the data
                    // are lensed and chi3 > chi1 for any event with real signal, but it also
                    // breaks the monotonicity H7 exists to establish: with fabs, an instrument
                    // whose epochs all fall outside the event contributes a small NEGATIVE
                    // difference from noise, and |dchi_L + dchi_R| can then fall below
                    // max(|dchi_L|, |dchi_R|). Signed, the sum is exact and monotone.
                    //
                    // dchiP and dchiA are not detection tests (nothing thresholds them); since
                    // Deviation 78 they are signed like dchiL, so a negative value (the simpler
                    // model fitting better) is visible rather than folded into a "size".
                    dchiL   = chi3  - chi1;             //lensing_effect (signed)
                    // Deviation 78: signed like dchiL (positive = the full model fits better).
                    dchiP   = chi2  - chi1;             //parallax_effect (signed)
                    dchiA   = chi2a - chi1a;            //deflection_effect (signed)
                    dchiL_L = chi3_L - chi1_L;
                    dchiP_L = (chi2_L  - chi1_L);   // signed (Deviation 78)
                    dchiA_L = (chi2a_L - chi1a_L);   // signed (Deviation 78)
                    dchiL_R = chi3_R - chi1_R;
                    dchiP_R = (chi2_R  - chi1_R);   // signed (Deviation 78)
                    dchiA_R = (chi2a_R - chi1a_R);   // signed (Deviation 78)
                    // flag_det is the JOINT run-test flag the table reports (Deviation 78): it was set in
                    // the Rubin branch only, so a Roman-only persistent signal left it at 0.
                    flag_det = (flag_det_L > 0 or flag_det_R > 0) ? 1 : 0;

                    // Three independent detection tests. detL/detR use each instrument's own
                    // epoch count (ndw_L/ndw_R) and run-test result — mixing in the joint ndw
                    // here would let Roman's dense epochs silently raise the bar for a purely
                    // Rubin-driven signal (and vice versa). detJ requires a persistent run in
                    // EITHER instrument's own cadence (flag_det_L or flag_det_R) rather than an
                    // interleaved run across two different cadences, which wouldn't mean anything
                    // (see Step B1 teaching brief). FFG[0] — the boolean that gates FisherM below,
                    // preserving the existing detection-then-Fisher ordering (plan §0.4) — is the
                    // union of all three, so an event Rubin alone clearly detects is never dropped
                    // from characterization just because the joint-scaled threshold happens to miss.
                    if (s->FWHM < Tobs and dchiL_L > cfg.dchiDet and flag_det_L > 0 and ndw_L > 10) detL = 1;
                    if (s->FWHM < Tobs and dchiL_R > cfg.dchiDet and flag_det_R > 0 and ndw_R > 10) detR = 1;
                    if (s->FWHM < Tobs and dchiL   > cfg.dchiDet and (flag_det_L > 0 or flag_det_R > 0) and ndw > 10) detJ = 1;

                    // Monotonicity. Adding data to an analysis cannot destroy information, so
                    // if either survey alone clears its detection bar, the combined stream --
                    // which contains that survey's data in full, plus more -- must clear its
                    // own.
                    //
                    // Since Step H7 this holds BY CONSTRUCTION rather than by patch. chi1 and
                    // chi3 are accumulated over both instruments' epochs in the same loops that
                    // fill chi1_L/chi1_R and chi3_L/chi3_R, so
                    //
                    //     chi1 = chi1_L + chi1_R   and   chi3 = chi3_L + chi3_R
                    //
                    // exactly, hence dchiL = dchiL_L + dchiL_R with the signed statistic above.
                    // All three tests now use the SAME fixed bar, so if either survey's own
                    // dchi clears it, the sum -- which is that dchi plus a quantity that is
                    // positive whenever the other survey sees any signal at all -- clears it
                    // too. The auxiliary gates agree: flag_det_L > 0 implies (flag_det_L or
                    // flag_det_R), and ndw_L > 10 implies ndw = ndw_L + ndw_R > 10.
                    //
                    // The old bar was 2*ndw, a threshold on the MEAN per-epoch chi-squared
                    // improvement rather than on total significance, which made pooling a
                    // survey with many low-signal epochs RAISE the joint bar without adding
                    // signal. In Roman's footprint that lifted the joint bar ~4,728 above
                    // Roman's own and vetoed 21.3% of all detections there.
                    //
                    // detJ_raw and the monotone patch are retained deliberately (H7 acceptance
                    // criterion 4) so the anomaly rate stays MEASURABLE and can be shown to
                    // have gone to zero rather than been hidden. If DET_ANOMALY is ever
                    // non-zero again, the construction above has been broken and the counter
                    // is how that gets noticed.
                    const int detJ_raw = detJ;
                    if (detL or detR) detJ = 1;

                    if (detL or detR or detJ) { //lensing — detected by Rubin, Roman, or the joint test
                        FFG[0] = 1; //Lensing
                        nlens += FFG[0];
                        FisherM(*s, *l, *as, *co, ndw);

                        if (co->flagi > 0) {
                            // flagi is always +1: FisherM's F*F^-1 checks are commented out
                            // (Deviation 40). The conditioning test that works is okA.
                            nerr += co->okA[SJOINT] ? 1.0 : 0.0;
                            ErrorCal(*co, *l, *s);

                            std::ofstream fil0_append(fnLDt, std::ios::app);
                            fil0_append << std::fixed << std::setprecision(5)
                                        << l->Ml       << " " << l->Dl       << " " << l->mul       << " " << s->fb[0]     << " " << s->mbs[0]   << " "
          <<std::setprecision(7)           << l->tE       << " " << l->murel    << " " << l->u0        << " " << s->lon       << " " << s->lat      << " "
                                        << l->piE      << " " << l->tetE     << " " << co->resu[1]  << " " << co->resu[2]  << " " << co->resu[3] << " "
                                        << co->resu[5] << " " << co->resu[9] << " " << co->resu[10] << " " << co->resu[13] << "\n";
                            fil0_append.close();
                        }

                        // -----------------------------------------------------------------
                        // Step H3. The satellite-parallax experiment as a CONTROLLED
                        // comparison on one event, instead of a difference between two runs.
                        //
                        // Why it must be done here. The obvious design -- two runs differing
                        // only in L2_OFFSET_AU, matched by row -- cannot work. The RNG is one
                        // stream and the per-event path draws conditionally on `acceptRubin
                        // or acceptRoman`, so moving the observer changes what is detectable,
                        // the streams fork at the first footprint sightline, and the two
                        // runs' DETECTED populations then differ. Measured on the v3 pair,
                        // that selection difference shifts the median sigma_tE -- which the
                        // satellite cannot physically touch -- by 15%, several times the
                        // effect being looked for. DEVIATIONS.md 35.
                        //
                        // Here the event is fixed: same draw, same epochs, same photometry,
                        // only the observer moves. satScale is the single knob carrying
                        // L2_OFFSET_AU into lightcurve(), and FisherM rebuilds the light
                        // curve from it, so flipping it and recomputing gives the forecast
                        // THIS event would have had with Roman at Earth.
                        //
                        // Re-evaluating is legitimate for an event whose data were generated
                        // with the offset on: the Fisher matrix is built from model
                        // derivatives at the true parameters, not from realised noise. The
                        // question is what each observing geometry can constrain.
                        if (cfg.pairSat) {
                            const double keepScale = as->satScale;
                            as->satScale = 0.0;
                            FisherM(*s, *l, *as, *coNS, ndw);
                            if (coNS->flagi > 0) ErrorCal(*coNS, *l, *s);
                            as->satScale = keepScale;

                            // Both halves of the forecast are recorded, because moving the
                            // observer acts on both and they carry different physics.
                            //
                            //   PHOTOMETRIC (okA, Era): the magnification depends on |u|, and
                            //   the two observers see different |u| at the same instant. That
                            //   difference is the simultaneous parallax baseline, and it lands
                            //   on piE.
                            //
                            //   ASTROMETRIC (okB, Erb): the centroid deflection is a vector,
                            //   theta_E * u/(u^2+2), so moving the observer changes its
                            //   direction as well as its size. Erb[0] is sigma(theta_E) and
                            //   Erb[3] the astrometric route to piE, which is independent of
                            //   the photometric one.
                            //
                            // relMl combines them: Ml = theta_E/(kappa piE) is the quantity a
                            // microlensing survey actually wants, and it needs one observable
                            // from each matrix, so it is the only place where a change in
                            // either shows up as a change in the science.
                            //
                            // SROMAN as well as SJOINT: the L2 offset is Roman's geometry, so
                            // Roman's own matrix is where any effect must appear first and
                            // undiluted by Rubin's epochs.
                            //
                            // Condition numbers for both matrices and both observers, so the
                            // conditioning question can be answered from the data instead of
                            // hypothesised. This is the diagnostic OPEN_ITEMS asked for.
                            okNS  = coNS->okA[SJOINT];
                            okNSb = coNS->okB[SJOINT];
                            sigtE_ns   = okNS  ? coNS->Era[SJOINT][1] : -1.0;
                            sigpiE_ns  = okNS  ? coNS->Era[SJOINT][3] : -1.0;
                            sigpiER_ns = coNS->okA[SROMAN] ? coNS->Era[SROMAN][3] : -1.0;
                            sigtetE_ns = okNSb ? coNS->Erb[SJOINT][0] : -1.0;
                            sigpiEb_ns = okNSb ? coNS->Erb[SJOINT][3] : -1.0;
                            relMl_ns   = coNS->relMl[SJOINT];
                            condA_ns   = coNS->condA[SJOINT];
                            condB_ns   = coNS->condB[SJOINT];
                        }
                    }

                    // Detection taxonomy (DetClass, Bulge.h). The joint fit is what makes a
                    // detection meaningful -- it sees strictly more data than either survey
                    // alone -- so the classes are distinguished by which telescopes ALSO
                    // detect the event unaided. DET_ANOMALY catches the case that should be
                    // impossible, a single-telescope detection the joint test misses.
                    dclsEvent = detClass(detL, detR, detJ);
                    const int dcls = dclsEvent;
                    nDetClass[dcls]    += 1;
                    NDetClassTot[dcls] += 1;
                    NDetClassTE[gg][dcls] += 1;
                    // Diagnostic only: how often the raw joint test contradicted a
                    // single-survey detection, before monotonicity was imposed.
                    // Both counters, and the run total is the one that was missing: it is
                    // read at the end of the run to decide whether to explain DET_ANOMALY,
                    // but was never incremented, so every run has reported exactly zero
                    // anomalies since the counter was introduced -- including the 2026-08-30
                    // production run, whose "0 ANOMALY" is this bug and not a measurement.
                    // The per-sightline counter was always right, which is why the anomaly
                    // was visible in the log all along and invisible in every summary.
                    // Deviation 79: the joint delta-chi2 is accumulated separately from the two
                    // per-survey ones; it must equal their sum on every event.
                    if (std::fabs((dchiL_L + dchiL_R) - dchiL) > 1e-6 * std::max(1.0, std::fabs(dchiL))) {
                        if (nDchiMismatch < 10)
                            std::cerr << "DCHI_MISMATCH sightline " << iScan << " joint " << dchiL
                                      << " L+R " << dchiL_L + dchiL_R << "\n";
                        nDchiMismatch += 1;
                    }
                    if (!detJ_raw and (detL or detR)) {
                        nDetClass[DET_ANOMALY]    += 1;
                        NDetClassTot[DET_ANOMALY] += 1;
                        // Deviation 79: log everything needed to explain it. The joint delta-chi2
                        // is the SUM of the two surveys' (each accumulated over its own epochs), and
                        // either term can be negative through noise -- a survey whose data happen to
                        // fit a flat baseline slightly better than the true model -- so a single-
                        // survey detection just over the bar can leave the sum just under it. That
                        // is physics, not a sign bug; "sum - joint" must be ~0 if the bookkeeping is
                        // right, and a non-zero value there WOULD be a bug.
                        std::cerr << std::setprecision(10)
                                  << "DET_ANOMALY_DETAIL sightline " << iScan << " lon " << s->lon
                                  << " lat " << s->lat << " | dchiL_L " << dchiL_L << " dchiL_R "
                                  << dchiL_R << " sum " << dchiL_L + dchiL_R << " joint " << dchiL
                                  << " (sum - joint " << (dchiL_L + dchiL_R) - dchiL << ") bar "
                                  << cfg.dchiDet << " | ndw_L " << ndw_L << " ndw_R " << ndw_R
                                  << " ndw " << ndw << " | flag_det_L " << flag_det_L
                                  << " flag_det_R " << flag_det_R << " | detL " << detL << " detR "
                                  << detR << " | tE " << l->tE << " u0 " << l->u0 << " t0 " << l->t0
                                  << "\n";
                    }

                    // Detected-event count for this tE bin. The joint detection is the one
                    // that defines "detected" here, per the taxonomy above; the per-class
                    // breakdown lives in nDetClass. gg was computed for every simulated event
                    // before the detection test, so nstE (the denominator) counts all draws
                    // and this counts the numerator -- which is what makes EFF an efficiency.
                    if (detJ) l->ndtE[gg] += 1.0;

                    // The numerators for the other six axes, using the bin indices computed
                    // for this draw before the detection test. Same definition of "detected"
                    // as the tE curve -- the joint test -- so all seven efficiencies describe
                    // the same thing and can be read side by side.
                    if (detJ) {
                        l->NdMl[ss] += 1.0;
                        l->Ndpi[qq] += 1.0;
                        l->Ndu0[ww] += 1.0;
                        l->Ndmu[vv] += 1.0;
                        l->Ndmb[zz] += 1.0;
                        l->Ndfb[pp] += 1.0;

                        // Nhalo and Nself were declared, written out, and never counted, so
                        // the EfLMC writer divided 0 by 0 and every row of every EfLMC file
                        // ever produced ends in two "-nan" columns. Counted here, with the
                        // definitions stated rather than guessed:
                        //   Nhalo[1]/Nhalo[0]  fraction of detections whose LENS is a halo
                        //                      star -- the population a MACHO search cares
                        //                      about, and negligible toward the bulge.
                        //   Nself[1]/Nself[0]  fraction that are bulge self-lensing: bulge
                        //                      lens AND bulge source, the dominant channel
                        //                      here and the one whose kinematics set tE.
                        l->Nhalo[0] += 1.0;
                        l->Nself[0] += 1.0;
                        if (l->struc == GalacticComponent::HALO) l->Nhalo[1] += 1.0;
                        if (l->struc == GalacticComponent::BULGE and
                            s->struc == GalacticComponent::BULGE) l->Nself[1] += 1.0;
                    }
                }
//            CHECK(gg >= 0);

            // Deviation 76: the observed peak, and the gap geometry measured from it.
            const auto [t0obs, uminObs] = observedPeak(*s, *l, *as);

            records.push_back(EventRecord{
                icon, static_cast<int>(FFG[0]),
                l->tE, l->RE/AU, l->piE, l->tetE, l->Vt, l->u0, l->Ml,
                s->opt*1.0e6, l->Dl, s->Ds, l->vl, s->vs,
                s->mbs[0], s->fb[0],
                gg, static_cast<int>(l->struc),
                s->FWHM/year, vsave/s->mus, l->DeltaT/s->errA, l->murel*year,
                co->resu[0], co->resu[1], co->resu[2], co->resu[3], co->resu[5],
                co->resu[9], co->resu[10], co->resu[13], co->resu[14],
                s->Map[2], s->nsbl[2],
                co->flagi, s->Ai[2],
                ndw_L, ndw_R,
                detL, detR, detJ,
                co->okA[SJOINT], co->okA[SRUBIN], co->okA[SROMAN],
                co->Era[SJOINT][1], co->Era[SRUBIN][1], co->Era[SROMAN][1],   //sigma(tE)
                co->Era[SJOINT][3], co->Era[SRUBIN][3], co->Era[SROMAN][3],   //sigma(piE)
                co->Erb[SJOINT][0], co->Erb[SRUBIN][0], co->Erb[SROMAN][0],   //sigma(tetE)
                dclsEvent,
                synergyClass(*co),
                co->condA[SJOINT], co->condA[SRUBIN], co->condA[SROMAN],
                // ---- the rest of the row (Step D1) ----
                l->t0, s->xi, s->lon, s->lat,
                s->mbs[1], s->fb[1],
                {s->magb[0], s->magb[1], s->magb[2], s->magb[3], s->magb[4], s->magb[5], s->magb[6]},
                {s->blend[0], s->blend[1], s->blend[2], s->blend[3], s->blend[4], s->blend[5], s->blend[6]},
                co->relMl[SJOINT], co->relMl[SRUBIN], co->relMl[SROMAN],
                co->okB[SJOINT], co->okB[SRUBIN], co->okB[SROMAN],
                co->condB[SJOINT], co->condB[SRUBIN], co->condB[SROMAN],
                sched.dtToSeasonEdge(t0obs), sched.zone(t0obs),
                nres5_L, nres20_L, nresPSF_L, nres5_R, nres20_R, nresPSF_R,
                dsepMax_L, dsepMax_R
            });
   
            // ------------------------------------------------------------------------------
            // Step H2. Two quantities that make the satellite-parallax effect visible in the
            // table instead of only implicit in the Fisher matrix.
            //
            // du_sat  the observer separation in Einstein radii at t0, piE * D_perp / AU.
            //         This is the amplitude of the satellite effect for this event and the
            //         natural x-axis of every Step H3 figure. It is NOT L2_OFFSET_AU * piE:
            //         D_perp is the separation projected perpendicular to the line of sight
            //         and runs ~0.87-0.99 of the full L2 offset around the year, so the
            //         simple product is a ceiling (DEVIATIONS.md 28.2). Computed by asking
            //         lightcurve() for both observers rather than re-deriving the projection
            //         here, so the two can never drift apart.
            //
            // nepL_pk, nepR_pk  epochs from each survey within +-2 tE of t0, i.e. while the
            //         event is actually magnified. Satellite parallax needs CONTEMPORANEOUS
            //         coverage: an event Roman saw in season 3 and Rubin saw in season 7 has
            //         none of it, however large ndw_L and ndw_R are. The plan asked for a
            //         single `nep_both` flag; two counts are the same cost and strictly more
            //         informative, and the flag is just (nepL_pk > 0 and nepR_pk > 0).
            //
            // Safe to call lightcurve() here: FisherM has already run, and the only state it
            // touches (s->ux/uy, s->def*, s->pos*, as->ue_n*) is not read by the row written
            // below and is recomputed from scratch by the next draw.
            // ------------------------------------------------------------------------------
            double duSat = 0.0;
            {
                lightcurve(*s, *l, *as, l->t0, 0);
                const double r1 = as->ue_n1, r2 = as->ue_n2;
                lightcurve(*s, *l, *as, l->t0, 1);
                const double dn1 = as->ue_n1 - r1, dn2 = as->ue_n2 - r2;
                duSat = l->piE * std::sqrt(dn1 * dn1 + dn2 * dn2);
            }
            int nepLpk = 0, nepRpk = 0;
            {
                const double win = 2.0 * l->tE;
                for (int i = 0; i < ndw; ++i) {
                    if (std::fabs(l->timn[i] - t0obs) > win) continue;   // Deviation 76
                    if (int(l->tele[i]) == 1) nepRpk += 1;
                    else                      nepLpk += 1;
                }
            }

            filg_in.open(testf, std::ios::app);
            filg_in << icon           << " " << FFG[0]         << " " << l->tE                      << " "
                    << l->RE / AU     << " " << l->piE         << " " << l->tetE                    << " "
                    << l->Vt          << " " << l->u0          << " " << l->Ml                      << " " 
                    << s->opt * 1.0e6 << " " << l->Dl          << " " << s->Ds                      << " "
                    << l->vl          << " " << s->vs          << " " << s->mbs[0]                  << " "
                    << s->fb[0]       << " " << gg             << " " << static_cast<int>(l->struc) << " "
                    << s->FWHM / year << " " << vsave / s->mus << " " << l->DeltaT / s->errA        << " " << l->murel * year << " "
                    << co->resu[0]    << " " << co->resu[1]    << " " << co->resu[2]                << " " 
                    << co->resu[3]    << " " << co->resu[5]    << " "
                    << co->resu[9]    << " " << co->resu[10]   << " " << co->resu[13]               << " " << co->resu[14]    << " "
                    << s->Map[2]      << " " << s->nsbl[2]     << " " << co->flagi                  << " " << s->Ai[2]        << " "
                    // per-survey bookkeeping (Step C5), appended so existing column indices hold
                    << ndw_L << " " << ndw_R << " "
                    << detL  << " " << detR  << " " << detJ << " "
                    << co->okA[SJOINT] << " " << co->okA[SRUBIN] << " " << co->okA[SROMAN] << " "
                    << co->Era[SJOINT][1] << " " << co->Era[SRUBIN][1] << " " << co->Era[SROMAN][1] << " "
                    << co->Era[SJOINT][3] << " " << co->Era[SRUBIN][3] << " " << co->Era[SROMAN][3] << " "
                    << co->Erb[SJOINT][0] << " " << co->Erb[SRUBIN][0] << " " << co->Erb[SROMAN][0] << " "
                    << dclsEvent << " " << synergyClass(*co) << " "
                    << co->condA[SJOINT] << " " << co->condA[SRUBIN] << " " << co->condA[SROMAN] << " "
                    // the rest of the row (Step D1) -- appended, so columns 1-57 keep their indices
                    << l->t0 << " " << s->xi << " " << s->lon << " " << s->lat << " "
                    << s->mbs[1] << " " << s->fb[1] << " ";
            for (int i = 0; i < M; ++i) filg_in << s->magb[i]  << " ";
            for (int i = 0; i < M; ++i) filg_in << s->blend[i] << " ";
            filg_in << co->relMl[SJOINT] << " " << co->relMl[SRUBIN] << " " << co->relMl[SROMAN] << " "
                    << co->okB[SJOINT]   << " " << co->okB[SRUBIN]   << " " << co->okB[SROMAN]   << " "
                    << co->condB[SJOINT] << " " << co->condB[SRUBIN] << " " << co->condB[SROMAN] << " "
                    // Gap geometry. dt_edge is NEGATIVE when t0 fell inside a Roman season;
                    // t0zone distinguishes a mid-mission gap (1) from before-launch/after-end
                    // (2), which must never be pooled -- only the former is gap-filling.
                    << sched.dtToSeasonEdge(t0obs) << " " << sched.zone(t0obs) << " "
                    // Sky area this event's sightline stands for, deg^2 (Step E1). Constant
                    // across an unstratified run; NOT constant once --stride-roman is used,
                    // and then any statistic pooled over sightlines must weight by it.
                    << wArea << " "
                    // Step H2: the satellite-parallax observable and contemporaneous coverage.
                    << duSat << " " << nepLpk << " " << nepRpk << " "
                    // Step R1: resolving the two images. Counts of qualifying epochs per
                    // survey, then the largest separation reached while both were detectable
                    // (-1 = never). The three bars differ only in what counts as "resolved".
                    << nres5_L << " " << nres20_L << " " << nresPSF_L << " " << dsepMax_L << " "
                    << nres5_R << " " << nres20_R << " " << nresPSF_R << " " << dsepMax_R << " "
                    // Deviation 71: astrometric noise variants N and P (joint, Roman).
                    << co->ErbV[AV_N][SJOINT][0] << " " << co->ErbV[AV_N][SROMAN][0] << " "
                    << co->ErbV[AV_P][SJOINT][0] << " " << co->ErbV[AV_P][SROMAN][0] << " "
                    << co->relMlV[AV_N][SJOINT] << " " << co->relMlV[AV_N][SROMAN] << " "
                    << co->relMlV[AV_P][SJOINT] << " " << co->relMlV[AV_P][SROMAN] << " "
                    << co->okBV[AV_N][SJOINT] << " " << co->okBV[AV_N][SROMAN] << " "
                    << co->okBV[AV_P][SJOINT] << " " << co->okBV[AV_P][SROMAN] << " "
                    << int(l->luminous) << " " << s->fLens[0] << " " << s->fLens[1] << " "
                    << t0obs << " " << uminObs << "\n";
            filg_in.close();

            // ------------------------------------------------------------------------------
            // Step S1. Commit this draw's buffered light curve if it fills a requested
            // sample class. Here, and not in the time loop, because detL/detR/detJ and the
            // Fisher sigmas -- which is what the classes are defined in terms of -- do not
            // exist until now.
            //
            // First match wins: an event is written once, under the first class in the spec
            // file whose quota is still open. Otherwise a `both` event would also land in
            // `any` and the same light curve would be drawn twice in one figure.
            // ------------------------------------------------------------------------------
            if (dumpSpec.on and not dumpBuf.empty()) {
                double maxShift = 0.0;
                for (const DumpEpoch& e : dumpBuf)
                    maxShift = std::max(maxShift,
                                        std::sqrt(e.def1c * e.def1c + e.def2c * e.def2c));

                const SampleFacts facts{detL, detR, detJ, l->tE, sched.zone(t0obs),
                                        nepLpk, nepRpk, ndw_L, ndw_R,
                                        co->okB[SROMAN], maxShift};

                for (SampleClass& c : dumpSpec.classes) {
                    if (c.kept >= c.quota)        continue;
                    if (!sampleMatch(c, facts))   continue;

                    std::ostringstream pr;
                    pr << std::setprecision(10)
                       << "# Sample event for class '" << c.name << "'.\n"
                       << "# One key per line. Angles in deg, times in days, masses in Msun,\n"
                       << "# distances in kpc, angular scales in mas. A sigma of -1 means the\n"
                       << "# parameter was NOT measured by that survey (inactive or no epochs),\n"
                       << "# never that it was measured to be -1.\n"
                       << "class "        << c.name        << "\n"
                       << "population "   << gPop->name    << "\n"
                       << "lon "          << s->lon        << "\n"
                       << "lat "          << s->lat        << "\n"
                       << "tE "           << l->tE         << "\n"
                       << "t0 "           << l->t0         << "\n"
                       << "u0 "           << l->u0         << "\n"
                       << "xi "           << s->xi         << "\n"
                       << "piE "          << l->piE        << "\n"
                       << "tetE "         << l->tetE       << "\n"
                       << "Ml "           << l->Ml         << "\n"
                       << "Dl "           << l->Dl         << "\n"
                       << "Ds "           << s->Ds         << "\n"
                       << "Vt "           << l->Vt         << "\n"
                       << "murel_yr "     << l->murel * year << "\n"
                       << "mus1 "         << s->mus1       << "\n"
                       << "mus2 "         << s->mus2       << "\n"
                       << "mul1 "         << l->mul1       << "\n"
                       << "mul2 "         << l->mul2       << "\n"
                       << "lens_struc "   << int(l->struc) << "\n"
                       << "mbs0 "         << s->mbs[0]     << "\n"
                       << "fb0 "          << s->fb[0]      << "\n"
                       << "mbs1 "         << s->mbs[1]     << "\n"
                       << "fb1 "          << s->fb[1]      << "\n";
                    pr << "magb";  for (int i = 0; i < M; ++i) pr << " " << s->magb[i];
                    pr << "\nblend"; for (int i = 0; i < M; ++i) pr << " " << s->blend[i];
                    pr << "\n"
                       << "ndw_L "        << ndw_L         << "\n"
                       << "ndw_R "        << ndw_R         << "\n"
                       << "nep_pk_L "     << nepLpk        << "\n"
                       << "nep_pk_R "     << nepRpk        << "\n"
                       << "detL "         << detL          << "\n"
                       << "detR "         << detR          << "\n"
                       << "detJ "         << detJ          << "\n"
                       << "dt_edge "      << sched.dtToSeasonEdge(t0obs) << "\n"
                       << "t0zone "       << sched.zone(t0obs)           << "\n"
                       << "du_sat "       << duSat         << "\n"
                       << "max_shift "    << maxShift      << "\n"
                       << "dsep_max_L "   << dsepMax_L     << "\n"
                       << "dsep_max_R "   << dsepMax_R     << "\n"
                       << "okA_J "  << co->okA[SJOINT] << " okA_L " << co->okA[SRUBIN]
                       << " okA_R " << co->okA[SROMAN] << "\n"
                       << "okB_J "  << co->okB[SJOINT] << " okB_L " << co->okB[SRUBIN]
                       << " okB_R " << co->okB[SROMAN] << "\n"
                       << "sigtE_J "   << co->Era[SJOINT][1] << " sigtE_L " << co->Era[SRUBIN][1]
                       << " sigtE_R "  << co->Era[SROMAN][1] << "\n"
                       << "sigpiE_J "  << co->Era[SJOINT][3] << " sigpiE_L " << co->Era[SRUBIN][3]
                       << " sigpiE_R " << co->Era[SROMAN][3] << "\n"
                       << "sigtetE_J "  << co->Erb[SJOINT][0] << " sigtetE_L " << co->Erb[SRUBIN][0]
                       << " sigtetE_R " << co->Erb[SROMAN][0] << "\n"
                       << "relMl_J "  << co->relMl[SJOINT] << " relMl_L " << co->relMl[SRUBIN]
                       << " relMl_R " << co->relMl[SROMAN] << "\n";

                    std::ostringstream idss;
                    idss << std::setw(3) << std::setfill('0') << (++dumpSeq);
                    writeSampleEvent(dumpSpec, c.name, idss.str(), dumpBuf, pr.str(),
                                     *s, *l, *as);
                    c.kept += 1;
                    std::cout << "  [sample] " << c.name << " " << idss.str()
                              << "  tE=" << l->tE << "d  Ml=" << l->Ml << "Msun"
                              << "  ndw_L=" << ndw_L << " ndw_R=" << ndw_R
                              << "  (" << c.kept << "/" << c.quota << ")" << std::endl;
                    break;
                }
            }

            // Step H3. One row per DETECTED event, carrying both forecasts for that same
            // event. Written here rather than beside the Fisher call because duSat and the
            // contemporaneous-coverage counts are computed above, and they are the axes
            // every H3 figure uses.
            if (cfg.pairSat and (detL or detR or detJ)) {
                std::ofstream fpair(fnPair, std::ios::app);
                const int okAs = co->okA[SJOINT], okBs = co->okB[SJOINT];
                fpair << std::setprecision(7)
                      << s->lon << " " << s->lat << " "
                      << l->tE  << " " << l->u0  << " " << l->piE << " " << l->tetE << " "
                      << duSat << " "
                      << okAs << " " << okNS << " " << okBs << " " << okNSb << " "
                      << (okAs ? co->Era[SJOINT][1] : -1.0) << " " << sigtE_ns   << " "
                      << (okAs ? co->Era[SJOINT][3] : -1.0) << " " << sigpiE_ns  << " "
                      << (co->okA[SROMAN] ? co->Era[SROMAN][3] : -1.0) << " " << sigpiER_ns << " "
                      << (okBs ? co->Erb[SJOINT][0] : -1.0) << " " << sigtetE_ns << " "
                      << (okBs ? co->Erb[SJOINT][3] : -1.0) << " " << sigpiEb_ns << " "
                      << co->relMl[SJOINT] << " " << relMl_ns << " "
                      << co->condA[SJOINT] << " " << condA_ns << " "
                      << co->condB[SJOINT] << " " << condB_ns << " "
                      << nepLpk << " " << nepRpk << " " << wArea << " "
                      << l->Ml << " " << l->Dl << " " << s->Ds << " " << l->Vt << "\n";
                fpair.close();
            }
//          


            //cout << "** End of saving in the file *********" << save << endl;
            //cout << "icon: " << icon << "\tnlens: " << nlens << "\tnerr: " << nerr << endl;
            // Per-sightline event budget -- see RunConfig. Was hardcoded to a stub
            // 20/5/1.0 with the production 850/150/2.0 commented out beside it.
            } while ((icon < cfg.iconTarget or nlens < cfg.nlensTarget or nerr < cfg.nerrTarget)
                     and nsim < cfg.maxDraws);

            // Did the sightline actually meet its budget, or did the cap stop it? The
            // distinction matters: a capped sightline's Poisson precision is whatever it
            // reached, not what was asked for, and averaging it in as though it were a
            // full sample would understate the error bars.
            const bool budgetMet = (icon  >= cfg.iconTarget and
                                    nlens >= cfg.nlensTarget and
                                    nerr  >= cfg.nerrTarget);
            if (!budgetMet) {
                nCapped += 1;
                cout << "  CAP: stopped at nsim = " << nsim << " with icon = " << icon
                     << "/" << cfg.iconTarget << ", nlens = " << nlens << "/" << cfg.nlensTarget
                     << ", nerr = " << nerr << "/" << cfg.nerrTarget << endl;
            }

            // Some epochs existed but nothing survived to be aggregated. Same reasoning as
            // the no-coverage skip above: the CHECK block below requires at least one
            // detected AND one characterised event, so this must not fall through.
            if (nlens < 1 or nerr <= 0.0) {
                nSkipBarren += 1;
                areaBarren  += wArea;
                cout << "  BARREN: " << nlens << " detected, " << nerr
                     << " characterised -- nothing to aggregate at this sightline" << endl;
                continue;
            }
            nAggregated    += 1;
            areaAggregated += wArea;
   
            for (int i = 0; i <= GG; ++i) {
                l->NstE[i] += l->nstE[i];
                l->NdtE[i] += l->ndtE[i];
                l->ndtE[i]  = double(l->ndtE[i] / (l->nstE[i] + eps)); // [0.0, 1.0]

                fil2 << std::fixed << std::setprecision(4)
                     << l->tEs[i]  << " "
                     << double(l->NdtE[i]  * 100.0 / (l->NstE[i] + eps)) << " "
                     << l->Mls[i]  << " "
                     << double(l->NdMl[i]  * 100.0 / (l->NsMl[i] + eps)) << " "
                     << l->pis[i]  << " "
                     << double(l->Ndpi[i]  * 100.0 / (l->Nspi[i] + eps)) << " "
                     << l->u0s[i]  << " "
                     << double(l->Ndu0[i]  * 100.0 / (l->Nsu0[i] + eps)) << " "
                     << l->mbs[i]  << " "
                     << double(l->Ndmb[i]  * 100.0 / (l->Nsmb[i] + eps)) << " "
                     << l->fbs[i]  << " "
                     << double(l->Ndfb[i]  * 100.0 / (l->Nsfb[i] + eps)) << " "
                     << l->mus[i]  << " "
                     << double(l->Ndmu[i]  * 100.0 / (l->Nsmu[i] + eps)) << " "
                     // The + eps every other column has, and these two lacked: before the
                     // counters above existed this was 0/0 and printed "-nan" on every row.
                     << double(l->Nhalo[1] * 100.0 / (l->Nhalo[0] + eps)) << " "
                     << double(l->Nself[1] * 100.0 / (l->Nself[0] + eps)) << "\n";

                fil2b << std::fixed  << std::setprecision(1)
                      << l->NdtE[i]  << " " << l->NstE[i]  << " "
                      << l->NdMl[i]  << " " << l->NsMl[i]  << " "
                      << l->Ndpi[i]  << " " << l->Nspi[i]  << " "
                      << l->Ndu0[i]  << " " << l->Nsu0[i]  << " "
                      << l->Ndmb[i]  << " " << l->Nsmb[i]  << " "
                      << l->Ndfb[i]  << " " << l->Nsfb[i]  << " "
                      << l->Ndmu[i]  << " " << l->Nsmu[i]  << " "
                      << l->Nhalo[1] << " " << l->Nhalo[0] << " "
                      << l->Nself[1] << " " << l->Nself[0] << "\n";
            }
   
            s->nstart = 0.0;
            for (int i = 0; i < Num; ++i) {
                test = double(s->nsdet[i]   / (s->nssim[i] + eps));
                s->nstart  += s->Rostari[i] * (s->Nstart   / s->Rostart) * test;
            } //Number/deg^{2}

            for (int i = 0; i < 2; ++i){
                tE[i]   = 0.0, RE[i]  = 0.0, piE[i]  = 0.0, tetE[i]  = 0.0; Vt[i]    = 0.0;
                u0[i]   = 0.0, Ml[i]  = 0.0, opd[i]  = 0.0; Dl[i]    = 0.0, Ds[i]    = 0.0;
                vl[i]   = 0.0, vs[i]  = 0.0, mbs[i]  = 0.0, fb[i]    = 0.0; numd[i]  = 0.0;
                fwhm[i] = 0.0; vsn[i] = 0.0; DelT[i] = 0.0; Struc[i] = 0.0; murel[i] = 0.0;
                Map[i]  = 0.0; nbl[i] = 0.0; Ext[i]  = 0.0;
            }

            EffiL = 0.0; EFF   = 0.0; Gamma  = 0.0; Neven = 0.0; EffiD = 0.0;
            ErtE  = 0.0; ErpiE = 0.0; ErtetE = 0.0; Erml  = 0.0; Erfb  = 0.0; nErAvg = 0.0;
            Erdl  = 0.0; Ermul = 0.0; Ermus  = 0.0; Eru0  = 0.0;
  
    int tempStruc;
    for (const auto& r : records) {
//        counter = r.counter;
        flagL = r.flagL;
        l->tE = r.tE;   l->RE = r.RE;   l->piE = r.piE;  l->tetE = r.tetE;
        l->Vt = r.Vt;   l->u0 = r.u0;   l->Ml  = r.Ml;
        s->opt = r.opt; l->Dl = r.Dl;   s->Ds  = r.Ds;   l->vl = r.vl; s->vs = r.vs;
        mbase = r.mbase; fblend = r.fblend; gg = r.gg; tempStruc = r.struc;
        s->FWHM = r.FWHM; vsave = r.vsave; l->DeltaT = r.DeltaT; l->murel = r.murel;
        co->resu[0]=r.resu0;  co->resu[1]=r.resu1;   co->resu[2]=r.resu2;
        co->resu[3]=r.resu3;  co->resu[5]=r.resu5;   co->resu[9]=r.resu9;
        co->resu[10]=r.resu10; co->resu[13]=r.resu13; co->resu[14]=r.resu14;
        s->Map[2]=r.Map2; s->nsbl[2]=r.nsbl2; co->flagi=r.flagi; s->Ai[2]=r.Ai2;
        // Must be replayed too: this loop runs after the whole field, so co->okA
        // otherwise holds whatever the LAST FisherM call left, not this event's.
        co->okA[SJOINT] = r.okJoint;

        l->struc = static_cast<GalacticComponent>(tempStruc);

        // ------------------ Update cumulative sums ------------------
        Struc[0] += 1.0;
        tE[0]    += l->tE;    RE[0]  += l->RE;     piE[0]  += l->piE;
        tetE[0]  += l->tetE;  Vt[0]  += l->Vt;     u0[0]   += l->u0;
        Ml[0]    += l->Ml;    opd[0] += s->opt;    Dl[0]   += l->Dl;
        Ds[0]    += s->Ds;    vl[0]  += l->vl;     vs[0]   += s->vs;
        mbs[0]   += mbase;    fb[0]  += fblend;    numd[0] += 1.0;
        fwhm[0]  += s->FWHM;  vsn[0] += vsave;     DelT[0] += l->DeltaT;
        murel[0] += l->murel; Map[0] += s->Map[2]; nbl[0]  += s->nsbl[2];
        Ext[0]   += s->Ai[2];

        // ------------------ Flagged case ---------------------------
        if (flagL > 0) {
            Struc[1] += 1.0;
            tE[1]    += l->tE;    RE[1]  += l->RE;     piE[1]  += l->piE;
            tetE[1]  += l->tetE;  Vt[1]  += l->Vt;     u0[1]   += l->u0;
            Ml[1]    += l->Ml;    opd[1] += s->opt;    Dl[1]   += l->Dl;
            Ds[1]    += s->Ds;    vl[1]  += l->vl;     vs[1]   += s->vs;
            mbs[1]   += mbase;    fb[1]  += fblend;    numd[1] += 1.0;
            fwhm[1]  += s->FWHM;  vsn[1] += vsave;     DelT[1] += l->DeltaT;
            murel[1] += l->murel; Map[1] += s->Map[2]; nbl[1]  += s->nsbl[2];
            Ext[1]   += s->Ai[2];

            EFF += static_cast<double>(l->ndtE[gg] / (l->tE / year)); // 1/years

            // Average only over events the joint fit could actually characterize. Step C4
            // reports sigma = -1 as an explicit "not characterizable" sentinel, and ErrorCal
            // divides it by the parameter value, so an unguarded sum pulls in large negative
            // numbers: one event in this field contributed resu[3] = -697, dragging the mean
            // fractional piE error negative and tripping CHECK(ErpiE > 0.0).
            //
            // Skipping those events is not the selection bias DEVIATIONS entry 8 warns about.
            // That rule is about never dropping an event from the joint-vs-single RATIO
            // statistics, where a missing single-survey sigma is itself the result. Here we are
            // forming a mean precision, and an event with no measurement has no precision to
            // average -- including it would be averaging a sentinel. The count of events the
            // mean is actually over is tracked separately so the denominator is honest.
            //
            // okA[SJOINT] is NOT sufficient on its own. It says the joint photometric
            // matrix inverted -- not that every parameter was in the fit. Since the joint
            // refactor gave each survey partition its own active parameter subset
            // (activePhotParams in Bulge.h), fb0 and mbs0 enter the joint fit only when the
            // event has Rubin epochs, and fb1/mbs1 only when it has Roman ones. An event
            // detected by Roman with no Rubin data therefore has a perfectly valid joint
            // fit in which Era[2] is still the -1.0 sentinel, and ErrorCal divides that by
            // fb0 regardless: one such event contributed resu[2] = -8499 at l=0.881,
            // b=-0.94 and dragged a 102-event mean to -83, tripping CHECK(Erfb > 0.0).
            //
            // So test the values themselves. An event is averaged only if all nine are real
            // measurements, which keeps every mean over the same event set and keeps
            // nErAvg meaningful as a single denominator. The cost is small and measured:
            // across 48,959 characterised events in the 2026-08-29 run exactly ONE was
            // excluded by this, and no column other than resu[2] was ever negative.
            const bool allMeasured = (co->flagi > 0 and co->okA[SJOINT]
                                      and co->resu[0]  >= 0.0 and co->resu[1]  >= 0.0
                                      and co->resu[2]  >= 0.0 and co->resu[3]  >= 0.0
                                      and co->resu[5]  >= 0.0 and co->resu[9]  >= 0.0
                                      and co->resu[10] >= 0.0 and co->resu[13] >= 0.0
                                      and co->resu[14] >= 0.0);
            if (allMeasured) {
                Eru0  += co->resu[0];  ErtE   += co->resu[1];  Erfb  += co->resu[2];
                ErpiE += co->resu[3];  ErtetE += co->resu[5];  Erml  += co->resu[9];
                Erdl  += co->resu[10]; Ermul  += co->resu[13]; Ermus += co->resu[14];
                nErAvg += 1.0;
            }
        }
    }

    for (int i = 0; i < 2; ++i) { // What is the purpose of this block?
        tE[i]    = double(tE[i]            / (numd[i] + eps)) / year; //[years]  
        RE[i]    = double(RE[i]            / (numd[i] + eps));//[AU]  
        piE[i]   = double(piE[i]           / (numd[i] + eps));//[]     
        tetE[i]  = double(tetE[i]          / (numd[i] + eps));//[mas]  
        Vt[i]    = double(Vt[i]            / (numd[i] + eps));//[km/s]  
        u0[i]    = double(u0[i]            / (numd[i] + eps));//[]  
        Ml[i]    = double(Ml[i]            / (numd[i] + eps));//[Msun]  
        opd[i]   = double(opd[i]           / (numd[i] + eps)) * u0m * u0m;//[] x10^{6}
        Dl[i]    = double(Dl[i]            / (numd[i] + eps));//[kpc]  
        Ds[i]    = double(Ds[i]            / (numd[i] + eps));//[kpc]  
        vl[i]    = double(vl[i]            / (numd[i] + eps));//[km/s]  
        vs[i]    = double(vs[i]            / (numd[i] + eps));//[km/s]  
        mbs[i]   = double(mbs[i]           / (numd[i] + eps));//[mag]  
        fb[i]    = double(fb[i]            / (numd[i] + eps));//[]
        fwhm[i]  = double(fwhm[i]          / (numd[i] + eps));//[years]
        vsn[i]   = double(vsn[i]           / (numd[i] + eps));//[]
        DelT[i]  = double(DelT[i]          / (numd[i] + eps));//[]
        murel[i] = double(murel[i]         / (numd[i] + eps));//[mas/year]
        Map[i]   = double(Map[i]           / (numd[i] + eps));//[mag]
        nbl[i]   = double(nbl[i]           / (numd[i] + eps));//[]
        Ext[i]   = double(Ext[i]           / (numd[i] + eps));//[mag]
        Struc[i] = double(Struc[i] * 100.0 / (numd[i] + eps));
    }

    EFF = double(EFF / (numd[1] + eps));//1/[years]
    Gamma = double(2.0 / M_PI * opd[0] * 1.0e-6 * EFF) / u0m; // 1/[star*year]  
    EffiL = double(numd[1] * 100.0 / (numd[0] + eps)); // probability of detecting lensing  
    EffiD = double(icon * 100.0 / (nsim    + eps)); // % of drawn stars that are visible (Deviation 78; was numd[0], i.e. 100% by construction)
    Neven = double(s->nstart * Gamma * 10.0);//deg^{-2}
   
    Eru0   = double(Eru0   / (nErAvg + eps));  
    Erfb   = double(Erfb   / (nErAvg + eps));  
    ErtE   = double(ErtE   / (nErAvg + eps));  
    ErpiE  = double(ErpiE  / (nErAvg + eps));  
    ErtetE = double(ErtetE / (nErAvg + eps));  
    Erml   = double(Erml   / (nErAvg + eps));  
    Erdl   = double(Erdl   / (nErAvg + eps));  
    Ermul  = double(Ermul  / (nErAvg + eps));  
    Ermus  = double(Ermus  / (nErAvg + eps));  
   
    //Maps
    fil3 << std::fixed << std::setprecision(6)
         << tE[0]    << " " << tE[1]    << " "
         << RE[0]    << " " << RE[1]    << " "
         << piE[0]   << " " << piE[1]   << " "
         << tetE[0]  << " " << tetE[1]  << " "
         << Vt[0]    << " " << Vt[1]    << " "
         << u0[0]    << " " << u0[1]    << " "
         << Ml[0]    << " " << Ml[1]    << " "
         << opd[0]   << " " << opd[1]   << " "
         << Dl[0]    << " " << Dl[1]    << " "
         << Ds[0]    << " " << Ds[1]    << " "
         << vl[0]    << " " << vl[1]    << " "
         << vs[0]    << " " << vs[1]    << " "
         << mbs[0]   << " " << mbs[1]   << " "
         << fb[0]    << " " << fb[1]    << " "
         << fwhm[0]  << " " << fwhm[1]  << " "
         << vsn[0]   << " " << vsn[1]   << " "
         << DelT[0]  << " " << DelT[1]  << " "
         << Struc[0] << " " << Struc[1] << " "
         << murel[0] << " " << murel[1] << " "
         << Map[0]   << " " << Map[1]   << " "
         << nbl[0]   << " " << nbl[1]   << " "
         << Ext[0]   << " " << Ext[1]   << " ";

    fil3 << std::setprecision(8)
         << EffiD      << " " << EffiL        << " "
         << std::log10(EFF) << " " << std::log10(Gamma) << " " << std::log10(Neven) << " "
         << Eru0       << " " << ErtE         << " " << Erfb         << " " << ErpiE << " " << ErtetE << " "
         << Erml       << " " << Erdl         << " " << Ermul        << " " << Ermus << " ";

    fil3 << std::setprecision(1)
         << nsim              << " " << numd[0]          << " " << numd[1] << " "
         << nerr              << " " << nri              << " " << nde     << " "
         << std::log10(s->Rostart) << " " << std::log10(s->Nstart) << " " << std::log10(s->nstart)
         // Step E1. Three columns appended, in this order:
         //   w_area  deg^2 of sky this sightline stands for -- no longer a run-wide constant
         //   lon,lat where it is. The map file had NO position column at all, so a row in it
         //           could not be tied to the events it produced, and the draw count `nsim`
         //           it records -- the denominator any pooled yield needs -- was unreachable
         //           from the event table. With these, an event joins its sightline on
         //           (lon, lat) and the correct pooled weight, w_area/nsim, is computable.
         << " " << std::setprecision(8) << wArea
         << " " << std::setprecision(6) << s->lon << " " << s->lat
         << "\n";

    // Flush the per-sightline outputs now rather than when the stream is destroyed. Every
    // production pause so far has been a kill, and a kill discards whatever is still
    // buffered: the 2026-09-06 chunk-1 stop lost the map rows of six completed sightlines
    // and left a half-written seventh, onto which the resuming run's first row was then
    // appended -- one 122-field line that made the whole file unreadable until the Python
    // reader learned to skip it (Deviation 41). The same stop cost EfLMC5/EfLMC5B the same
    // six blocks, which is why fil2/fil2b are flushed here too. A sightline costs minutes
    // of CPU, so three flushes per sightline are free, and what reaches disk is then what
    // the log says was finished.
    fil3.flush();
    fil2.flush();
    fil2b.flush();

    cout << "nsim:  "  << nsim    << "\t Ndetected:  " << icon    << "\t Nlensing:  " << nlens << "\t NError:  " << nerr << endl;
    cout << "Detection classes:";
    for (int c = 0; c < NDETCLASS; ++c)
        cout << "  " << detClassName(c) << ": " << nDetClass[c];
    cout << endl;
    if (nDetClass[DET_ANOMALY] > 0) {
        cout << "  WARNING: " << nDetClass[DET_ANOMALY] << " event(s) detected by one telescope "
             << "but NOT by the joint test. Adding data cannot destroy signal, so this is a "
             << "threshold inconsistency -- see DetClass in include/fisher/fisher.h." << endl;
    }
    cout << "numd0:  " << numd[0] << "\t numd1:  "     << numd[1] << endl;
    cout << "EFF:  "   << EFF     << "\t Gamma:  "     << Gamma   << "\t Neven:  "    << Neven << endl;
    cout << "l.tE:  "  << l->tE   << "\t l.Ml:  "      << l->Ml   << endl;
   
    CHECK(EFF > 0.0);
    CHECK(Gamma > 0.0);
    CHECK(EffiL > 0.0);
    CHECK(EffiD > 0.0);
    CHECK(Neven > 0.0);
    
    CHECK(DelT[0] > 0.0);
    CHECK(DelT[1] > 0.0);
    CHECK(vsn[0] > 0.0);
    CHECK(vsn[1] > 0.0);
    CHECK(fwhm[0] > 0.0);
    CHECK(fwhm[1] > 0.0);
    
    CHECK(fb[0] > 0.0);
    CHECK(fb[0] <= 1.0);
    CHECK(fb[1] > 0.0);
    CHECK(fb[1] <= 1.0);
    
    CHECK(mbs[0] > 0.0);
    CHECK(mbs[1] > 0.0);
    
    CHECK(vs[0] > 0.0);
    CHECK(vs[1] > 0.0);
    CHECK(vl[0] > 0.0);
    CHECK(vl[1] > 0.0);
    
    CHECK(Ds[0] > 0.0);
    CHECK(Ds[1] > 0.0);
    CHECK(Dl[0] > 0.0);
    CHECK(Dl[1] >= 0.0);
    
    CHECK(piE[1] >= 0.0);
    CHECK(tetE[1] > 0.0);
    
    CHECK(murel[0] > 0.0);
    CHECK(murel[1] > 0.0);
    
    if (nErAvg > 0.0) CHECK(ErtE > 0.0);
    if (nErAvg > 0.0) CHECK(ErpiE > 0.0);
    if (nErAvg > 0.0) CHECK(ErtetE > 0.0);
    if (nErAvg > 0.0) CHECK(Erml > 0.0);
    if (nErAvg > 0.0) CHECK(Erdl > 0.0);
    if (nErAvg > 0.0) CHECK(Ermul > 0.0);
    if (nErAvg > 0.0) CHECK(Ermus > 0.0);
    if (nErAvg > 0.0) CHECK(Eru0 > 0.0);
    if (nErAvg > 0.0) CHECK(Erfb > 0.0);
    
    CHECK(nerr != 0.0);
    CHECK(numd[0] != 0.0);
    CHECK(numd[1] != 0.0);
    CHECK(nsim != 0.0);
    
    // NOT an equality. `icon` counts stars that were OBSERVABLE (flagf > 0 and ndw > 2);
    // `numd[0]` counts every record pushed, and the push site sits OUTSIDE that gate, so a
    // star that was drawn but never observable still gets a record. The two are equal only
    // where every draw is observable, which is true in the dense stub patch that every run
    // before commit 81a6b04 used and false as soon as the scan reaches sparse sky -- at
    // l=-3.499, b=-1.98 it is 5 observable out of 15 drawn.
    //
    // Loosening this assertion changes no computed value. It does expose a real question
    // about what the per-sightline denominators mean -- EffiD = numd[0]/nsim is 100% by
    // construction, and the [0] means average over drawn rather than observed stars. That
    // is a science decision, recorded in OPEN_ITEMS.md, not something to change here.
    CHECK(icon <= numd[0]);
    CHECK(numd[1] == nlens);
     
    cout << "==============================================================" << endl;  
  
   }}//end of the sightline loop (was: right_accention and declinaton)

    // ---------------------------------------------------------------------------------------
    // Run-wide detection statistics.
    //
    // Counts are Poisson, so each is quoted with its sqrt(N) uncertainty -- the joint-only
    // class is expected to be rare, and a rare count without an error bar cannot be argued
    // from. Percentages are of all simulated events, which is the denominator that makes them
    // comparable between runs of different length.
    // ---------------------------------------------------------------------------------------
    cout << "\n================ RUN TOTALS ================" << endl;
    cout << "dchiL bookkeeping: " << nDchiMismatch << " event(s) with joint != Rubin + Roman"
         << (nDchiMismatch ? "  <-- a BUG, see DCHI_MISMATCH lines" : " (as it must be)") << endl;
    // Sightline accounting. The first THREE partition the grid and must sum to its size;
    // if they do not, a sightline left the loop by a path nobody wrote down. `nCapped` is
    // not part of that partition -- it overlaps both aggregated and barren, since hitting
    // the draw cap says how a sightline stopped, not what it yielded.
    cout << "Sightlines: " << nAggregated << " aggregated, "
         << nSkipNoCoverage << " skipped (no coverage), "
         << nSkipBarren << " skipped (barren), "
         << nCapped << " hit --maxdraws" << endl;
    // Summed weights, not count x cell area: under stratification the sightlines stand for
    // different pieces of sky and the two stopped being the same number (Step E1).
    cout << "  scanned area " << (areaAggregated + areaNoCoverage + areaBarren)
         << " deg^2, of which " << areaAggregated
         << " deg^2 produced events" << endl;
    if (nCapped > 0)
        cout << "  WARNING: " << nCapped << " sightline(s) stopped on the draw cap with their "
             << "budget unmet. Their Poisson precision is lower than requested -- raise "
             << "--maxdraws or accept the larger error bars." << endl;
    {
        // Appended, not written with the startup block: these counts are only known now,
        // and the area weighting downstream depends on them.
        std::ofstream fprov(std::string(PATH_OUT_DIR) + "run_provenance.txt", std::ios::app);
        if (fprov) {
            fprov << "# ---- sightline outcome (written at end of run) ----\n"
                  << "# sightlines_aggregated   " << nAggregated << "\n"
                  << "# sightlines_no_coverage  " << nSkipNoCoverage << "\n"
                  << "# sightlines_barren       " << nSkipBarren << "\n"
                  << "# sightlines_capped       " << nCapped << "\n"
                  << "# area_with_events_deg2   " << areaAggregated << "\n"
                  << "# area_no_coverage_deg2   " << areaNoCoverage << "\n"
                  << "# area_barren_deg2        " << areaBarren << "\n";
        }
    }
    // Draws that never produced a light curve never reach the detection test, so they are not
    // in NDetClassTot at all. Deriving DET_NONE by subtraction rather than counting it there
    // keeps the classes summing to the number of simulated events, which is what makes the
    // percentages below meaningful as detection rates.
    {
        long det = 0;
        for (int c = 1; c < DET_ANOMALY; ++c) det += NDetClassTot[c];
        NDetClassTot[DET_NONE] = nSimTot - det;
    }
    cout << "simulated events: " << nSimTot << endl;
    for (int c = 0; c < NDETCLASS; ++c) {
        const long n = NDetClassTot[c];
        cout << "  " << std::left << std::setw(26) << detClassName(c) << std::right
             << std::setw(9) << n << " +/- " << std::setw(7) << std::fixed
             << std::setprecision(1) << std::sqrt(double(n))
             << "   (" << std::setprecision(4)
             << (nSimTot > 0 ? 100.0 * double(n) / double(nSimTot) : 0.0) << "%)" << endl;
    }
    if (NDetClassTot[DET_ANOMALY] > 0) {
        cout << "  WARNING: DET_ANOMALY counts events where the RAW joint test contradicted a "
             << "single-survey\n           detection. Since Step H7 all three tests share one "
             << "fixed bar (--dchi-det " << cfg.dchiDet << "), and\n           because chi2 "
             << "accumulates over both instruments this count should be ZERO by\n           "
             << "construction. A non-zero value means that construction is broken -- most "
             << "likely the\n           signed lensing statistic has picked up a sign "
             << "convention it should not have.\n           detJ is still monotone downstream, "
             << "so no output table is wrong, but investigate." << endl;
    }

    // Per-tE breakdown, printed only for bins that contain something, since most are empty.
    cout << "\n--- detections by tE bin ---" << endl;
    cout << std::left << std::setw(12) << "tE [d]" << std::right;
    for (int c = 1; c < DET_ANOMALY; ++c) cout << std::setw(14) << detClassName(c);
    cout << endl;
    for (int i = 0; i <= GG; ++i) {
        long tot = 0;
        for (int c = 1; c < DET_ANOMALY; ++c) tot += NDetClassTE[i][c];
        if (tot == 0) continue;
        cout << std::left << std::setw(12) << std::setprecision(2) << l->tEs[i] << std::right;
        for (int c = 1; c < DET_ANOMALY; ++c) cout << std::setw(14) << NDetClassTE[i][c];
        cout << endl;
    }
    cout << "============================================" << endl;

   //fclose(filh);  
   return(0);
}
