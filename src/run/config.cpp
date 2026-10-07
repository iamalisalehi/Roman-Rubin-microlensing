// Usage text for the command-line flags.
#include "run/config.h"
#include "surveys/footprints.h"

void printUsage(const char* prog) {
    std::cout
        << "usage: " << prog << " [options]\n"
        << "  --stride N     sightline grid step, in units of dd=" << dd << " deg (default 10)\n"
        << "                 grid spacing is N*dd deg. If the footprint grid (see\n"
        << "                 --stride-roman) would be coarser than one Roman detector, it\n"
        << "                 is refined automatically when --stride-roman is not given.\n"
        << "  --seed S       base random seed (default " << seed << "); every sightline is\n"
        << "                 re-seeded from (S, its scan index), so runs are reproducible\n"
        << "                 sightline by sightline, however the scan is split\n"
        << "  --end-index N  stop before sightline N (with --start-index: one chunk of a split\n"
        << "                 scan; chunks [0,a), [a,b), ... together equal the full run)\n"
        << "  --start-index N   skip the first N sightlines of the scan and resume there.\n"
        << "                 The scan vector is deterministic for a given --stride/\n"
        << "                 --stride-roman/--stub, and the output files are opened in\n"
        << "                 append mode, so a run interrupted after N sightlines is\n"
        << "                 continued exactly by re-running with the same flags plus\n"
        << "                 --start-index N, where N counts sightlines the scan ENTERED:\n"
        << "                     grep -c 'NEW STEP' <the interrupted run's stdout log>\n"
        << "                 NOT the line count of MapLMC5.dat. A sightline with no\n"
        << "                 coverage, or a barren one, is entered and counted but never\n"
        << "                 writes a map row, so the map file undercounts -- by 86 of 710\n"
        << "                 on the 2026-09-05 v2 run. Resuming at the map count would\n"
        << "                 re-simulate that gap and append duplicate rows to test5.dat.\n"
        << "                 If the run resumed an earlier one, add its --start-index too.\n"
        << "  --pair-satellite\n"
        << "                 characterise each detected event twice, with Roman at L2 and\n"
        << "                 with the offset zeroed, writing both forecasts to h3_pair.dat.\n"
        << "                 Step H3's experiment, done on one event at a time because two\n"
        << "                 runs cannot be compared event-for-event. Costs about one extra\n"
        << "                 Fisher call per detection. Incompatible with\n"
        << "                 --no-satellite-parallax, which leaves nothing to compare.\n"
        << "  --dchi-det X   delta-chi2 a lensing model must beat a flat baseline by for\n"
        << "                 an event to count as detected (default "
        << DCHI_DET_DEFAULT << ", Penny+2019).\n"
        << "                 A FIXED bar, applied identically to the Rubin-only, Roman-only\n"
        << "                 and joint tests, which is what makes the joint test monotone:\n"
        << "                 chi2 accumulates over both instruments, so dchi = dchi_L +\n"
        << "                 dchi_R and either survey clearing the bar alone forces the sum\n"
        << "                 over it. Every yield is conditioned on this number; it is\n"
        << "                 recorded in run_provenance.txt.\n"
        << "  --stride-roman N  sightline grid step inside Roman's footprint, same units\n"
        << "                 (default: same as --stride, i.e. an unstratified scan). Must\n"
        << "                 divide --stride. Sightlines outside the footprint keep the\n"
        << "                 coarse step and carry the sky area they stand for.\n"
        << "  --events N     per-sightline detected-event target, icon (default 850)\n"
        << "  --lenses N     per-sightline characterised-event target, nlens (default 150)\n"
        << "  --nerr X       per-sightline Fisher-error target (default 2.0)\n"
        << "  --maxdraws X   per-sightline cap on drawn stars (default 5e4)\n"
        << "  --stub         scan a 0.1x0.1 deg test patch inside GBTDS field 3 (both rolls)\n"
        << "  --population N which lens population to simulate (default 'bulge').\n"
        << "                 Every output file is named from the population's tag, so two\n"
        << "                 populations never overwrite each other: 'bulge' writes\n"
        << "                 test5.dat as before, 'bh' writes testbh.dat, 'ns' testns.dat.\n"
        << "                 bulge = Kroupa IMF + remnants, 0.01-30 Msun\n"
        << "                 besancon = ordinary stars, mass and light of a random member of\n"
        << "                            the lens's component Besancon list (no BH/NS/BD)\n"
        << "                 bh    = flat in log M, 3-1000 Msun\n"
        << "                 ns    = neutron stars, Gaussian about 1.35 Msun\n"
        << "                 macho-uniform/-m05/-m1/-m2 = the legacy 3-5000 Msun options\n"
        << "  --dry-run      build the sightline grid, report the strata and the sky-area\n"
        << "                 weights, then exit without drawing any stars\n"
        << "  --no-satellite-parallax   put Roman at the centre of the Earth, removing the\n"
        << "                 Earth-L2 spatial baseline while leaving the timing alone. The\n"
        << "                 'off' run of the satellite-parallax experiment (PHASE_H_PLAN H3)\n"
        << "  --dump-samples F  write full light curves and astrometric tracks for a few\n"
        << "                 illustrative events, as described by spec file F. Off by\n"
        << "                 default. Consumes no RNG, so the set of simulated events is\n"
        << "                 identical with and without it. See the spec-file format in\n"
        << "                 the Step S1 block of include/run/sample_dump.h.\n"
        << "  --help         this message\n";
}

bool parseCommandLine(int argc, char** argv, RunConfig& cfg, int& exitCode) {
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
        else if (arg == "--help")  { printUsage(argv[0]); exitCode = 0; return false; }
        else {
            std::cerr << "ERROR: unknown option '" << arg << "'\n";
            printUsage(argv[0]);
            exitCode = 2; return false;
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
        exitCode = 1; return false;
    }
    if (not (cfg.dchiDet > 0.0) or not std::isfinite(cfg.dchiDet)) {
        std::cerr << "ERROR: --dchi-det (" << cfg.dchiDet << ") must be finite and positive. "
                  << "It is a delta-chi2 detection bar; a non-positive bar would declare every "
                  << "event detected.\n";
        exitCode = 1; return false;
    }
    if (cfg.stride < 1) {
        std::cerr << "ERROR: --stride must be >= 1\n";
        exitCode = 2; return false;
    }
    if (cfg.iconTarget < 1 or cfg.nlensTarget < 0 or cfg.nerrTarget < 0.0) {
        std::cerr << "ERROR: --events must be >= 1, --lenses and --nerr >= 0\n";
        exitCode = 2; return false;
    }
    // A cap below the event budget would stop every sightline early, which is not a cap
    // but a silent redefinition of the budget.
    if (cfg.endIndex >= 0 and cfg.endIndex <= cfg.startIndex) {
        std::cerr << "ERROR: --end-index (" << cfg.endIndex << ") must exceed --start-index ("
                  << cfg.startIndex << ").\n";
        exitCode = 2; return false;
    }
    if (cfg.startIndex < 0) {
        std::cerr << "ERROR: --start-index (" << cfg.startIndex << ") cannot be negative.\n";
        exitCode = 1; return false;
    }
    if (cfg.maxDraws < double(cfg.iconTarget)) {
        std::cerr << "ERROR: --maxdraws (" << cfg.maxDraws << ") is below --events ("
                  << cfg.iconTarget << "); no sightline could reach its budget\n";
        exitCode = 2; return false;
    }
    return true;
}

int resolveGridSteps(RunConfig& cfg, const GbtdsLayout& gl, GridSteps& steps) {
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
    steps.gridStep = gridStep;
    steps.kSub     = kSub;
    steps.fineStep = fineStep;
    return 0;
}
