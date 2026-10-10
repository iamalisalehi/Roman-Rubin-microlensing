// The run's output files, the per-event table header and the provenance block.
#include "run/outputs.h"
#include <sstream>   // provenance block

// The Makefile makes this object depend on the git stamp: GIT_COMMIT is printed below.
#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif

// Column names of the per-event table, in exactly the order the `filg_in <<` block writes
// them. Anything added to the row must be appended both here and there; check with
// `head -1 | wc -w` against a data line's `wc -w`.
//
// Two indexing systems that look alike: magb_*/blend_* are per FILTER (ugrizy, F146);
// mbs0/fb0 and mbs1/fb1 are per TELESCOPE (0 = Rubin, 1 = Roman).
const char* eventTableHeader()
{
    return
        "# icon FFG0 tE RE_AU piE tetE Vt u0 Ml opt_1e6 Dl Ds vl vs mbs0 fb0 gg struc "
        "FWHM_yr vsave_mus DeltaT_errA murel_yr "
        "rel_u0 rel_tE rel_fb0 rel_piE rel_tetE rel_Ml rel_Dl rel_mul rel_mus "
        "Map_r nsbl_r flagi Ai_r "
        "ndw_L ndw_R detL detR detJ okA_J okA_L okA_R "
        "sigtE_J sigtE_L sigtE_R sigpiE_J sigpiE_L sigpiE_R sigtetE_J sigtetE_L sigtetE_R "
        "detCls synClass condA_J condA_L condA_R "
        "t0 xi lon lat mbs1 fb1 "
        "magb_u magb_g magb_r magb_i magb_z magb_y magb_F146 "
        "blend_u blend_g blend_r blend_i blend_z blend_y blend_F146 "
        "relMl_J relMl_L relMl_R okB_J okB_L okB_R condB_J condB_L condB_R "
        "dt_edge t0zone w_area du_sat nepL_pk nepR_pk "
        // Epoch counts at which the two lensing-induced images were separately detectable and
        // separated by more than the bar named in the suffix; dsep_max is the largest
        // separation reached at such an epoch, -1 if there was none.
        "nres5_L nres20_L nresPSF_L dsep_max_L nres5_R nres20_R nresPSF_R dsep_max_R "
        // The astrometric forecast under the N (nominal) and P (pessimistic) noise variants, joint and Roman partitions; the main sigtetE_*/relMl_*/okB_* columns
        // are variant W (white). Rubin's partition is the same in all three. See AST_SIGC.
        "sigtetE_NJ sigtetE_NR sigtetE_PJ sigtetE_PR relMl_NJ relMl_NR relMl_PJ relMl_PR "
        "okB_NJ okB_NR okB_PJ okB_PR "
        // Luminous lens (1 = a main-sequence star whose light is blended) and the lens's
        // share of the baseline flux in Rubin's reference band and in F146.
        "lensLum fLens_L fLens_R "
        // The observed (Earth-frame, parallax-bent) peak time and impact parameter.
        // t0zone, dt_edge and nep_pk_* are measured from t0obs, not from t0.
        "t0obs umin_obs "
        // Distance of the unresolved neighbours' light centroid from the source [mas], Rubin's
        // reference band and F146 (0 = no neighbour in the disc).
        "blendOff_L blendOff_R "
        // Finite-source size: the source radius [Rsun], its angular radius theta* [mas], and
        // rho = theta* / thetaE, the source radius in Einstein radii.
        "Rstar thetaStar rho";
}

int openOutputs(const RunConfig& cfg, RunOutputs& o) {
    // A fresh clone has no output directories; create them (a no-op when they exist).
    {
        std::error_code ec;
        std::filesystem::create_directories(PATH_OUT_DIR, ec);
        std::filesystem::create_directories(PATH_DENSITY_DIR, ec);
    }
    // File names, named from the population's tag so a black-hole run cannot append to the
    // bulge run's table (the default population's tag is "5").
    const std::string tag(gPop->tag);
    std::string fnLDt   = std::string(PATH_OUT_DIR) + "LpLMC"  + tag +  ".dat";
    std::string fnEff   = std::string(PATH_OUT_DIR) + "EfLMC"  + tag +  ".dat";
    std::string fnEffB  = std::string(PATH_OUT_DIR) + "EfLMC"  + tag + "B.dat";
    std::string fnGam   = std::string(PATH_OUT_DIR) + "MapLMC" + tag +  ".dat";
    std::string testf   = "./test"                       + tag +  ".dat";
    std::string fnPair  = "./h3_pair.dat";   //written only under --pair-satellite

    // Every output accumulates across the scan, so a continuation (--start-index > 0) must
    // append to them; truncating would discard everything the interrupted run wrote.
    const bool resuming = (cfg.startIndex > 0);
    // A --dry-run must not touch the previous run's outputs, so it opens them in append mode
    // (and writes nothing); a fresh run truncates every output.
    const bool keepOld = resuming or cfg.dryRun;
    const std::ios::openmode accumulate =
        std::ios::out | (keepOld ? std::ios::app : std::ios::trunc);

    // LpLMC<tag>.dat is appended to by the sightline loop. Create it if missing so a population
    // whose tag has never been run can start; on a resume it is left untouched.
    { std::ofstream ensureLp(fnLDt, accumulate); }   // create; truncated on a fresh run
    std::ofstream& fil2  = o.fil2;   std::ofstream& fil2b = o.fil2b;   std::ofstream& fil3 = o.fil3;
    fil2.open(fnEff,   accumulate);
    fil2b.open(fnEffB, accumulate);
    fil3.open(fnGam, accumulate);

    // The per-event table: write the column header once, in a scope of its own, and leave
    // `filg_in` itself closed. The per-event write calls filg_in.open(), and open() on an
    // already-open ofstream sets failbit and is silently discarded. The open/append/close per
    // event is deliberate: each row reaches disk as it is produced, so an interrupted run
    // keeps what it computed.
    //
    // The header is written when, and only when, the table has no content yet. That covers
    // both uses of --start-index: continuing an interrupted run (table already holds the
    // header and earlier rows, which must be kept) and starting a fresh run at an offset in an
    // empty directory (nothing to preserve, header needed).
    std::streamoff tableBytes = -1;
    {
        std::ifstream probe(testf, std::ios::ate | std::ios::binary);
        tableBytes = probe ? static_cast<std::streamoff>(probe.tellg()) : std::streamoff(-1);
    }
    if (cfg.dryRun) {
        // dry run: leave the previous run's table alone.
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
        // Re-running the scan from the start in a directory that already holds a table.
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
    // The --pair-satellite side file. Same rule as the event table: header only when empty.
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
            // Ml, Dl, Ds and Vt are appended so older files still parse by position. The pooled
            // event-rate weight needs sqrt(Ml)*Vt*Z(Ds) per event, and pi_rel = 1/Dl - 1/Ds
            // cannot be solved for Dl and Ds from the other columns.
              << "nepL_pk nepR_pk w_area Ml Dl Ds Vt\n";
        }
    }

    // o.filg_in is opened in append mode per event (see above).

    // Sample-dump state: `dumpBuf` holds the current draw's recorded epochs and is cleared at
    // the top of every draw; `dumpSeq` numbers the events actually written out.
    SampleSpec&            dumpSpec = o.dumpSpec;
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

    if (!fil2 || !fil2b || !fil3) {
        std::cerr << "Cannot open one or more files!" << std::endl;
        return 1;
    }
    o.fnLDt = fnLDt;
    o.testf = testf;
    o.fnPair = fnPair;
    return 0;
}

int writeRunProvenance(const RunConfig& cfg, const GbtdsLayout& gl, const GridSteps& steps,
                       const SightlineGrid& grid, const RomanSchedule& sched, const extin& ex) {
    const double gridStep = steps.gridStep;
    const int    kSub     = steps.kSub;
    const double fineStep = steps.fineStep;
    const double scanReach = grid.scanReach;
    const double lonMin = grid.lonMin, lonMax = grid.lonMax, latMin = grid.latMin, latMax = grid.latMax;
    const long   nSightlines = grid.nSightlines, nSightlinesRoman = grid.nSightlinesRoman;
    const double areaFootprint = grid.areaFootprint, areaOutside = grid.areaOutside;
    const double areaScanned = grid.areaScanned;
    const auto&  areaOnDetector = grid.areaOnDetector;
    const auto&  classGrid = grid.classGrid;
    const auto&  classExact = grid.classExact;
    const auto&  romanFields = grid.romanFields;
    const size_t nFieldsCovered = grid.nFieldsCovered;

    // Provenance, written before any science output so an interrupted run still records what
    // it was. Neven is a surface density in deg^-2 and survey-wide yields are summed downstream,
    // each sightline weighted by the sky area it stands for: (N*dd)^2 at stride N, not dd^2.
    // With the stratified grid that weight differs between sightlines, so area_per_sightline is
    // meaningful only when stratified = 0; otherwise use the per-row `w_area` column.
    const double areaPerSightline = gridStep * gridStep; // deg^2; unstratified runs only

    {
        std::ostringstream prov;
        prov << "# Roman+Rubin microlensing forecast -- run provenance\n"
             << "# git_commit          " << GIT_COMMIT << "\n"
             << "# built               " << __DATE__ << " " << __TIME__ << "\n"
             // The pooled event-rate weight carries a sqrt(Ml) factor that depends on the
             // mass function sampled, so the analysis layer reads the population from here.
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
             << "   # base; each sightline re-seeded from (seed, index)\n"
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
             << "# roman_noise         phot: Penny+2019 curve anchored to 5-sigma "
             << ROMAN_DEPTH5_AB << " AB (66 s); depth " << thre[6] << ", saturation "
             << satu[6] << " AB\n"
             << "# astrometric_noise   k/SNR (+) floor, SNR from the epoch's photometric error; Roman k "
             << ROMAN_AST_K << " mas, floor " << ROMAN_AST_FLOOR << " mas; Rubin kappa "
             << LSST_AST_KAPPA << " x visit FWHM_geom, floor " << LSST_AST_FLOOR << " mas; blend offset "
             << (AST_BLEND_OFFSET_FREE ? "fitted (marginalised)" : "known") << "\n"
             << "# astrometric_variants W: one frame, white; N: frame per roll + " << cfg.astSigcN
             << " mas/day; P: frame per roll + " << cfg.astSigcP << " mas/day\n"
             << "# extinction          files/ext/ext_tables.dat: " << ex.nTables << " x "
             << ex.nDist << ", k " << ex.k << " --" << ex.built << "\n"

             << "# rng_seed            " << cfg.seedBase << "\n"
             // 1 = Roman at Sun-Earth L2; 0 = Roman at the centre of the Earth, so the piE
             // forecast contains only the annual Earth-orbit parallax.
             << "# satellite_parallax  " << (cfg.noSatPar ? 0 : 1)
             << "   # 0 = Roman forced to Earth's position\n"
             << "# L2_offset_AU        " << (cfg.noSatPar ? 0.0 : L2_OFFSET_AU) << "\n"
             << "# dump_samples        " << (cfg.dumpSpec.empty() ? std::string("none")
                                                                    : cfg.dumpSpec) << "\n"
             << "# pair_satellite      " << (cfg.pairSat ? 1 : 0)
             << "   # every detection characterised at L2 AND at Earth\n"
             << "# dchi_det            " << cfg.dchiDet
             << "   # fixed detection bar; every yield is conditioned on it\n";
        if (!cfg.dryRun) {   // a dry run must not overwrite the last run's provenance
            std::ofstream fprov(std::string(PATH_OUT_DIR) + "run_provenance.txt");
            if (!fprov) {
                std::cerr << "Cannot write run_provenance.txt\n";
                return 1;
            }
            fprov << prov.str();
        }
        std::cout << prov.str() << std::flush;
    }
    return 0;
}
