// The run's output files, the per-event table header and the provenance block.
#include "run/outputs.h"
#include <sstream>   // provenance block

// The Makefile makes this object depend on the git stamp: GIT_COMMIT is printed below.
#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif

// ---------------------------------------------------------------------------
// Column names of the per-event table, in exactly the order the `filg_in <<` block in
// main() writes them (Step D1).
//
// Written once at the top of the file so a table is self-describing: an unlabelled
// 88-column matrix is unusable six months later, and mis-numbering a column by one is
// the kind of error that produces a plausible plot of the wrong quantity. Anything
// added to the row MUST be appended both here and there, in the same place -- the
// verification step is `head -1 | wc -w` against a data line's `wc -w`.
//
// Note the two indexing systems, which look alike and are not (see CLAUDE.md):
// magb_*/blend_* are per FILTER (ugrizy, F146); mbs0/fb0 and mbs1/fb1 are per TELESCOPE
// (0 = Rubin, 1 = Roman).
// ---------------------------------------------------------------------------
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
        // Step R1. Epoch counts at which the two lensing-induced images were separately
        // detectable AND separated by more than the bar named in the suffix; dsep_max is the
        // largest separation reached at such an epoch, -1 if there was none.
        "nres5_L nres20_L nresPSF_L dsep_max_L nres5_R nres20_R nresPSF_R dsep_max_R "
        // Deviation 71. The astrometric forecast under the N (nominal) and P (pessimistic)
        // noise variants, joint and Roman partitions; the main sigtetE_*/relMl_*/okB_* columns
        // are variant W (white). Rubin's partition is the same in all three. See AST_SIGC.
        "sigtetE_NJ sigtetE_NR sigtetE_PJ sigtetE_PR relMl_NJ relMl_NR relMl_PJ relMl_PR "
        "okB_NJ okB_NR okB_PJ okB_PR "
        // Deviation 74: luminous lens (1 = a main-sequence star whose light is blended) and the
        // lens's share of the baseline flux in Rubin's reference band and in F146.
        "lensLum fLens_L fLens_R "
        // Deviation 76: the observed (Earth-frame, parallax-bent) peak time and impact parameter.
        // t0zone, dt_edge and nep_pk_* are now measured from t0obs, not from t0.
        "t0obs umin_obs";
}

int openOutputs(const RunConfig& cfg, RunOutputs& o) {
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
    std::ofstream& fil2  = o.fil2;   std::ofstream& fil2b = o.fil2b;   std::ofstream& fil3 = o.fil3;
    fil2.open(fnEff,   accumulate);
    fil2b.open(fnEffB, accumulate);
    fil3.open(fnGam, accumulate);   // was always ios::app (Deviation 78)

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

    // o.filg_in: opened in append mode per event -- see above

    // Step S1 state. `dumpBuf` holds the CURRENT draw's recorded epochs and is cleared at
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

    // Check all
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
    // column of the event table (and the w_area column of the map file). The header
    // says so, so that a downstream script cannot quietly use the wrong one.
    // ----------------------------------------------------------------------
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
             << "# extinction          files/ext/ext_tables.dat: " << ex.nTables << " x "
             << ex.nDist << ", k " << ex.k << " --" << ex.built << "\n"

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
    return 0;
}
