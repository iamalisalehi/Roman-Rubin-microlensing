// The end-of-run report. Moved out of main() verbatim.
#include "run/summary.h"

void writeRunSummary(const RunConfig& cfg, const RunTotals& totals, const lens& l) {
    const long nDchiMismatch = totals.nDchiMismatch;
    const int  nAggregated = totals.nAggregated, nSkipNoCoverage = totals.nSkipNoCoverage;
    const int  nSkipBarren = totals.nSkipBarren, nCapped = totals.nCapped;
    const double areaAggregated = totals.areaAggregated, areaNoCoverage = totals.areaNoCoverage;
    const double areaBarren = totals.areaBarren;
    const long nSimTot = totals.nSimTot;
    std::array<long, NDETCLASS> NDetClassTot = totals.NDetClassTot;   // DET_NONE is derived below
    const auto& NDetClassTE = totals.NDetClassTE;

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
        cout << std::left << std::setw(12) << std::setprecision(2) << l.tEs[i] << std::right;
        for (int c = 1; c < DET_ANOMALY; ++c) cout << std::setw(14) << NDetClassTE[i][c];
        cout << endl;
    }
    cout << "============================================" << endl;

}
