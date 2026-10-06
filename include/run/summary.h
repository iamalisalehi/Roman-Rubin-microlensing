// The end-of-run report: the totals the sightline loop accumulated, printed and appended to the
// provenance file.
#ifndef ROMAN_RUN_SUMMARY_H
#define ROMAN_RUN_SUMMARY_H

#include "common.h"
#include "types.h"
#include "run/config.h"
#include "fisher/fisher.h"

// Run-wide tallies, filled by main() from the loop's counters once the scan is done.
struct RunTotals {
    long   nDchiMismatch = 0;     // events where dchiL != dchiL_L + dchiL_R (should stay 0)
    int    nAggregated = 0;       // sightlines that reached the aggregation block
    int    nSkipNoCoverage = 0;   // no Rubin AND no Roman epochs
    int    nSkipBarren = 0;       // epochs, but nothing characterised
    int    nCapped = 0;           // stopped by --maxdraws with the budget unmet
    double areaAggregated = 0.0, areaNoCoverage = 0.0, areaBarren = 0.0;   // deg^2
    long   nSimTot = 0;           // simulated events
    std::array<long, NDETCLASS> NDetClassTot{};                   // detection-class counts
    std::vector<std::array<long, NDETCLASS>> NDetClassTE;         // ... per tE bin (GG + 1 bins)
};

// Print the run totals and the detection-class tables, and append the sightline outcome to
// run_provenance.txt.
void writeRunSummary(const RunConfig& cfg, const RunTotals& totals, const lens& l);

#endif // ROMAN_RUN_SUMMARY_H
