// The list of sightlines to simulate, each with the sky area it stands for.
#ifndef ROMAN_RUN_SIGHTLINES_H
#define ROMAN_RUN_SIGHTLINES_H

#include "common.h"
#include "run/config.h"
#include "surveys/visits.h"
#include "surveys/footprints.h"

struct Sightline { double lon, lat, area; bool inFootprint; int col; int romanClass; };

// What the sightline build produces: the scan list and the bookkeeping that the provenance
// block, the dry-run summary and the loop need afterwards.
struct SightlineGrid {
    std::vector<Sightline>      scan;                 // the sightlines, in scan order
    std::vector<FieldPlacement> romanFields;          // distinct (centre, layout) Roman placements
    double scanReach  = 0.0;                          // scan region: within this of a field centre [deg]
    double lonMin = 0.0, lonMax = 0.0, latMin = 0.0, latMax = 0.0;
    double cellArea   = 0.0;                          // deg^2, one fine cell
    long   nSightlines = 0, nSightlinesRoman = 0;
    double areaFootprint = 0.0, areaOutside = 0.0, areaScanned = 0.0;   // deg^2
    std::array<double, GBTDS_NLAYOUT> areaOnDetector{};  // point-sampled, per roll
    std::array<double, 4> classGrid{}, classExact{};     // none/spring/autumn/both, deg^2
    size_t nFieldsCovered = 0;                        // placements with a sightline on a detector
};

// Field placements, coverage guard, the stratified grid with its areas, the area-weight
// post-stratification and the printed summary. Returns 0, or the exit code of the error printed.
int buildSightlines(const RunConfig& cfg, const GbtdsLayout& gl, const GridSteps& steps,
                    const lsst& ls, const roman& ro, SightlineGrid& grid);

// The --dry-run report: strata, areas and shares. Prints only.
void printDryRunStrata(const SightlineGrid& grid, const GridSteps& steps);

#endif // ROMAN_RUN_SIGHTLINES_H
