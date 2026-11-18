// The epochs of an event's light curve: the visits the two surveys take on a sightline, merged into one
// time-ordered grid. Separate from observe.cpp so that the grid can be tested without any data files.
#pragma once
#include <vector>

// One epoch of the time grid, and what the two surveys do there.
struct TimeStep {
    double tim;        // [days]: the visit's own time (both visits' if they coincide); the seed point's time
    double dt;         // tim minus the previous point's tim; GRID_SEED_LEAD_DAYS at the seed point [days], > 0
    int    sq, sqR;    // the visit index (lsCt[k] / roCt[k]) taken here, -1 if that survey takes none
    bool   rubin, roman;
};

// The sightline's epochs: the sorted union of the matched visits' times, coincident times merged into one point,
// preceded by a seed point GRID_SEED_LEAD_DAYS before the first visit. Empty if ndd == nddR == 0.
// lsCt/roCt are read only at [0, ndd) / [0, nddR); times must be strictly increasing within each list (CHECK).
void buildTimeGrid(int ndd, const std::vector<int>& lsCt, const std::vector<double>& lsTim,
                   int nddR, const std::vector<int>& roCt, const std::vector<double>& roTim,
                   std::vector<TimeStep>& grid);
