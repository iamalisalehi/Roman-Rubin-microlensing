// The merged time grid of a sightline (see include/sim/timegrid.h).
//
// The two surveys' matched visit lists are each strictly increasing in time, so the grid is their merge:
// at every step the earlier of the two next visits is taken, or both if they fall at the same instant.
// Every point is a visit's own time, so no epoch is recorded earlier or later than the sky was observed;
// the one extra point, GRID_SEED_LEAD_DAYS before the first visit, gives the source-velocity diagnostic
// a previous deflection to difference against at the first visit.
#include "sim/timegrid.h"
#include "common.h"

void buildTimeGrid(int ndd, const std::vector<int>& lsCt, const std::vector<double>& lsTim,
                   int nddR, const std::vector<int>& roCt, const std::vector<double>& roTim,
                   std::vector<TimeStep>& grid)
{
    grid.clear();
    if (ndd <= 0 and nddR <= 0) return;
    grid.reserve(std::size_t(std::max(ndd, 0)) + std::size_t(std::max(nddR, 0)) + 1);

    // The lists are read only below their counts, so an empty survey's vectors may be empty.
    const double firstRubin = ndd  > 0 ? lsTim[lsCt[0]] : 0.0;
    const double firstRoman = nddR > 0 ? roTim[roCt[0]]  : 0.0;
    const double first = (ndd > 0 and nddR > 0) ? std::min(firstRubin, firstRoman)
                                                : (ndd > 0 ? firstRubin : firstRoman);

    grid.push_back(TimeStep{first - GRID_SEED_LEAD_DAYS, GRID_SEED_LEAD_DAYS, -1, -1, false, false});

    int    gi = 0, giR = 0;                    // the next visit due in each survey's list
    double prevRubin = -1.0, prevRoman = -1.0; // the last time taken from each list, for the ordering check
    while (gi < ndd or giR < nddR) {
        const bool haveL = gi < ndd, haveR = giR < nddR;
        const double tL = haveL ? lsTim[lsCt[gi]] : 0.0;
        const double tR = haveR ? roTim[roCt[giR]] : 0.0;
        const bool takeL = haveL and (not haveR or tL <= tR);
        const bool takeR = haveR and (not haveL or tR <= tL);

        TimeStep step;
        step.rubin = takeL;
        step.roman = takeR;
        step.tim   = takeL ? tL : tR;
        step.sq    = takeL ? lsCt[gi]  : -1;
        step.sqR   = takeR ? roCt[giR] : -1;
        if (takeL) {
            CHECK(tL >= 0.0 and tL <= Tobs);
            CHECK(gi == 0 or tL > prevRubin);
            prevRubin = tL;
            ++gi;
        }
        if (takeR) {
            CHECK(tR >= 0.0 and tR <= Tobs);
            CHECK(giR == 0 or tR > prevRoman);
            prevRoman = tR;
            ++giR;
        }
        step.dt = step.tim - grid.back().tim;
        CHECK(step.dt > 0.0);
        grid.push_back(step);
    }
}
