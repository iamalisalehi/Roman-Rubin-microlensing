// Roman season clustering (buildRomanSchedule) and the RomanSchedule queries.
#include "surveys/schedule.h"

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///                 Roman observing-season geometry (Step D1)
///
/// Roman's visit list is a comb: ~70-day observing windows when the bulge is far enough
/// from the Sun, separated by ~110-day gaps when it is not. Everything downstream that
/// asks "did Roman have data near this event's peak?" needs those windows, and the only
/// authoritative statement of them is the epoch times themselves. So they are recovered
/// from the data rather than restated -- see the note above SEASON_GAP_MIN_DAYS in config/parameters.h.
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
RomanSchedule buildRomanSchedule(const roman& ro)
{
    // A sorted, de-duplicated copy: the file lists every FIELD at every epoch, so the
    // same instant appears once per field and the raw column is neither unique nor
    // guaranteed monotonic. ~2.4 MB and one sort, once per run.
    std::vector<double> t(ro.tim.begin(), ro.tim.end());
    std::sort(t.begin(), t.end());
    t.erase(std::unique(t.begin(), t.end()), t.end());

    RomanSchedule sch;
    if (t.empty()) return sch;

    sch.missionStart = t.front();
    sch.missionEnd   = t.back();
    sch.minSeasonGap = std::numeric_limits<double>::infinity();

    double start = t.front();
    for (size_t i = 1; i < t.size(); ++i) {
        const double d = t[i] - t[i - 1];
        if (d > SEASON_GAP_MIN_DAYS) {
            sch.seasons.emplace_back(start, t[i - 1]);
            start = t[i];
            // Track the two spacing populations the threshold is meant to separate, so
            // main() can verify the separation is real instead of trusting it.
            sch.minSeasonGap = std::min(sch.minSeasonGap, d);
        } else {
            sch.maxInSeasonSpacing = std::max(sch.maxInSeasonSpacing, d);
        }
    }
    sch.seasons.emplace_back(start, t.back());
    if (sch.seasons.size() == 1) sch.minSeasonGap = 0.0; //nothing was ever classed as a gap

    // Shortest season. A season holding a single epoch has zero length, which is the signature
    // of a schedule sampled more coarsely than the threshold: every epoch becomes its own
    // "season", maxInSeasonSpacing never gets updated at all (so it stays a healthy-looking 0)
    // and minSeasonGap stays above the threshold. Both other margins pass; only this catches it.
    sch.minSeasonLength = std::numeric_limits<double>::infinity();
    for (const auto& sea : sch.seasons)
        sch.minSeasonLength = std::min(sch.minSeasonLength, sea.second - sea.first);

    return sch;
}

double RomanSchedule::dtToSeasonEdge(double t0) const
{
    if (seasons.empty()) return 0.0;

    double best   = std::numeric_limits<double>::infinity();
    bool   inside = false;
    for (const auto& s : seasons) {
        if (t0 >= s.first and t0 <= s.second) inside = true;
        best = std::min(best, std::min(std::fabs(t0 - s.first), std::fabs(t0 - s.second)));
    }
    // Negative inside a season, positive outside. See the header comment on the sign.
    return inside ? -best : best;
}

int RomanSchedule::seasonOf(double t) const
{
    for (size_t k = 0; k < seasons.size(); ++k)
        if (t >= seasons[k].first and t <= seasons[k].second) return int(k);
    return -1;
}

int RomanSchedule::zone(double t0) const
{
    if (seasons.empty())                      return T0_OFF_MISSION;
    // Off-mission is tested FIRST and wins: a peak before launch or after the mission ends
    // is not a gap in any useful sense, however close it happens to sit to an edge.
    if (t0 < missionStart or t0 > missionEnd) return T0_OFF_MISSION;
    for (const auto& s : seasons)
        if (t0 >= s.first and t0 <= s.second) return T0_IN_SEASON;
    return T0_IN_GAP;
}
