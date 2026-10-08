// Roman observing-season geometry: season windows derived from the visit list, and where an event's t0 falls.
#ifndef ROMAN_SURVEYS_SCHEDULE_H
#define ROMAN_SURVEYS_SCHEDULE_H

#include "common.h"
#include "surveys/visits.h"

// Roman observing-season geometry.
//
// Roman looks at the bulge only when the Sun angle permits, so its visit list is a comb of
// ~70-day seasons separated by ~110-day gaps. Where an event's peak falls relative to them is the
// independent variable of the gap-filling result. The windows are derived from the epoch times in
// RomanBaseline.dat, never restated as constants, so they cannot disagree with the visit list.
// Clustering rule and SEASON_GAP_MIN_DAYS: see config/parameters.h, section 5.

// Where t0 sits relative to Roman's mission. Three states, not two: an event peaking before
// launch or after the end is Rubin-only by construction, whereas one peaking in a mid-mission gap
// is bracketed by Roman photometry and its wings are still measured. The two must not be pooled.
enum T0Zone {
    T0_IN_SEASON   = 0, //Roman was observing at t0
    T0_IN_GAP      = 1, //between two Roman seasons -- the gap-filling regime
    T0_OFF_MISSION = 2, //before Roman's first epoch or after its last
};

struct RomanSchedule {
    std::vector<std::pair<double,double>> seasons; //[start, end] in simulation days
    double missionStart       = 0.0;
    double missionEnd         = 0.0;
    double maxInSeasonSpacing = 0.0; //largest spacing kept INSIDE a season
    double minSeasonGap       = 0.0; //smallest spacing treated as a gap
    double minSeasonLength    = 0.0; //shortest season found, end - start [d]. Zero means some
                                     //"season" holds a single epoch (a schedule sampled more
                                     //coarsely than SEASON_GAP_MIN_DAYS).

    // Signed days from t0 to the nearest season boundary.
    //   negative -> t0 is INSIDE a season; |value| is how deep into it the peak sits
    //   positive -> t0 is outside every season; value is the distance to the nearest edge
    // Use zone() to tell a mid-mission gap from off-mission; both give a positive value.
    double dtToSeasonEdge(double t0) const;
    int    zone(double t0) const;
    int    seasonOf(double t) const;   // index into seasons containing t, or -1
};

// Cluster ro.tim into seasons (sorted, de-duplicated copy of the epoch times). Called once per run.
RomanSchedule buildRomanSchedule(const roman& ro);

#endif // ROMAN_SURVEYS_SCHEDULE_H
