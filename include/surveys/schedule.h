// Roman observing-season geometry: season windows derived from the visit list, and where an event's t0 falls.
#ifndef ROMAN_SURVEYS_SCHEDULE_H
#define ROMAN_SURVEYS_SCHEDULE_H

#include "common.h"
#include "surveys/visits.h"

// ---------------------------------------------------------------------------
// Roman observing-season geometry (Step D1).
//
// Roman can only look at the bulge when the Sun angle permits, so its ten-year visit
// list is a comb of ~70-day observing seasons separated by ~110-day gaps. Where an
// event's peak falls relative to those seasons is the independent variable of the
// gap-filling result, so it has to be computed at simulation time and stored.
//
// The windows are DERIVED from the epoch times actually present in RomanBaseline.dat,
// never restated as constants here. The schedule already lives in
// Baseline/generateRomanBaseline.py; a second copy would drift silently, and it is
// expected to change (OPEN_ITEMS.md: the GBTDS footprint and cadence are still being
// reconciled against STScI's current pages). Deriving means the C++ can never disagree
// with the visit list it is actually integrating.
//
// Clustering rule and the SEASON_GAP_MIN_DAYS constant: see config/parameters.h, section 5.

// Where t0 sits relative to Roman's mission. Three states, not two.
//
// An event peaking before Roman launches, or after it ends, is Rubin-only BY
// CONSTRUCTION: there is no Roman data anywhere near it and nothing for a joint fit to
// rescue. An event peaking in a MID-MISSION gap is a completely different object --
// Roman brackets it, with dense photometry on both sides, so a long-tE event's wings
// are still measured even though its peak was missed. That second case is the one the
// joint-fit science claim is about. Collapsing the two into a single "Roman had no data
// at t0" boolean would dilute the headline result with events that were never
// candidates, which is the easiest available way to wash the effect out.
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
                                     //"season" holds a single epoch, which is not a season at
                                     //all -- see the guard in main(). A schedule sampled more
                                     //coarsely than SEASON_GAP_MIN_DAYS degenerates that way
                                     //while leaving the other two margins looking healthy.

    // Signed days from t0 to the nearest season boundary.
    //   negative -> t0 is INSIDE a season; |value| is how deep into it the peak sits
    //   positive -> t0 is outside every season; value is the distance to the nearest edge
    // Signed this way so the gap-filling plot reads left to right: x < 0 is "Roman was
    // watching", x > 0 is "Roman was not", and the joint-over-Roman precision gain is
    // expected to grow with x. Use zone() to tell a mid-mission gap from off-mission;
    // both give a positive dtToSeasonEdge and they must not be pooled.
    double dtToSeasonEdge(double t0) const;
    int    zone(double t0) const;
    int    seasonOf(double t) const;   // index into seasons containing t, or -1 (Deviation 71)
};

// Cluster ro.tim into seasons. One pass over a sorted, de-duplicated copy of the epoch
// times; called once per run, not per event.
RomanSchedule buildRomanSchedule(const roman& ro);

#endif // ROMAN_SURVEYS_SCHEDULE_H
