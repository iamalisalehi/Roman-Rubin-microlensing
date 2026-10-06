// Step S1: the illustrative-event dump (--dump-samples): which events to keep and how to write them.
#ifndef ROMAN_RUN_SAMPLE_DUMP_H
#define ROMAN_RUN_SAMPLE_DUMP_H

#include "common.h"
#include "types.h"

// ---------------------------------------------------------------------------
// Step S1. Sample-event light-curve and astrometry dump.
//
// The per-event table is one row per event. The per-epoch light curve and the
// astrometric track exist only inside the Monte Carlo's time loop and are
// overwritten by the next draw. That is enough for every pooled statistic this
// project reports, and not enough to DRAW one event -- which is what an
// illustrative figure in a paper needs.
//
// Two rules this dump obeys. Both are load-bearing.
//
//  1. IT CONSUMES NO RNG. Every value written is one the simulation had already
//     computed, including the noisy magnitude `magnio`, which the detection
//     chi-square has already drawn. The mt19937_64 stream is therefore identical
//     with and without --dump-samples, so turning the dump on cannot change WHICH
//     events get simulated. The legacy magC0/datC0 dump this supersedes called
//     RandN() inline and did not have that property. PROGRESS.md records what
//     stream divergence cost once already: the H3 run, started mid-stream, matched
//     0 of 2,673 v3 events.
//
//  2. IT BUFFERS, THEN COMMITS. Whether an event is "Rubin-only" or "both" is not
//     known until the time loop has ended and FisherM has run. So epochs are
//     accumulated in memory for EVERY draw and written only for the few that fill a
//     requested sample class. One event's buffer is a few thousand records.
// ---------------------------------------------------------------------------

// One recorded observation, in whichever observer frame the telescope that took it
// lives in: Rubin's values come from lightcurve(..., 0) (geocentric), Roman's from
// lightcurve(..., 1) (L2, ~0.01 AU further out). Keeping the two frames distinct is
// the point rather than an inconvenience -- the two observers see different impact
// parameters at the same instant, and that difference IS the satellite parallax.
struct DumpEpoch {
    double t;             //time [days]; day 0 = the first Rubin bulge visit
    int    tele;          //0 = Rubin, 1 = Roman. The same tag as lens::tele[]
    int    filt;          //0-5 = LSST ugrizy, 6 = Roman F146. ONE filter per visit
    double magObs;        //the simulated measurement: model + Gaussian noise [mag]
    double magMod;        //model magnitude, blended, WITH microlensing parallax [mag]
    double magMod0;       //the same model WITHOUT the parallax term [mag]
    double errMag;        //1-sigma photometric error from the instrument model [mag]
    double u;             //lens-source separation WITH parallax [Einstein radii]
    double u0;            //lens-source separation WITHOUT parallax [Einstein radii]
    double A;             //source magnification from u  [dimensionless, >= 1]
    double A0;            //source magnification from u0 [dimensionless, >= 1]
    double def1c, def2c;  //centroid deflection WITH parallax [mas]
    double def1a, def2a;  //centroid deflection WITHOUT parallax [mas]
    double pos1b, pos2b;  //unlensed source position: proper motion + its own parallax [mas]
    double pos1c, pos2c;  //observed lensed centroid, = pos*b + def*c [mas]
    double lens1, lens2;  //lens position: its proper motion + its own parallax [mas]
    double errAst;        //1-sigma astrometric error from the instrument model [mas]
};

// One requested sample class: a named predicate over a FINISHED event, plus how many
// of them to keep. The numeric cuts live here rather than in the predicate so that one
// spec file can ask for two different tE windows without a recompile.
struct SampleClass {
    std::string name;
    int    quota    = 0;
    int    kept     = 0;
    double teMin    = -1.0;   //days; negative = no cut
    double teMax    = -1.0;   //days; negative = no cut
    double shiftMin = -1.0;   //mas;  negative = no cut
};

struct SampleSpec {
    bool        on       = false;
    std::string dir      = "samples";
    // Dense model-curve sampling. stepTE is a FRACTION OF tE, not a number of days,
    // because the sampling density that draws one event well is set by that event's own
    // timescale: 0.2 d steps are wasteful for a 900 d black-hole event and too coarse for
    // a 5 d one. This keeps the peak window at a fixed ~2*spanTE/stepTE points whatever
    // the event, so the file size cannot blow up on the long-tE population.
    double      stepTE   = DEFAULT_S2_STEP_TE;  //peak-window step, in units of tE
    double      spanTE   = DEFAULT_S2_SPAN_TE;   //peak-window half-width, in units of tE
    double      dtCoarse = DEFAULT_S2_DT_COARSE;   //step over the rest of the mission [days]
    std::vector<SampleClass> classes;
};

// Everything a selector is allowed to look at. All of it is already in the per-event
// table, so a class can be checked against a FINISHED run's table before spending a
// run on it -- which is how you find out a cut is empty without waiting for the run.
struct SampleFacts {
    int    detL, detR, detJ;  //per-survey and joint detection booleans
    double tE;                //Einstein crossing time [days]
    int    t0zone;            //0 = t0 in a Roman season, 1 = mid-mission gap, 2 = off-mission
    int    nepLpk, nepRpk;    //epochs within +-2 tE of the OBSERVED peak t0obs, per survey (Dev. 76)
    int    ndwL, ndwR;        //epochs over the WHOLE mission, per survey. ndwR == 0 means
                              //Roman never observed this sightline at all: it lies outside
                              //the GBTDS footprint
    int    okBRoman;          //did Roman's ASTROMETRIC Fisher matrix invert?
    double maxShift;          //largest |centroid deflection| at any recorded epoch [mas]
};

// Parse a sample-spec file (--dump-samples). Returns false (after printing why) on a malformed file.
bool parseSampleSpec(const std::string& path, SampleSpec& spec);
// Does a finished event satisfy a sample class?
bool sampleMatch(const SampleClass& c, const SampleFacts& f);
// Write one buffered event's light curve and astrometric track to the dump files.
void writeSampleEvent(const SampleSpec& spec, const std::string& cls,
                      const std::string& id, const std::vector<DumpEpoch>& buf,
                      const std::string& params,
                      source& s, lens& l, astromet& as);

#endif // ROMAN_RUN_SAMPLE_DUMP_H
