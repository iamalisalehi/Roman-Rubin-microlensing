// The illustrative-event dump (--dump-samples): which events to keep and how to write them.
#ifndef ROMAN_RUN_SAMPLE_DUMP_H
#define ROMAN_RUN_SAMPLE_DUMP_H

#include "common.h"
#include "types.h"

// Sample-event light-curve and astrometry dump. The per-event table has one row per event; the
// per-epoch light curve and astrometric track exist only inside the Monte Carlo time loop. This
// dump keeps them for a few chosen events so they can be plotted. Two rules:
//
//  1. It consumes no RNG. Every value written was already computed, including the noisy
//     magnitude `magnio`, so the random stream is identical with and without --dump-samples.
//
//  2. It buffers, then commits. Whether an event is "Rubin-only" or "both" is not known until
//     the time loop has ended and FisherM has run, so epochs are accumulated for every draw
//     and written only for events that fill a requested sample class.

// One recorded observation, in the observer frame of the telescope that took it: Rubin's values
// come from evaluateTrack(..., 0, ...) (geocentric), Roman's from evaluateTrack(..., 1, ...) (L2). The two
// observers see different impact parameters at the same instant; that difference is the
// satellite parallax.
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

// One requested sample class: a named predicate over a finished event, plus how many of them
// to keep. The numeric cuts live here so one spec file can ask for different tE windows.
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
    // Dense model-curve sampling. stepTE is a fraction of tE, not a number of days, so the
    // peak window has a fixed ~2*spanTE/stepTE points whatever the event's timescale.
    double      stepTE   = DEFAULT_S2_STEP_TE;  //peak-window step, in units of tE
    double      spanTE   = DEFAULT_S2_SPAN_TE;   //peak-window half-width, in units of tE
    double      dtCoarse = DEFAULT_S2_DT_COARSE;   //step over the rest of the mission [days]
    std::vector<SampleClass> classes;
};

// Everything a selector may look at. All of it is in the per-event table, so a class can be
// checked against a finished run's table before spending a run on it.
struct SampleFacts {
    int    detL, detR, detJ;  //per-survey and joint detection booleans
    double tE;                //Einstein crossing time [days]
    int    t0zone;            //0 = t0 in a Roman season, 1 = mid-mission gap, 2 = off-mission
    int    nepLpk, nepRpk;    //epochs within +-2 tE of the OBSERVED peak t0obs, per survey
    int    ndwL, ndwR;        //epochs over the whole mission, per survey. ndwR == 0 means
                              //Roman never observed this sightline (outside the GBTDS footprint)
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
