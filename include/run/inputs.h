// Reading the input files: the Rubin and Roman visit lists, the noise tables, the extinction
// tables, the CMD, the luminous-lens table and the camera footprint.
#ifndef ROMAN_RUN_INPUTS_H
#define ROMAN_RUN_INPUTS_H

#include "common.h"
#include "types.h"
#include "surveys/visits.h"
#include "surveys/footprints.h"
#include "surveys/schedule.h"
#include "galaxy/catalogue.h"
#include "galaxy/extinction.h"

// Each reader prints its progress line and returns 0, or the process exit code of the error it
// reported (the run is over when it is non-zero).

// The GBTDS detector layout, read first: the footprint grid is checked against the
// real detector size, and the scan region is built from the field outline. Exits on bad input.
GbtdsLayout loadGbtdsLayout();

int readRubinVisits(lsst& ls);            // BulgeBaseline.dat
int readRubinAstromTable(lsst& ls);       // sigmaA_LSST.txt
int readRomanErrorTable(roman& ro);       // sigma_Roman.txt, rescaled to ROMAN_DEPTH5_AB
int readRomanVisits(roman& ro);           // RomanBaseline.dat
int buildRomanSeasons(const roman& ro, RomanSchedule& sched);   // season geometry + refuse-to-run guard
int readSkyTables(extin& ex, const lsst& ls);   // extinction, luminous lenses, Rubin camera map, galToIcrs check
int readCmdTables(CMD& cm);               // read_cmd

// Everything above, in this order (the order is the order of the printed log).
int loadInputs(lsst& ls, roman& ro, extin& ex, CMD& cm, RomanSchedule& sched);

#endif // ROMAN_RUN_INPUTS_H
