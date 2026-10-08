// The run's output files and the provenance block.
#ifndef ROMAN_RUN_OUTPUTS_H
#define ROMAN_RUN_OUTPUTS_H

#include "common.h"
#include "run/config.h"
#include "run/sample_dump.h"
#include "run/sightlines.h"
#include "surveys/footprints.h"
#include "surveys/schedule.h"
#include "galaxy/extinction.h"

// Column names of the per-event table, in write order. Kept next to the writer's data so
// the two cannot drift: a header that disagrees with its columns is worse than none.
const char* eventTableHeader();

// The output streams and file names the sightline loop writes to, plus the sample-dump state.
// filg_in is deliberately left CLOSED: the loop opens it per event (see openOutputs).
struct RunOutputs {
    std::ofstream fil2, fil2b, fil3;   // EfLMC<tag>.dat, EfLMC<tag>B.dat, MapLMC<tag>.dat
    std::ofstream filg_in;             // per-event table, opened in append mode per event
    std::string   fnLDt;               // LpLMC<tag>.dat (appended by the loop)
    std::string   testf;               // the per-event table, test<tag>.dat
    std::string   fnPair;              // paired (satellite-parallax) side file
    SampleSpec             dumpSpec;   // the parsed --dump-samples spec
    std::vector<DumpEpoch> dumpBuf;    // the CURRENT draw's recorded epochs
    long                   dumpSeq = 0;// events written out so far
};

// Create/open the output files, write the table headers, set up the sample
// dump. Returns 0, or the exit code of the error printed.
int openOutputs(const RunConfig& cfg, RunOutputs& outs);

// Write run_provenance.txt (unless --dry-run) and print the same block to stdout.
// Returns 0, or the exit code of the error printed.
int writeRunProvenance(const RunConfig& cfg, const GbtdsLayout& gl, const GridSteps& steps,
                       const SightlineGrid& grid, const RomanSchedule& sched, const extin& ex);

#endif // ROMAN_RUN_OUTPUTS_H
