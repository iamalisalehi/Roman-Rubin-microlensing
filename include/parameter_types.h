#ifndef PARAMETER_TYPES_H
#define PARAMETER_TYPES_H

// The few TYPES that the hand-edited parameter table (config/parameters.h) is written in.
// Types only, no values: edit the values in config/parameters.h.

// ======================= WHICH LENS POPULATION IS BEING SIMULATED =======================
//
// This was a compile-time constant (`IMnum`), which meant a second population needed a
// rebuild and silently reused the first one's output filenames if you forgot to change it.
// It is now a runtime object chosen by `--population`, because the whole point of a
// population study is to run several and compare them.
//
// A population owns FOUR things, and they must travel together or a run is mislabelled:
//   the mass function it draws from, the mass range that function is defined on, the grid
//   spacing the mass-efficiency histogram uses, and the tag every output file carries.
//
// THE TAG IS THE SAFETY CATCH. `test<tag>.dat`, `MapLMC<tag>.dat` and the rest are named
// from it, so two populations cannot overwrite each other, and a table always says which
// population produced it. The default keeps tag "5", so existing filenames and every
// analysis command that names them are unchanged.
enum class MassFunction {
    KROUPA_REMNANTS,   // Kroupa IMF, then the initial-final mass relation: the bulge today
    LOG_UNIFORM,       // flat in log M -- the honest prior when the mass function is unknown
    NEUTRON_STAR,      // a measured NS mass distribution (Ozel & Freire 2016)
    UNIFORM,           // legacy MACHO-search options, inherited from the LMC simulation
    POWER_LAW_05,      //   dN/dM ~ M^-0.5
    POWER_LAW_10,      //   dN/dM ~ M^-1
    POWER_LAW_20       //   dN/dM ~ M^-2
};

struct LensPopulation {
    const char*  name;      // what --population takes
    const char*  tag;       // output-file suffix
    MassFunction mf;
    double       mlMin;     // Msun; also the low edge of the Mls efficiency grid
    double       mlMax;
    bool         logGrid;   // log-spaced mass bins: mandatory once the range spans decades
    int          legacyId;  // the old IMnum, for the two IMnum==1 debug dumps
    const char*  note;
};

// Astrometric noise-model variants (Deviation 71): W white, N nominal, P pessimistic. Their
// definitions and the AST_SIGC_* values are documented in config/parameters.h, section 5.
enum AstroVariant { AV_W = 0, AV_N = 1, AV_P = 2, NAVAR = 3 };

#endif // PARAMETER_TYPES_H
