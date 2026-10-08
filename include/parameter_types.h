#ifndef PARAMETER_TYPES_H
#define PARAMETER_TYPES_H

// The few TYPES that the hand-edited parameter table (config/parameters.h) is written in.
// Types only, no values: edit the values in config/parameters.h.

// ======================= WHICH LENS POPULATION IS BEING SIMULATED =======================
//
// A runtime object chosen by `--population`. It owns four things that must travel together:
// the mass function, the mass range it is defined on, the grid spacing of the mass-efficiency
// histogram, and the tag every output file carries.
//
// The tag names `test<tag>.dat`, `MapLMC<tag>.dat` and the rest, so two populations cannot
// overwrite each other and a table always says which population produced it.
enum class MassFunction {
    KROUPA_REMNANTS,   // Kroupa IMF, then the initial-final mass relation: the bulge today
    LOG_UNIFORM,       // flat in log M -- the honest prior when the mass function is unknown
    NEUTRON_STAR,      // a measured NS mass distribution (Ozel & Freire 2016)
    UNIFORM,           // legacy MACHO-search options, inherited from the LMC simulation
    POWER_LAW_05,      //   dN/dM ~ M^-0.5
    POWER_LAW_10,      //   dN/dM ~ M^-1
    POWER_LAW_20,      //   dN/dM ~ M^-2
    BESANCON_CATALOGUE // the lens is a random member of its component's Besancon list: mass and light
};

struct LensPopulation {
    const char*  name;      // what --population takes
    const char*  tag;       // output-file suffix
    MassFunction mf;
    double       mlMin;     // Msun; also the low edge of the Mls efficiency grid
    double       mlMax;
    bool         logGrid;   // log-spaced mass bins: mandatory once the range spans decades
    int          legacyId;  // legacy mass-function id, used by the two debug dumps
    const char*  note;
};

// Astrometric noise-model variants: W white, N nominal, P pessimistic. Their
// definitions and the AST_SIGC_* values are documented in config/parameters.h, section 5.
enum AstroVariant { AV_W = 0, AV_N = 1, AV_P = 2, NAVAR = 3 };

#endif // PARAMETER_TYPES_H
