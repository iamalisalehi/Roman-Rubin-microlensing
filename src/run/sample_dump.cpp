// Sample-event dump: spec parsing, class matching, and writing the light curves.
#include "run/sample_dump.h"
#include "events/lightcurve.h"

static bool knownSampleClass(const std::string& n)
{
    return n == "any" or n == "ns_typical" or n == "bh_short" or n == "bh_long"
        or n == "both" or n == "rubin_only" or n == "roman_only"
        or n == "gap_filler" or n == "astrometric";
}

// Does this finished event belong to class `c`?
//
// `any`, `ns_typical`, `bh_short` and `bh_long` are the same predicate (the joint fit detected
// it), separated only by the tE cuts in the spec file; the lens population comes from
// --population, and the names only label the output files.
//
// The "only" classes require epochs from the other survey within +-2 tE of the peak
// (nepR_pk or nepL_pk > 0); otherwise "Rubin-only" would just mean the event lay outside
// Roman's footprint or in a season gap.
bool sampleMatch(const SampleClass& c, const SampleFacts& f)
{
    if (c.teMin    > 0.0 and f.tE       < c.teMin)    return false;
    if (c.teMax    > 0.0 and f.tE       > c.teMax)    return false;
    if (c.shiftMin > 0.0 and f.maxShift < c.shiftMin) return false;

    // These classes show what both telescopes see of one kind of event, so both must have
    // observed the field.
    if (c.name == "any" or c.name == "ns_typical"
        or c.name == "bh_short" or c.name == "bh_long")
        return f.detJ == 1 and f.ndwL > 0 and f.ndwR > 0;
    if (c.name == "both")        return f.detL == 1 and f.detR == 1;
    if (c.name == "rubin_only")  return f.detL == 1 and f.detR == 0 and f.nepRpk > 0;
    if (c.name == "roman_only")  return f.detR == 1 and f.detL == 0 and f.nepLpk > 0;
    // Gap filling: the peak fell in a mid-mission Roman gap (t0zone 1, not 2) and Rubin alone
    // caught it. ndwR > 0 is required because t0zone comes from the mission-wide season
    // schedule, which ignores pointing: it reads 1 for a sightline Roman never visits too.
    if (c.name == "gap_filler")
        return f.detL == 1 and f.detR == 0 and f.t0zone == 1 and f.ndwR > 0;
    // Roman's astrometric matrix inverted, so the centroid ellipse is a measurement.
    if (c.name == "astrometric") return f.detJ == 1 and f.okBRoman == 1;
    return false;
}

// Reads the spec file. One directive per line; '#' starts a comment; blank lines ignored.
//
//     dir        samples/bh     # where the per-event files go
//     step_te    0.01           # peak-window step, as a FRACTION of this event's tE
//     span_te    3.0            # peak-window half-width, in units of tE
//     dt_coarse  2.0            # step over the rest of the mission [days]
//     class  <name>  <quota>  [te_min=X] [te_max=X] [shift_min=X]
//
// An unknown class name or key is a hard error, so a typo does not silently produce no samples.
bool parseSampleSpec(const std::string& path, SampleSpec& spec)
{
    std::ifstream in(path);
    if (!in) { std::cerr << "ERROR: cannot open sample spec '" << path << "'\n"; return false; }

    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        const std::size_t hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        std::istringstream ls(line);
        std::string key;
        if (!(ls >> key)) continue;

        if (key == "dir")            { ls >> spec.dir;      continue; }
        if (key == "step_te")        { ls >> spec.stepTE;   continue; }
        if (key == "span_te")        { ls >> spec.spanTE;   continue; }
        if (key == "dt_coarse")      { ls >> spec.dtCoarse; continue; }
        if (key != "class") {
            std::cerr << "ERROR: " << path << ":" << lineNo << ": unknown directive '"
                      << key << "'\n";
            return false;
        }
        SampleClass c;
        if (!(ls >> c.name >> c.quota)) {
            std::cerr << "ERROR: " << path << ":" << lineNo << ": 'class' needs a name "
                      << "and a quota\n";
            return false;
        }
        if (!knownSampleClass(c.name)) {
            std::cerr << "ERROR: " << path << ":" << lineNo << ": unknown class '"
                      << c.name << "'. Known: any ns_typical bh_short bh_long both "
                      << "rubin_only roman_only gap_filler astrometric\n";
            return false;
        }
        std::string kv;
        while (ls >> kv) {
            const std::size_t eq = kv.find('=');
            if (eq == std::string::npos) {
                std::cerr << "ERROR: " << path << ":" << lineNo << ": expected key=value, got '"
                          << kv << "'\n";
                return false;
            }
            const std::string k = kv.substr(0, eq);
            const double      v = std::atof(kv.substr(eq + 1).c_str());
            if      (k == "te_min")    c.teMin    = v;
            else if (k == "te_max")    c.teMax    = v;
            else if (k == "shift_min") c.shiftMin = v;
            else {
                std::cerr << "ERROR: " << path << ":" << lineNo << ": unknown key '" << k
                          << "'. Known: te_min te_max shift_min\n";
                return false;
            }
        }
        spec.classes.push_back(c);
    }
    if (spec.classes.empty()) {
        std::cerr << "ERROR: " << path << " requested no classes\n";
        return false;
    }
    if (!(spec.stepTE > 0.0) or !(spec.spanTE > 0.0) or !(spec.dtCoarse > 0.0)) {
        std::cerr << "ERROR: " << path << ": step_te, span_te and dt_coarse must be > 0\n";
        return false;
    }
    spec.on = true;
    return true;
}

// Writes the three files that describe one sample event:
//
//   <class>_<id>_epochs.dat  what the two surveys actually recorded, one row per visit
//   <class>_<id>_model.dat   a dense noise-free model curve, in BOTH observer frames
//   <class>_<id>_params.dat  the event's true parameters and its forecast sigmas
//
// The dense curve is regenerated here rather than taken from the time loop, whose adaptive
// step goes coarse in the wings. It costs one lightcurve() call per grid point and consumes
// no RNG.
void writeSampleEvent(const SampleSpec& spec, const std::string& cls,
                             const std::string& id, const std::vector<DumpEpoch>& buf,
                             const std::string& params,
                             source& s, lens& l, astromet& as)
{
    const std::string stem = spec.dir + "/" + cls + "_" + id;

    {
        std::ofstream fp(stem + "_params.dat");
        fp << params;
    }

    {
        std::ofstream fe(stem + "_epochs.dat");
        fe << "# Recorded observations. tele: 0 = Rubin, 1 = Roman. filt: 0-5 = ugrizy, "
              "6 = F146.\n"
           << "# mag_obs is the simulated measurement the detection test used; mag_mod and\n"
           << "# mag_mod_noplx are the blended model with and without microlensing parallax.\n"
           << "# Positions and deflections are in mas, in the frame of the telescope that\n"
           << "# took the epoch (Roman's is the L2 frame).\n"
           << "# t tele filt mag_obs mag_mod mag_mod_noplx err_mag u u_noplx A A_noplx "
              "def1c def2c def1a def2a pos1b pos2b pos1c pos2c lens1 lens2 err_ast\n";
        fe << std::setprecision(8);
        for (const DumpEpoch& e : buf)
            fe << e.t     << " " << e.tele  << " " << e.filt  << " "
               << e.magObs<< " " << e.magMod<< " " << e.magMod0 << " " << e.errMag << " "
               << e.u     << " " << e.u0    << " " << e.A      << " " << e.A0     << " "
               << e.def1c << " " << e.def2c << " " << e.def1a  << " " << e.def2a  << " "
               << e.pos1b << " " << e.pos2b << " " << e.pos1c  << " " << e.pos2c  << " "
               << e.lens1 << " " << e.lens2 << " " << e.errAst << "\n";
    }

    // Two-tier grid: a fine window around the peak, and a coarse grid over the whole mission
    // so the year-scale parallax wobble in the wings is kept.
    std::vector<double> grid;
    const double tStart = 0.0  * year - 100.0;
    const double tEnd   = 10.0 * year + 100.0;
    for (double t = tStart; t <= tEnd; t += spec.dtCoarse) grid.push_back(t);
    const double half = spec.spanTE * l.tE;
    const double fine = spec.stepTE * l.tE;
    for (double t = l.t0 - half; t <= l.t0 + half; t += fine)
        if (t >= tStart and t <= tEnd) grid.push_back(t);
    std::sort(grid.begin(), grid.end());

    std::ofstream fm(stem + "_model.dat");
    fm << "# Dense noise-free model curve. frame: 0 = geocentric (Rubin), 1 = L2 (Roman).\n"
       << "# Both frames are written for the WHOLE grid, so the satellite-parallax\n"
       << "# difference can be read off directly as frame 1 minus frame 0 at equal t.\n"
       << "# mag_* are blended model magnitudes per FILTER (u g r i z y F146); mag0_* are\n"
       << "# the same without the microlensing parallax term. Positions in mas.\n"
       << "# t frame u u_noplx A A_noplx mag_u..mag_F146 mag0_u..mag0_F146 "
          "def1c def2c def1a def2a pos1b pos2b pos1c pos2c lens1 lens2\n";
    fm << std::setprecision(8);
    for (int frame = 0; frame < 2; ++frame) {
        for (double t : grid) {
            lightcurve(s, l, as, t, frame);
            const double A  = magnifOf(s.ut);
            const double A0 = magnifOf(s.ut0);
            fm << t << " " << frame << " "
               << s.ut << " " << s.ut0 << " " << A << " " << A0;
            // All M filters are safe here: they are built from A directly, not read from the
            // time loop's magni[]/magni0[] scratch arrays, which at a Roman epoch still hold
            // Rubin-frame values in slots 0-5.
            for (int i = 0; i < M; ++i)
                fm << " " << s.magb[i] - 2.5 * std::log10(A  * s.blend[i] + 1.0 - s.blend[i]);
            for (int i = 0; i < M; ++i)
                fm << " " << s.magb[i] - 2.5 * std::log10(A0 * s.blend[i] + 1.0 - s.blend[i]);
            fm << " " << s.def1c << " " << s.def2c << " " << s.def1a << " " << s.def2a
               << " " << s.pos1b << " " << s.pos2b << " " << s.pos1c << " " << s.pos2c
               << " " << l.pos1  << " " << l.pos2  << "\n";
        }
    }
}
