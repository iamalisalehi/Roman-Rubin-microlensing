// The light-curve model's interface to the vendored VBMicrolensing library; see events/vbm_model.h.
// The only file that includes the library.
#include "events/vbm_model.h"
#include "VBMicrolensingLibrary.h"

namespace {

// One library object for the whole run. It holds ~0.5 MB of tables, so it lives on the heap; it is built
// only when the first sightline needs it.
std::unique_ptr<VBMicrolensing> gVbm;
bool gSightlineSet = false;
// First and last row of the Sun table, as JD - 2450000 (the library's clock). The library reads the Earth's
// position and velocity from rows around the requested time without checking the range.
double gSunFirst = 0.0, gSunLast = 0.0;

VBMicrolensing& vbm() {
    CHECK(gVbm != nullptr);
    return *gVbm;
}

[[noreturn]] void fail(const std::string& what) { throw std::runtime_error("VBMicrolensing setup: " + what); }

// A file the library would only warn about (and then compute with an empty table).
void requireReadable(const std::string& path, const char* what) {
    std::ifstream f(path);
    if (!f) fail("cannot read " + path + " (" + what + "; run from the repository root)");
}

// First and last JD - 2450000 between the $$SOE and $$EOE markers of a Horizons-format table, the layout
// the library reads: one row per line, the Julian date first.
void tableRange(const std::string& path, double& first, double& last) {
    std::ifstream f(path);
    std::string line;
    bool inside = false;
    first = last = 0.0;
    int rows = 0;
    while (std::getline(f, line)) {
        if (line.compare(0, 5, "$$SOE") == 0) { inside = true; continue; }
        if (line.compare(0, 5, "$$EOE") == 0) break;
        if (!inside) continue;
        double jd;
        std::istringstream row(line);
        if (!(row >> jd)) fail(path + ": a row between $$SOE and $$EOE does not start with a Julian date");
        if (rows++ == 0) first = jd - 2450000.0;
        last = jd - 2450000.0;
    }
    if (rows == 0) fail(path + ": no rows between $$SOE and $$EOE");
}

// ICRS -> Galactic rotation of the Hipparcos catalogue (ESA 1997, The Hipparcos and Tycho Catalogues, vol. 1,
// sec. 1.5.3, eq. 1.5.11): r_gal = A r_icrs, with the north Galactic pole at (RA, Dec) = (192.85948, +27.12825)
// deg and the Galactic longitude of the north celestial pole 122.93192 deg, defined in the ICRS itself. The
// rows are the Galactic x (toward the centre), y (toward l = 90) and z (north) axes in the ICRS.
constexpr double A_ICRS_TO_GAL[3][3] = {
    {-0.0548755604162154, -0.8734370902348850, -0.4838350155487132},
    {+0.4941094278755837, -0.4448296299600112, +0.7469822444972189},
    {-0.8676661490190047, -0.1980763734312015, +0.4559837761750669}};

// A vector's Galactic Cartesian components -> ICRS: r_icrs = A^T g.
void galToIcrs(const double g[3], double r[3]) {
    for (int i = 0; i < 3; ++i)
        r[i] = A_ICRS_TO_GAL[0][i] * g[0] + A_ICRS_TO_GAL[1][i] * g[1] + A_ICRS_TO_GAL[2][i] * g[2];
}

double dot(const double a[3], const double b[3]) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

}  // namespace

void vbmInit() {
    if (gVbm) return;

    const std::string esplPath = PATH_VBM_ESPL_TABLE, sunPath = PATH_VBM_SUN_TABLE;
    const std::string targetPath = PATH_VBM_TARGET_FILE;
    const std::string romanPath = std::string(PATH_ROMAN_EPHEM_DIR) + "/satellite1.txt";
    requireReadable(esplPath, "the finite-source magnification table");
    requireReadable(sunPath, "the Sun ephemeris");
    requireReadable(targetPath, "the placeholder target; python tools/make_roman_l2_ephemeris.py writes it");
    requireReadable(romanPath, "Roman's ephemeris; python tools/make_roman_l2_ephemeris.py writes it");

    // The ephemerides must cover the simulation, [0, Tobs] days. The library extrapolates silently past a
    // satellite table, and reads outside the Sun table, which holds the rows around the time it is asked
    // for (up to two beyond it).
    double sunFirst, sunLast, romanFirst, romanLast;
    tableRange(sunPath, sunFirst, sunLast);
    tableRange(romanPath, romanFirst, romanLast);
    if (sunFirst + 1.0 > SIM_DAY0_JD2450000 || sunLast - 2.0 < SIM_DAY0_JD2450000 + Tobs)
        fail(sunPath + " does not cover the simulation");
    if (romanFirst > SIM_DAY0_JD2450000 || romanLast < SIM_DAY0_JD2450000 + Tobs)
        fail(romanPath + " does not cover the simulation; regenerate it with tools/make_roman_l2_ephemeris.py");

    auto v = std::make_unique<VBMicrolensing>();
    v->LoadESPLTable(PATH_VBM_ESPL_TABLE);
    std::string sun = PATH_VBM_SUN_TABLE;    // the library's setters take non-const char*
    v->LoadSunTable(sun.data());
    // Reads the target (the placeholder: the real one is set per sightline) and every satellite<k>.txt in
    // the directory; satellite = k then selects table k, 0 being the Earth.
    std::string target = PATH_VBM_TARGET_FILE, dir = PATH_ROMAN_EPHEM_DIR;
    v->SetObjectCoordinates(target.data(), dir.data());
    if (v->ESPLoff || !v->suntable || !v->AreCoordinatesSet() || v->nsat != 1)
        fail("the library did not load its tables (ESPL " + std::string(v->ESPLoff ? "no" : "yes") +
             ", Sun table " + (v->suntable ? "yes" : "no") + ", target " + (v->AreCoordinatesSet() ? "yes" : "no") +
             ", satellites " + std::to_string(v->nsat) + " of 1)");

    // Times are JD - 2450000 as given, not heliocentric (the library would subtract the light-travel time
    // to the Sun and then evaluate the ephemeris there).
    v->t_in_HJD = false;
    v->Tol = VBM_TOL;

    gSunFirst = sunFirst;
    gSunLast = sunLast;
    gSightlineSet = false;
    gVbm = std::move(v);
}

SkyFrame vbmSetSightline(double lDeg, double bDeg) {
    vbmInit();

    // The Galactic unit vectors toward the sightline, increasing l and increasing b, rotated to the ICRS.
    const double l = lDeg / RAa, b = bDeg / RAa;
    const double gn[3] = {std::cos(b) * std::cos(l), std::cos(b) * std::sin(l), std::sin(b)};
    const double gl[3] = {-std::sin(l), std::cos(l), 0.0};
    const double gb[3] = {-std::sin(b) * std::cos(l), -std::sin(b) * std::sin(l), std::cos(b)};
    double n[3], el[3], eb[3];
    galToIcrs(gn, n);
    galToIcrs(gl, el);
    galToIcrs(gb, eb);

    SkyFrame f;
    f.lDeg = lDeg;
    f.bDeg = bDeg;
    double ra = std::atan2(n[1], n[0]);
    if (ra < 0.0) ra += 2.0 * pi;
    const double dec = std::atan2(n[2], std::hypot(n[0], n[1]));
    f.raDeg = ra * RAa;
    f.decDeg = dec * RAa;

    // Celestial North and East at the sightline, and the Galactic axes on them.
    const double eN[3] = {-std::sin(dec) * std::cos(ra), -std::sin(dec) * std::sin(ra), std::cos(dec)};
    const double eE[3] = {-std::sin(ra), std::cos(ra), 0.0};
    f.toNE[0][0] = dot(eN, el);  f.toNE[0][1] = dot(eN, eb);
    f.toNE[1][0] = dot(eE, el);  f.toNE[1][1] = dot(eE, eb);

    // The library's coordinate string, "hh:mm:ss.ssssss +dd:mm:ss.sssss", built from integers (millionths
    // of a second of time, hundred-thousandths of an arcsecond) so that rounding carries into the minutes.
    // Its parser takes the sign from the degrees field, which for -00 is lost, so a declination within a
    // degree of the equator cannot be written.
    if (std::fabs(f.decDeg) < 1.0) fail("sightline too close to the celestial equator for the library's coordinate string");
    const long long perDay = 24LL * 3600 * 1000000;
    const long long ras = std::llround(f.raDeg / 15.0 * 3600.0 * 1e6) % perDay;
    const long long decs = std::llround(std::fabs(f.decDeg) * 3600.0 * 1e5);
    char text[96];
    std::snprintf(text, sizeof text, "%02lld:%02lld:%02lld.%06lld %c%02lld:%02lld:%02lld.%05lld",
                  ras / 3600000000LL, ras / 60000000LL % 60, ras / 1000000LL % 60, ras % 1000000LL,
                  f.decDeg < 0.0 ? '-' : '+',
                  decs / 360000000LL, decs / 6000000LL % 60, decs / 100000LL % 60, decs % 100000LL);
    // Every light-curve call of the library drops its parallax caches (Earth's position and velocity at t0,
    // which belong to the old target) before it starts, so changing the target between calls is safe.
    vbm().SetObjectCoordinates(text);
    if (!vbm().AreCoordinatesSet()) fail(std::string("the library rejected the coordinates \"") + text + "\"");
    gSightlineSet = true;
    return f;
}

void earthVelocityNE(double tSimDay, double v[2]) {
    CHECK(gSightlineSet);
    const double t = tSimDay + SIM_DAY0_JD2450000;
    CHECK(t >= gSunFirst + 1.0 && t <= gSunLast - 2.0);
    double south_west[2];
    vbm().EarthVelocityOnSky(t, south_west);
    v[0] = -south_west[0];
    v[1] = -south_west[1];
}

int vbmTestSatelliteShift(double tSimDay, double piEN, double piEE, double dy[2]) {
    CHECK(gSightlineSet);
    VBMicrolensing& V = vbm();
    double t = tSimDay + SIM_DAY0_JD2450000;     // the library's array arguments are non-const
    // u0, ln tE, t0, piE_N, piE_E. Any trajectory will do: both observers share it, only the difference
    // of their tracks is returned.
    double pr[5] = {0.3, std::log(50.0), t - 30.0, piEN, piEE};
    double mag, y1[2], y2[2];
    int flag = 0;
    for (int sat = 0; sat < 2; ++sat) {
        V.satellite = sat;
        V.PSPLLightCurveParallax(pr, &t, &mag, &y1[sat], &y2[sat], 1);
        flag = std::max(flag, V.parallaxextrapolation);
    }
    V.satellite = 0;
    dy[0] = y1[1] - y1[0];
    dy[1] = y2[1] - y2[0];
    return flag;
}
