// Test of LOCAL ADDITION 2 to the vendored VBMicrolensing (external/VBMicrolensing/): the observers' positions
// computed once per sightline epoch (SetObserverGrid, SetParallaxCache) instead of on every call of
// ComputeParallax. The addition must change no result of any function, bit for bit.
//
// Two library objects are set up the way src/events/vbm_model.cpp sets its one up (same tables, JD - 2450000
// times, t_in_HJD = false, Tol, a1 = 0); one never hears of the cache (upstream's code path), the other has
// it on and a grid. For three sightlines in a row, each with an epoch grid laid out from the real Rubin and
// Roman visit times, both observers (Earth, Roman at L2), 50 random events each at t0 and t0 +- 1e-4 d, the
// point-source and finite-source light curves with and without centroids are evaluated by both on the grid
// interleaved with off-grid times, and every output array (magnification, the four centroid arrays, y1, y2)
// is compared with memcmp, together with the library's extrapolation flag. Also covered: a grid left over
// from the previous sightline, the cache switched off, a table reload (which must empty the cache), a
// descending (binary-search) access order, and grid times outside the Sun and Roman tables.
//
// Needs Baseline/BulgeBaseline.dat and Baseline/RomanBaseline.dat (git-ignored inputs) and the vendored files;
// run from the repo root:   make vbmcachetest && ./vbmcachetest
// Exit status 0 iff every comparison was identical and the cache hit and missed where it should.
#include "parameters.h"
#include "VBMicrolensingLibrary.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

static long long nCompareArrays = 0, nCompareValues = 0, nFlagChecks = 0, nBad = 0;
static int nfail = 0;

static void expect(bool ok, const char* what) {
    std::printf("%-86s %s\n", what, ok ? "ok" : "FAIL");
    if (!ok) ++nfail;
}

// The two objects.
static std::unique_ptr<VBMicrolensing> makeInstance() {
    auto v = std::make_unique<VBMicrolensing>();
    v->LoadESPLTable(PATH_VBM_ESPL_TABLE);
    std::string sun = PATH_VBM_SUN_TABLE;
    v->LoadSunTable(sun.data());
    std::string target = PATH_VBM_TARGET_FILE, dir = PATH_ROMAN_EPHEM_DIR;
    v->SetObjectCoordinates(target.data(), dir.data());
    v->t_in_HJD = false;
    v->Tol = VBM_TOL;
    v->a1 = 0.0;
    return v;
}

static void reloadSun(VBMicrolensing& v) {
    std::string sun = PATH_VBM_SUN_TABLE;
    v.LoadSunTable(sun.data());
}
static void reloadSatellites(VBMicrolensing& v) {
    std::string target = PATH_VBM_TARGET_FILE, dir = PATH_ROMAN_EPHEM_DIR;
    v.SetObjectCoordinates(target.data(), dir.data());
}
static void setTarget(VBMicrolensing& v, const char* text) {
    std::string s = text;
    v.SetObjectCoordinates(s.data());
}

// Column `col` of a whitespace table with a '#' header, as doubles.
static std::vector<double> readColumn(const char* path, int col) {
    std::ifstream f(path);
    std::vector<double> out;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() or line[0] == '#') continue;
        std::istringstream is(line);
        double x = 0.0;
        for (int i = 0; i <= col; ++i) is >> x;
        out.push_back(x);
    }
    return out;
}

// First and last JD - 2450000 between the $$SOE and $$EOE markers of a Horizons-format table.
static void tableRange(const char* path, double& first, double& last) {
    std::ifstream f(path);
    std::string line;
    bool in = false, any = false;
    while (std::getline(f, line)) {
        if (line.compare(0, 5, "$$SOE") == 0) { in = true; continue; }
        if (line.compare(0, 5, "$$EOE") == 0) break;
        if (!in) continue;
        const double jd = std::stod(line) - 2450000.0;
        if (!any) first = jd;
        last = jd;
        any = true;
    }
}

// An epoch grid the way the simulator lays one out (src/sim/observe.cpp): strictly increasing simulation
// days, here in the library's clock (t + SIM_DAY0_JD2450000, the expression vbmEvaluate uses). Rubin's visit
// days and Roman's epochs, thinned differently for each sightline, plus the padding steps either side.
static std::vector<double> makeGrid(const std::vector<double>& rubin, const std::vector<double>& roman, int k) {
    std::vector<double> t;
    for (double x : {-100.0, -40.0, 20.0}) t.push_back(x);
    for (std::size_t i = std::size_t(k); i < rubin.size(); i += 29 + 7 * std::size_t(k))
        if (rubin[i] >= 0.0 and rubin[i] <= Tobs) t.push_back(rubin[i]);
    for (std::size_t i = std::size_t(5 * k); i < roman.size(); i += 1500 + 211 * std::size_t(k))
        if (roman[i] >= 0.0 and roman[i] <= Tobs) t.push_back(roman[i]);
    t.push_back(Tobs + 40.0);
    t.push_back(Tobs + 100.0);
    std::sort(t.begin(), t.end());
    t.erase(std::unique(t.begin(), t.end()), t.end());
    for (double& x : t) x = x + SIM_DAY0_JD2450000;
    return t;
}

// The grid interleaved with its midpoints (off-grid times between every pair), still increasing.
static std::vector<double> withMidpoints(const std::vector<double>& g) {
    std::vector<double> t;
    for (std::size_t i = 0; i < g.size(); ++i) {
        t.push_back(g[i]);
        if (i + 1 < g.size()) t.push_back(0.5 * (g[i] + g[i + 1]));
    }
    return t;
}

struct Event {
    double u0, lntE, t0, lnrho, piN, piE, muN, muE, piS, tetE;
};

static std::vector<Event> makeEvents(std::mt19937_64& g, int n, double tFirst, double tLast) {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::vector<Event> ev;
    for (int i = 0; i < n; ++i) {
        Event e;
        e.u0 = (U(g) < 0.5 ? -1.0 : 1.0) * std::exp(std::log(1e-3) + U(g) * (std::log(1.5) - std::log(1e-3)));
        e.lntE = std::log(5.0) + U(g) * (std::log(300.0) - std::log(5.0));
        e.t0 = tFirst + U(g) * (tLast - tFirst);
        e.lnrho = std::log(1e-4) + U(g) * (std::log(0.05) - std::log(1e-4));
        const double pie = std::exp(std::log(3e-3) + U(g) * (std::log(0.5) - std::log(3e-3)));
        const double ang = 2.0 * M_PI * U(g);
        e.piN = pie * std::cos(ang);
        e.piE = pie * std::sin(ang);
        e.muN = -10.0 + 20.0 * U(g);
        e.muE = -10.0 + 20.0 * U(g);
        e.piS = 0.05 + 0.15 * U(g);
        e.tetE = 0.1 + 0.9 * U(g);
        ev.push_back(e);
    }
    return ev;
}

// The output arrays of one call.
struct Out {
    std::vector<double> mag, y1, y2, c1s, c2s, c1l, c2l;
    int flag = 0;
    void size(std::size_t n) {
        for (auto* v : {&mag, &y1, &y2, &c1s, &c2s, &c1l, &c2l}) v->assign(n, 0.0);
    }
};

enum Kind { PSPL, ESPL, PSPL_ASTRO, ESPL_ASTRO, NKIND };
static const char* kindName[NKIND] = {"PSPL", "ESPL", "PSPL astro", "ESPL astro"};

static void evaluate(VBMicrolensing& V, Kind kind, const Event& e, double t0, int sat, std::vector<double> t, Out& o) {
    const int n = int(t.size());
    o.size(std::size_t(n));
    V.satellite = sat;
    V.parallaxextrapolation = 0;
    const bool astro = (kind == PSPL_ASTRO or kind == ESPL_ASTRO);
    V.astrometry = astro;
    switch (kind) {
    case PSPL: {
        double pr[9] = {e.u0, e.lntE, t0, e.piN, e.piE, e.muN, e.muE, e.piS, e.tetE};
        V.PSPLLightCurveParallax(pr, t.data(), o.mag.data(), o.y1.data(), o.y2.data(), n);
        break; }
    case ESPL: {
        double pr[10] = {e.u0, e.lntE, t0, e.lnrho, e.piN, e.piE, e.muN, e.muE, e.piS, e.tetE};
        V.ESPLLightCurveParallax(pr, t.data(), o.mag.data(), o.y1.data(), o.y2.data(), n);
        break; }
    case PSPL_ASTRO: {
        double pr[9] = {e.u0, e.lntE, t0, e.piN, e.piE, e.muN, e.muE, e.piS, e.tetE};
        V.PSPLAstroLightCurve(pr, t.data(), o.mag.data(), o.c1s.data(), o.c2s.data(), o.c1l.data(), o.c2l.data(),
                              o.y1.data(), o.y2.data(), n);
        break; }
    case ESPL_ASTRO: {
        double pr[10] = {e.u0, e.lntE, t0, e.lnrho, e.piN, e.piE, e.muN, e.muE, e.piS, e.tetE};
        V.ESPLAstroLightCurve(pr, t.data(), o.mag.data(), o.c1s.data(), o.c2s.data(), o.c1l.data(), o.c2l.data(),
                              o.y1.data(), o.y2.data(), n);
        break; }
    default: break;
    }
    V.astrometry = false;
    V.satellite = 0;
    o.flag = V.parallaxextrapolation;
}

// memcmp of every array (the centroid ones only for the astrometric kinds) and of the flag.
static bool identical(Kind kind, const Out& a, const Out& b) {
    bool same = (a.flag == b.flag);
    ++nFlagChecks;
    const bool astro = (kind == PSPL_ASTRO or kind == ESPL_ASTRO);
    auto cmp = [&](const std::vector<double>& x, const std::vector<double>& y) {
        ++nCompareArrays;
        nCompareValues += (long long)x.size();
        if (x.size() != y.size() or std::memcmp(x.data(), y.data(), x.size() * sizeof(double)) != 0) same = false;
    };
    cmp(a.mag, b.mag);
    cmp(a.y1, b.y1);
    cmp(a.y2, b.y2);
    if (astro) { cmp(a.c1s, b.c1s); cmp(a.c2s, b.c2s); cmp(a.c1l, b.c1l); cmp(a.c2l, b.c2l); }
    return same;
}

// One case: the same call on the plain and on the cached object, compared. Returns the cached object's hits
// and misses in the call through the references.
static void compareCall(VBMicrolensing& plain, VBMicrolensing& cached, Kind kind, const Event& e, double t0, int sat,
                        const std::vector<double>& t, long long& hits, long long& misses, const char* label) {
    Out a, b;
    evaluate(plain, kind, e, t0, sat, t, a);
    const long long h0 = cached.pc_hits, m0 = cached.pc_misses;
    evaluate(cached, kind, e, t0, sat, t, b);
    hits += cached.pc_hits - h0;
    misses += cached.pc_misses - m0;
    if (!identical(kind, a, b)) {
        ++nBad;
        if (nBad < 10) std::printf("  MISMATCH: %s, %s, satellite %d, t0 %.6f, flags %d vs %d\n", label, kindName[kind], sat, t0, a.flag, b.flag);
    }
}

// All four kinds, both observers, t0 and t0 +- 1e-4, for every event, over the array t.
static void sweepEvents(VBMicrolensing& plain, VBMicrolensing& cached, const std::vector<Event>& ev,
                        const std::vector<double>& t, long long& hits, long long& misses, const char* label) {
    for (const Event& e : ev)
        for (int sat = 0; sat <= 1; ++sat)
            for (double d : {0.0, 1e-4, -1e-4})
                for (int k = 0; k < NKIND; ++k)
                    compareCall(plain, cached, Kind(k), e, e.t0 + d, sat, t, hits, misses, label);
}

int main() {
    for (const char* p : {PATH_VBM_ESPL_TABLE, PATH_VBM_SUN_TABLE, PATH_VBM_TARGET_FILE, PATH_BULGE_BASELINE, PATH_ROMAN_BASELINE}) {
        std::ifstream f(p);
        if (!f) { std::printf("cannot read %s (run from the repository root)\n", p); return 1; }
    }
    const std::vector<double> rubin = readColumn(PATH_BULGE_BASELINE, 5), roman = readColumn(PATH_ROMAN_BASELINE, 5);

    auto plain = makeInstance();
    auto cached = makeInstance();
    cached->SetParallaxCache(true);
    expect(plain->nsat == 1 and cached->nsat == 1 and plain->suntable and cached->suntable, "both objects set up (Roman as satellite 1)");

    const char* targets[3] = {"17:53:50.78 -28:50:55.64", "17:41:20.00 -29:44:10.00", "18:06:30.50 -26:55:20.00"};
    std::mt19937_64 rng(20261113);
    long long hitsAll = 0, missesAll = 0;

    // (1) three sightlines in a row, each with its own grid.
    std::vector<double> grids[3];
    for (int k = 0; k < 3; ++k) {
        grids[k] = makeGrid(rubin, roman, k);
        const std::vector<double> t = withMidpoints(grids[k]);
        setTarget(*plain, targets[k]);
        setTarget(*cached, targets[k]);
        cached->SetObserverGrid(grids[k].data(), int(grids[k].size()));
        const std::vector<Event> ev = makeEvents(rng, 50, grids[k].front(), grids[k].back());
        const long long bad0 = nBad;
        long long hits = 0, misses = 0;
        sweepEvents(*plain, *cached, ev, t, hits, misses, "grid + midpoints");
        // A call is one array of grid + midpoints: every other entry is on the grid, except Earth-only...
        // both observers use the same grid, so exactly the grid entries hit.
        const long long calls = (long long)ev.size() * 2 * 3 * NKIND;
        char what[160];
        std::snprintf(what, sizeof what, "sightline %d: %zu grid times, %lld calls identical (hits %lld, misses %lld)", k, grids[k].size(), calls, hits, misses);
        expect(nBad == bad0, what);
        // Not the exact count: the astrometric functions also evaluate at t + dt for the lens velocity and
        // the t0 caches add calls; the hit count must be a large share of the evaluations and nothing else.
        expect(hits > 0 and hits > misses / 2, "  the cache is used on the grid times");
        hitsAll += hits;
        missesAll += misses;

        // Descending order: the binary-search path.
        std::vector<double> rev(t.rbegin(), t.rend());
        for (std::size_t i = 0; i < ev.size(); i += 5)
            for (int sat = 0; sat <= 1; ++sat)
                compareCall(*plain, *cached, PSPL, ev[i], ev[i].t0, sat, rev, hits, misses, "descending");
        expect(nBad == bad0, "  descending access order identical");

        // A grid from the previous sightline still set: sightline B evaluated with A's grid.
        if (k == 1) {
            cached->SetObserverGrid(grids[0].data(), int(grids[0].size()));
            long long h = 0, m = 0;
            sweepEvents(*plain, *cached, makeEvents(rng, 10, grids[1].front(), grids[1].back()), t, h, m, "stale grid");
            std::snprintf(what, sizeof what, "  stale grid of sightline 0 under sightline 1: identical (hits %lld, misses %lld)", h, m);
            expect(nBad == bad0, what);
            cached->SetObserverGrid(grids[1].data(), int(grids[1].size()));
        }
    }

    // (2) the cache switched off with a grid set: no hit at all, same outputs.
    {
        const std::vector<double> t = withMidpoints(grids[2]);
        const long long bad0 = nBad;
        cached->SetParallaxCache(false);
        long long h = 0, m = 0;
        sweepEvents(*plain, *cached, makeEvents(rng, 10, grids[2].front(), grids[2].back()), t, h, m, "cache off");
        char what[160];
        std::snprintf(what, sizeof what, "cache off, grid set: identical, hits %lld (must be 0), misses %lld", h, m);
        expect(nBad == bad0 and h == 0, what);
        cached->SetParallaxCache(true);
        h = m = 0;
        sweepEvents(*plain, *cached, makeEvents(rng, 5, grids[2].front(), grids[2].back()), t, h, m, "cache back on");
        std::snprintf(what, sizeof what, "cache on again: identical, hits %lld (> 0), misses %lld", h, m);
        expect(nBad == bad0 and h > 0, what);
    }

    // (3) a table reload empties the cache (both tables), outputs unchanged; a new grid fills it again.
    {
        const std::vector<double> t = withMidpoints(grids[2]);
        const long long bad0 = nBad;
        char what[160];
        long long h = 0, m = 0;
        reloadSun(*plain);
        reloadSun(*cached);
        sweepEvents(*plain, *cached, makeEvents(rng, 5, grids[2].front(), grids[2].back()), t, h, m, "after LoadSunTable");
        std::snprintf(what, sizeof what, "after LoadSunTable: identical, hits %lld (must be 0), misses %lld", h, m);
        expect(nBad == bad0 and h == 0, what);
        cached->SetObserverGrid(grids[2].data(), int(grids[2].size()));
        h = m = 0;
        sweepEvents(*plain, *cached, makeEvents(rng, 5, grids[2].front(), grids[2].back()), t, h, m, "new grid");
        std::snprintf(what, sizeof what, "new grid after the reload: identical, hits %lld (> 0), misses %lld", h, m);
        expect(nBad == bad0 and h > 0, what);
        reloadSatellites(*plain);
        reloadSatellites(*cached);
        setTarget(*plain, targets[2]);
        setTarget(*cached, targets[2]);
        h = m = 0;
        sweepEvents(*plain, *cached, makeEvents(rng, 5, grids[2].front(), grids[2].back()), t, h, m, "after SetObjectCoordinates(file, dir)");
        std::snprintf(what, sizeof what, "after SetObjectCoordinates(file, dir): identical, hits %lld (must be 0), misses %lld", h, m);
        expect(nBad == bad0 and h == 0, what);
        cached->SetObserverGrid(grids[2].data(), int(grids[2].size()));
        h = m = 0;
        sweepEvents(*plain, *cached, makeEvents(rng, 5, grids[2].front(), grids[2].back()), t, h, m, "new grid 2");
        std::snprintf(what, sizeof what, "new grid after the second reload: identical, hits %lld (> 0), misses %lld", h, m);
        expect(nBad == bad0 and h > 0, what);
    }

    // (4) grid times outside the tables: outputs and flags must match too. One case per table end, so that a
    // wrong flag of one end cannot hide behind another's. The cache is compared with itself switched off (same
    // object): past the end of the Sun table the library reads the row beyond it, which holds whatever the
    // allocator left, different in two objects.
    {
        double sunFirst = 0, sunLast = 0, romFirst = 0, romLast = 0;
        tableRange(PATH_VBM_SUN_TABLE, sunFirst, sunLast);
        tableRange((std::string(PATH_ROMAN_EPHEM_DIR) + "/satellite1.txt").c_str(), romFirst, romLast);
        struct Case { const char* name; std::vector<double> out; int flagEarth, flagRoman; };
        const Case cases[4] = {
            {"before the Sun table", {sunFirst - 30.0, sunFirst - 0.5}, 1, 2},
            {"after the Sun table", {sunLast + 0.5, sunLast + 30.0}, 1, 2},
            {"before Roman's table", {romFirst - 20.0, romFirst - 1e-3}, 0, 2},
            {"after Roman's table", {romLast + 1e-3, romLast + 20.0}, 0, 2}};
        std::mt19937_64 r2(7);
        const std::vector<Event> ev = makeEvents(r2, 4, grids[0].front(), grids[0].back());
        for (const Case& c : cases) {
            std::vector<double> g = c.out;
            for (double x : grids[0]) if (x > romFirst and x < romLast and g.size() < 40) g.push_back(x);
            g.push_back(romFirst);
            g.push_back(romLast);
            std::sort(g.begin(), g.end());
            g.erase(std::unique(g.begin(), g.end()), g.end());
            const std::vector<double>& t = g;   // grid times only: the flags can then come from nothing but the cache
            const long long bad0 = nBad;
            cached->SetObserverGrid(g.data(), int(g.size()));
            long long h = 0, m = 0;
            bool flagsAsExpected = true;
            for (const Event& e : ev)
                for (int sat = 0; sat <= 1; ++sat)
                    for (int k = 0; k < NKIND; ++k) {
                        Out a, b;
                        cached->SetParallaxCache(false);
                        evaluate(*cached, Kind(k), e, e.t0, sat, t, a);
                        cached->SetParallaxCache(true);
                        const long long h0 = cached->pc_hits, m0 = cached->pc_misses;
                        evaluate(*cached, Kind(k), e, e.t0, sat, t, b);
                        h += cached->pc_hits - h0;
                        m += cached->pc_misses - m0;
                        // Roman (satellite 1) is outside its own table wherever the Sun table is, and its flag, 2, overrides the Earth's
                        const int want = sat == 0 ? c.flagEarth : c.flagRoman;
                        if (a.flag != want or b.flag != want) flagsAsExpected = false;
                        if (!identical(Kind(k), a, b)) { ++nBad; std::printf("  MISMATCH %s: %s satellite %d flags %d vs %d\n", c.name, kindName[k], sat, a.flag, b.flag); }
                    }
            char what[200];
            std::snprintf(what, sizeof what, "grid times %s: identical, flags %d (Earth) and %d (Roman) as upstream (hits %lld, misses %lld)", c.name, c.flagEarth, c.flagRoman, h, m);
            expect(nBad == bad0 and h > 0 and flagsAsExpected, what);
        }
    }

    std::printf("\ncomparisons: %lld output arrays (%lld values) memcmp-identical, %lld extrapolation flags compared, %lld mismatching calls\n",
                nCompareArrays, nCompareValues, nFlagChecks, nBad);
    std::printf("cache lookups of the cached object in the sightline sweeps: %lld hits, %lld misses\n", hitsAll, missesAll);
    expect(nBad == 0, "every comparison identical");
    return nfail == 0 ? 0 : 1;
}
