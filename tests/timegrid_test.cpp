// Unit test for the merged time grid (src/sim/timegrid.cpp buildTimeGrid). Needs no data files.
// Exit 0 = all held.
//
// Every case builds the matched-visit lists the way matchVisibleEpochs leaves them: ct[k] indexes a larger
// tim array whose other entries are unmatched visits, so a builder that read tim[k] instead of tim[ct[k]]
// would fail. Times are compared with == (the grid copies the lists' times; no arithmetic is done on them);
// the step lengths are differences and are compared to 1e-12.
#include "common.h"
#include "sim/timegrid.h"
#include <cstdio>
#include <cmath>
#include <random>
#include <string>

static int nfail = 0;
static void pass(const char* name) { std::printf("PASS %s\n", name); }
static void fail(const char* name, const std::string& what) {
    std::printf("FAIL %s: %s\n", name, what.c_str());
    ++nfail;
}

// A survey's matched list: the full tim array (matched times at odd positions of ct's targets, unmatched
// decoys in between) and ct pointing at the matched ones.
struct List {
    std::vector<int>    ct;
    std::vector<double> tim;
    int n() const { return int(ct.size()); }
};
static List makeList(const std::vector<double>& matched)
{
    List L;
    L.tim.push_back(-1.0e6);                 // decoy before everything
    for (double t : matched) {
        L.ct.push_back(int(L.tim.size()));
        L.tim.push_back(t);
        L.tim.push_back(-5.0e5);             // decoy after each matched visit
    }
    return L;
}

static bool allStepsPositive(const std::vector<TimeStep>& g, std::string& why)
{
    for (std::size_t i = 0; i < g.size(); ++i)
        if (not (g[i].dt > 0.0)) { why = "dt <= 0 at point " + std::to_string(i); return false; }
    for (std::size_t i = 1; i < g.size(); ++i)
        if (not (g[i].tim > g[i - 1].tim)) { why = "times not strictly increasing at " + std::to_string(i); return false; }
    return true;
}

static std::vector<TimeStep> build(const List& a, const List& b)
{
    std::vector<TimeStep> g;
    buildTimeGrid(a.n(), a.ct, a.tim, b.n(), b.ct, b.tim, g);
    return g;
}

static void unionOfVisits()
{
    const char* name = "unionOfVisits";
    const List a = makeList({0.0, 2.0, 2.02, 10.0}), b = makeList({1.0, 1.0084, 400.0});
    const auto g = build(a, b);
    const double tim[8] = {-3.0, 0.0, 1.0, 1.0084, 2.0, 2.02, 10.0, 400.0};
    const int    rub[8] = {0, 1, 0, 0, 1, 1, 1, 0};
    const int    rom[8] = {0, 0, 1, 1, 0, 0, 0, 1};
    const double dt[8]  = {3.0, 3.0, 1.0, 0.0084, 0.9916, 0.02, 7.98, 390.0};
    if (g.size() != 8) { fail(name, "size " + std::to_string(g.size())); return; }
    int iL = 0, iR = 0;
    for (int i = 0; i < 8; ++i) {
        if (g[i].tim != tim[i]) { fail(name, "tim at " + std::to_string(i)); return; }
        if (g[i].rubin != bool(rub[i]) or g[i].roman != bool(rom[i])) { fail(name, "flags at " + std::to_string(i)); return; }
        if (std::fabs(g[i].dt - dt[i]) > 1e-12) { fail(name, "dt at " + std::to_string(i)); return; }
        const int wantL = rub[i] ? a.ct[iL++] : -1, wantR = rom[i] ? b.ct[iR++] : -1;
        if (g[i].sq != wantL or g[i].sqR != wantR) { fail(name, "indices at " + std::to_string(i)); return; }
    }
    std::string why;
    if (not allStepsPositive(g, why)) { fail(name, why); return; }
    pass(name);
}

static void coincidentVisits()
{
    const char* name = "coincidentVisits";
    const List a = makeList({5.0, 7.0}), b = makeList({5.0, 6.0});
    const auto g = build(a, b);
    if (g.size() != 4) { fail(name, "size " + std::to_string(g.size())); return; }
    const double tim[4] = {2.0, 5.0, 6.0, 7.0};
    for (int i = 0; i < 4; ++i)
        if (g[i].tim != tim[i]) { fail(name, "tim at " + std::to_string(i)); return; }
    if (not (g[1].rubin and g[1].roman and g[1].sq == a.ct[0] and g[1].sqR == b.ct[0])) { fail(name, "coincident point"); return; }
    if (g[2].rubin or not g[2].roman or g[2].sq != -1 or g[2].sqR != b.ct[1]) { fail(name, "Roman-only point"); return; }
    if (not g[3].rubin or g[3].roman or g[3].sq != a.ct[1] or g[3].sqR != -1) { fail(name, "Rubin-only point"); return; }
    std::string why;
    if (not allStepsPositive(g, why)) { fail(name, why); return; }
    pass(name);
}

static void seasonStartsExact()
{
    const char* name = "seasonStartsExact";
    const List a = makeList({0.0, 1.0, 200.0}), b = makeList({150.0, 150.0084, 150.0168});
    const auto g = build(a, b);
    int firstRoman = -1, nRoman = 0, iRubin200 = -1;
    for (std::size_t i = 0; i < g.size(); ++i) {
        if (g[i].roman) { if (firstRoman < 0) firstRoman = int(i); ++nRoman; }
        if (g[i].rubin and g[i].tim == 200.0) iRubin200 = int(i);
    }
    if (firstRoman < 0 or g[firstRoman].tim != 150.0) { fail(name, "first Roman point not at 150 exactly"); return; }
    if (nRoman != 3 or iRubin200 < 0 or firstRoman + 2 >= iRubin200) { fail(name, "Roman points do not precede 200"); return; }
    if (g[firstRoman + 1].tim != 150.0084 or g[firstRoman + 2].tim != 150.0168) { fail(name, "Roman times"); return; }
    std::string why;
    if (not allStepsPositive(g, why)) { fail(name, why); return; }
    pass(name);
}

static void romanOnlySightline()
{
    const char* name = "romanOnlySightline";
    const std::vector<int> emptyCt;
    const std::vector<double> emptyTim;                // a correct builder never indexes these
    const List b = makeList({306.0, 306.0084});
    std::vector<TimeStep> g;
    buildTimeGrid(0, emptyCt, emptyTim, b.n(), b.ct, b.tim, g);
    if (g.size() != 3) { fail(name, "size " + std::to_string(g.size())); return; }
    const double tim[3] = {303.0, 306.0, 306.0084};
    for (int i = 0; i < 3; ++i) {
        if (g[i].tim != tim[i]) { fail(name, "tim at " + std::to_string(i)); return; }
        if (g[i].rubin or g[i].sq != -1) { fail(name, "Rubin flag/index at " + std::to_string(i)); return; }
    }
    if (g[0].roman or not g[1].roman or not g[2].roman) { fail(name, "Roman flags"); return; }
    std::string why;
    if (not allStepsPositive(g, why)) { fail(name, why); return; }
    pass(name);
}

static void rubinOnlySightline()
{
    const char* name = "rubinOnlySightline";
    const List a = makeList({1.0, 2.0, 3.0});
    const std::vector<int> roCt(3, -1);                // the matcher's fill value; never read
    const std::vector<double> roTim(3, 0.0);
    std::vector<TimeStep> g;
    buildTimeGrid(a.n(), a.ct, a.tim, 0, roCt, roTim, g);
    if (g.size() != 4) { fail(name, "size " + std::to_string(g.size())); return; }
    for (const TimeStep& s : g)
        if (s.roman or s.sqR != -1) { fail(name, "Roman flag or index set"); return; }
    for (int i = 1; i < 4; ++i)
        if (not g[i].rubin or g[i].sq != a.ct[i - 1]) { fail(name, "Rubin visit " + std::to_string(i)); return; }
    std::string why;
    if (not allStepsPositive(g, why)) { fail(name, why); return; }
    pass(name);
}

static void emptySightline()
{
    const char* name = "emptySightline";
    const std::vector<int> ct;
    const std::vector<double> tim;
    std::vector<TimeStep> g(5);                        // stale content must be cleared
    buildTimeGrid(0, ct, tim, 0, ct, tim, g);
    if (not g.empty()) { fail(name, "grid not empty"); return; }
    pass(name);
}

static void everyVisitOnce()
{
    const char* name = "everyVisitOnce";
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> U(0.0, 3470.0);
    auto draw = [&](int n) {
        std::vector<double> t;
        while (int(t.size()) < n) {                    // sorted-unique
            t.push_back(U(rng));
            if (int(t.size()) == n) { std::sort(t.begin(), t.end()); t.erase(std::unique(t.begin(), t.end()), t.end()); }
        }
        return t;
    };
    const std::vector<double> tl = draw(2000), tr = draw(50000);
    const List a = makeList(tl), b = makeList(tr);
    const auto g = build(a, b);
    int nL = 0, nR = 0;
    for (const TimeStep& s : g) {
        if (s.rubin) { if (s.sq != a.ct[nL] or a.tim[s.sq] != s.tim) { fail(name, "Rubin visit not at its own time"); return; } ++nL; }
        if (s.roman) { if (s.sqR != b.ct[nR] or b.tim[s.sqR] != s.tim) { fail(name, "Roman visit not at its own time"); return; } ++nR; }
    }
    if (nL != a.n() or nR != b.n()) { fail(name, "visit counts " + std::to_string(nL) + ", " + std::to_string(nR)); return; }
    std::string why;
    if (not allStepsPositive(g, why)) { fail(name, why); return; }
    pass(name);
}

int main()
{
    unionOfVisits();
    coincidentVisits();
    seasonStartsExact();
    romanOnlySightline();
    rubinOnlySightline();
    emptySightline();
    everyVisitOnce();
    std::printf("%s\n", nfail == 0 ? "all held" : "FAILURES");
    return nfail == 0 ? 0 : 1;
}
