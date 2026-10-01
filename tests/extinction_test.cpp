// Unit test for the CCM89 reddening law (helper.cpp AlAv), added with Deviation 53 after the
// law had run inverted for two months without anything noticing. Needs no data files.
//
// Pinned: (1) A_V/A_V = 1 at V (0.549 um) for every R_V the populations use -- the law's
// definition, so it cannot pass with the wavelength handled wrongly; (2) extinction falls
// monotonically from u to F146; (3) A_lambda/A_V at the seven survey bands for R_V = 2.5,
// against an independent Python evaluation of the CCM89 polynomials (O'Donnell 1994 is NOT used;
// the code implements the original 1989 optical coefficients). Exit 0 = all held.
//
// Deviation 70 added the extinction TABLE reader, readExtinction(): a small fixture file is read
// and its interpolation and nearest-table choice checked; then five malformed files (a NaN, a
// decreasing profile, a short row, an extra value, a wrong row count) must each make it exit
// non-zero -- run in a forked child, since refusing means exiting. The old reader accepted all of
// these silently.
#include "Bulge.h"
#include <cstdio>
#include <cmath>
#include <sys/wait.h>
#include <unistd.h>

static const char* FIX = "extinction_test_fixture.tmp";

static void writeFixture(const std::string& body, int nTables) {
    std::ofstream f(FIX);
    f << "# ext_tables v1 -- test fixture\n# k 0.08\n# n_tables " << nTables
      << "\n# n_dist 3\n# dist 1.0 2.0 4.0\n" << body;
}

static bool readerRefuses(const std::string& body, int nTables) {
    writeFixture(body, nTables);
    std::fflush(nullptr);
    pid_t pid = fork();
    if (pid == 0) {
        std::freopen("/dev/null", "w", stdout);
        std::freopen("/dev/null", "w", stderr);
        extin ex;
        readExtinction(ex, FIX);
        std::_Exit(0);                         // accepted: the test fails
    }
    int status = 0;
    waitpid(pid, &status, 0);
    return WIFEXITED(status) and WEXITSTATUS(status) != 0;
}

int main()
{
    int fail = 0;
    auto check = [&](bool ok, const char* what, double got, double want) {
        std::printf("%-4s %-34s got %.4f  want %.4f\n", ok ? "ok" : "FAIL", what, got, want);
        if (!ok) ++fail;
    };

    const double Rvs[] = {2.5, 3.1};
    for (double Rv : Rvs) {
        char buf[64]; std::snprintf(buf, sizeof buf, "A_V/A_V at 0.549 um, R_V=%.1f", Rv);
        double v = AlAv(0.549, Rv);
        check(std::fabs(v - 1.0) < 0.01, buf, v, 1.0);
    }

    // Survey bands in Bulge.h order: ugrizy, F146.
    const double lam[7]   = {0.367, 0.482, 0.622, 0.755, 0.869, 0.971, 1.464};
    const char*  name[7]  = {"u", "g", "r", "i", "z", "y", "F146"};
    const double want25[7] = {1.694, 1.216, 0.854, 0.627, 0.460, 0.381, 0.197};
    for (int i = 0; i < 7; ++i) {
        char buf[64]; std::snprintf(buf, sizeof buf, "A/A_V %s, R_V=2.5", name[i]);
        double v = AlAv(lam[i], 2.5);
        check(std::fabs(v - want25[i]) < 5e-3, buf, v, want25[i]);
        if (i > 0) {
            std::snprintf(buf, sizeof buf, "falls from %s to %s", name[i-1], name[i]);
            double prev = AlAv(lam[i-1], 2.5);
            check(v < prev, buf, v, prev);
        }
    }
    // ---- the table reader ----
    writeFixture("0.0 -1.0 1.0 2.0 6.0\n0.5 -1.0 0.5 0.5 1.5\n", 2);
    {
        extin ex;
        readExtinction(ex, FIX);
        check(ex.nTables == 2 and ex.nDist == 3, "reader: 2 tables x 3 distances", ex.nTables * 10 + ex.nDist, 23);
        check(std::fabs(interpExtinctionAlongSightline(ex, 0, 2.0) - 2.0) < 1e-6, "reader: A_V at a grid distance",
              interpExtinctionAlongSightline(ex, 0, 2.0), 2.0);
        check(std::fabs(interpExtinctionAlongSightline(ex, 0, 3.0) - 4.0) < 1e-6, "reader: A_V interpolated",
              interpExtinctionAlongSightline(ex, 0, 3.0), 4.0);
        check(std::fabs(interpExtinctionAlongSightline(ex, 1, 9.0) - 1.5) < 1e-6, "reader: held beyond the grid",
              interpExtinctionAlongSightline(ex, 1, 9.0), 1.5);
        check(std::fabs(interpExtinctionAlongSightline(ex, 0, 0.1) - 1.0) < 1e-6, "reader: held before the grid",
              interpExtinctionAlongSightline(ex, 0, 0.1), 1.0);
        check(nearestSightline(ex, 0.4, -1.1) == 1, "reader: nearest table", nearestSightline(ex, 0.4, -1.1), 1);
    }
    const struct { const char* what; const char* body; int n; } bad[] = {
        {"refuses: a NaN",             "0.0 -1.0 1.0 nan 6.0\n", 1},
        {"refuses: decreasing A_V",    "0.0 -1.0 1.0 3.0 2.0\n", 1},
        {"refuses: a short row",       "0.0 -1.0 1.0 2.0\n", 1},
        {"refuses: an extra value",    "0.0 -1.0 1.0 2.0 3.0 4.0\n", 1},
        {"refuses: wrong row count",   "0.0 -1.0 1.0 2.0 3.0\n", 2},
    };
    for (const auto& c : bad) {
        const bool refused = readerRefuses(c.body, c.n);
        check(refused, c.what, refused, 1);
    }
    std::remove(FIX);

    std::printf("%s\n", fail ? "extinction_test: FAILED" : "extinction_test: all held");
    return fail ? 1 : 0;
}
