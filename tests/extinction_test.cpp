// Unit test for the CCM89 reddening law (helper.cpp AlAv), added with Deviation 53 after the
// law had run inverted for two months without anything noticing. Needs no data files.
//
// Pinned: (1) A_V/A_V = 1 at V (0.549 um) for every R_V the populations use -- the law's
// definition, so it cannot pass with the wavelength handled wrongly; (2) extinction falls
// monotonically from u to F146; (3) A_lambda/A_V at the seven survey bands for R_V = 2.5,
// against an independent Python evaluation of the CCM89 polynomials (O'Donnell 1994 is NOT used;
// the code implements the original 1989 optical coefficients). Exit 0 = all held.
#include "Bulge.h"
#include <cstdio>
#include <cmath>

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
    std::printf("%s\n", fail ? "extinction_test: FAILED" : "extinction_test: all held");
    return fail ? 1 : 0;
}
