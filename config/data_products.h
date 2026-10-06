#ifndef DATA_PRODUCTS_H
#define DATA_PRODUCTS_H

// Numbers that DESCRIBE THE DATA FILES on disk: row counts and catalogue mean masses.
// These are measurements of the data, not choices. A later step will make a script generate
// this file; for now they are the current values and must be updated by hand whenever the
// corresponding data file is regenerated (the read guards in Bulge_LSST.cpp fire on a mismatch).
// Model and survey choices live in config/parameters.h instead.

constexpr int Na = 96;     //rows in "sigmaA_LSST.txt"
constexpr int NaRoman = 123;  // rows in sigma_roman.txt

// CMD/components/*.dat row counts (CMD_BESANCON: ThinDisk, Bulge, ThickDisk, Halo). The lists are the
// COMPLETE Besancon population (no visibility filter; the bulge is a 3.5M random subsample of
// 7,256,344). Since bos10 (Deviation 88) they come from the noise-free catalogue bos10 with Besancon's
// own stellar types; Deviation 81 had built them from bos9, and those lists (889406 / 3500000 / 1058765 /
// 5025 rows) are kept in CMD/components_v2_bos9dev82/. CMD/components/provenance.txt has the counts.
constexpr int N1 = 1288584, N2 = 3500000, N3 = 2008646, N4 = 18541;

// Mean stellar mass of each POPULATION (provenance.txt, mean_mass_population: before the bulge subsample;
// dark entries included). Disk_model's star count Nstart = rho / <m>, so these make Nstart count exactly
// the population a draw comes from. bos10 values (Deviation 88); they replace the bos9 values
// 0.4212 / 0.4199 / 0.4594 / 0.3774 of Deviation 81 and, before that, the legacy 0.403445 (thin), 0.4542
// (thick, halo) and 0.308571 (bulge) of an unrecorded "mass_averaged.cpp". Mirrored in
// analysis/galaxy_model.py (MBAR_*) and analysis/besancon_sample.py (MEANMASS): change all three together.
constexpr double MEANMASS_THIN  = 0.3664;
constexpr double MEANMASS_BULGE = 0.4148;
constexpr double MEANMASS_THICK = 0.4849;
constexpr double MEANMASS_HALO  = 0.4224;

// Data rows in BulgeBaseline.dat, EXCLUDING the header. Regenerated 2026-10-01 (Deviation 69:
// every pointing that can image the distance-rule scan region; was 3686 from a box) from
// baseline_v5.1.0_10yrs.db; readbaselineBulge.py prints the value to use here. The
// previous 7373 counted a doubled file (append-mode bug) and read 3687 phantom rows.
constexpr int Nl = 12915;   // Deviation 80: pointings within scan reach + 1.94 deg (was 12308 at + 1.75)

// Data rows in RomanBaseline.dat, EXCLUDING the header. generateRomanBaseline.py prints
// the value to use here; it must be updated whenever the season pattern, cadence or
// mission start day changes, or the read guard in Bulge_LSST.cpp will fire.
// Current: 6 high-cadence seasons (F146 every 12.1 min) + 4 low-cadence (every 5 days),
// on STScI's real alternating spring/fall visibility windows.
constexpr int NlRoman = 302406;

#endif // DATA_PRODUCTS_H
