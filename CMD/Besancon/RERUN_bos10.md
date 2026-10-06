# Besancon re-run: parameters for bos10 (replacement catalogue) and bos11a-e (validation cones)

**STATUS 2026-10-05: run by the user; files `bos10`, `bos11b-e` (+ `-head`) are here and passed the
acceptance checks (DEVIATIONS 84). `bos11a` is EMPTY: its V range was entered as `-99.00--99.00`;
re-run it with `-99.00- 99.00` -- done: `bos11a-redo`. bos10 has 4 corrupted rows (field count != 38) that readers must skip.**

Written 2026-10-04 from the bos9 audit (OPEN_ITEMS "BESANCON SAMPLE AUDIT", DEVIATIONS 83).
The blocks below are in the layout of `bos9-head.dat`, so each line maps to one field of the
Besancon web form. Lines marked `<-- CHANGED` differ from bos9; everything else is copied from
bos9-head.dat unchanged. Model version must stay **1612** (the header's first block, the evolution
model, density laws and metallicity table are the model's, not form inputs: keep them as in bos9).

## Why a re-run (what bos9 gets wrong for our use)

1. **Injected parameter noise.** bos9 was generated with Besancon's default "errors" on the stellar
   parameters, and the OUTPUT values are perturbed by them: Mass by sigma = 0.10 Msun (bulge red
   giants, a single-age population, show a mass spread of exactly 0.101), Age by 0.10 Gyr (every
   single-age population has an age sd of exactly 0.100), Teff 50 K, logg 0.1, [M/H] 0.1, [a/Fe] 0.02,
   and V/B/U/I/K by 0.02 mag. The catalogue mass is therefore a noisy label: lens_ml.dat and the
   Deviation 82 shift, which bin by it, are biased. Stars whose noisy mass came out negative were
   dropped, and half the halo (age 14 Gyr + noise > the 14-Gyr range) was dropped too.
   **Fix: every error law = 0.**
2. **V <= 29 apparent-magnitude limit.** With Av = 0 this cuts M_V > 29 - DM(d). It does not touch
   the bulge (its main sequence ends at M_V ~ 12.9), but it truncates the thin disc beyond a few kpc
   (fraction of thin-disc stars with M_V > 13 falls from 0.45 within 1 kpc to 0.28 at 7.5-8 kpc) and
   every component's white dwarfs. **Fix: no apparent-magnitude limit.**
3. **Colour cuts** B-V, U-B, V-I, V-K >= -1.00 remove the hottest stars (U-B reaches exactly -1.000
   in bos9). Negligible for the bulge, but there is no reason to keep them. **Fix: -99 to 99.**
4. **Distance range 0-8 kpc** covers only the near half of the bulge. Harmless for the bulge's
   stellar types (its luminosity function is identical in every distance bin in bos9), but a run to
   the simulator's own MaxD = 12 kpc lets the Besancon density along the whole sightline be compared
   with the C++ `Disk_model` (they already disagree by factors 0.5-2 per component at 0-8 kpc).
   **Fix: 0-12 kpc, with half the solid angle to keep the file near bos9's size.**
5. **Age range upper bound.** With zero errors the halo is exactly 14 Gyr and the range is written
   as half-open `[0, 14[`; widen to 15 so the halo is not dropped. Same reasoning for the [a/Fe]
   and AbsMag ranges (widened only for safety).

Extinction stays OFF (Av = 0): the simulator applies its own dust at each star's simulated
distance. Kinematics stay off (the simulator draws its own). Photometric system unchanged (our
pipeline uses Teff, logg, Mbol, [M/H], [a/Fe] with MIST bolometric corrections, not Besancon's bands).

Expected size of bos10: ~6 M bulge + ~2 M thin + ~1 M thick rows, i.e. ~3 GB, similar to bos9.

## Run A -- bos10 (the replacement catalogue). REQUIRED.

```
  (l =     0.5; b =   -1.4)  Distance  0.0 a  12.0 kpc;  Solid angle      0.050  square degree   <-- CHANGED (distance 0-12, solid angle 0.05)
  Adaptive step in distance  (20,30,50)
      Apparent magnitudes Band V    :-99.00- 99.00                                                <-- CHANGED (was 10.00- 29.00)
                          Band B    :-99.00- 99.00
                          Band U    :-99.00- 99.00
                          Band R    :-99.00- 99.00
                          Band I    :-99.00- 99.00
                          Band J    :-99.00- 99.00
                          Band H    :-99.00- 99.00
                          Band K    :-99.00- 99.00
                          Band L    :-99.00- 99.00

      Colours:  B-V        -99.00 -  99.00                                                        <-- CHANGED (was -1.00)
                U-B        -99.00 -  99.00                                                        <-- CHANGED (was -1.00)
                V-I        -99.00 -  99.00                                                        <-- CHANGED (was -1.00)
                V-K        -99.00 -  99.00                                                        <-- CHANGED (was -1.00)

      AbsMag   range [      -10.00-       30.00[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)   <-- CHANGED (range; was -7/24)
      Teff     range [      300.00-   100000.00[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)   <-- CHANGED (error was 50)
      Logg     range [       -1.00-        9.50[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)   <-- CHANGED (error was 0.1)
      Met      range [       -4.20-        1.50[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)   <-- CHANGED (error was 0.1)
      [Alp/Fe] range [       -1.00-        1.00[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)   <-- CHANGED (range; error was 0.02)
      Age      range [        0.00-       15.00[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)   <-- CHANGED (range; error was 0.1)
      Mass     range [        0.00-       90.00[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)   <-- CHANGED (error was 0.1)
      Selection in age: from .00E9Gyr to OldTh
      Spectral types: O0.0 - D5.0  Luminosity class    :
      No kinematics
   Computation mode :           1
      Poisson drawings (ran1) for N <=  6
      Photometric error laws: V   : 0.000000E+00+mag*0.000000E+00+mag^2*0.000000E+00             <-- CHANGED (was 0.02)
                              B   : 0.000000E+00+mag*0.000000E+00+mag^2*0.000000E+00             <-- CHANGED (was 0.02)
                              U   : 0.000000E+00+mag*0.000000E+00+mag^2*0.000000E+00             <-- CHANGED (was 0.02)
                              I   : 0.000000E+00+mag*0.000000E+00+mag^2*0.000000E+00             <-- CHANGED (was 0.02)
                              K   : 0.000000E+00+mag*0.000000E+00+mag^2*0.000000E+00             <-- CHANGED (was 0.02)
   Atmosphere models: Basel 3.1
      Rv :   3.100
      System:   V      B      U      R      I      J      H      K      L      12
      Al/Av:    1.009  1.332  1.555  0.841  0.600  0.293  0.184  0.118  0.000  0.000
      Creation of a catalogue of simulated stars
   ...
   INTERSTELLAR MATTER AND EXTINCTION
      Interstellar matter: Einasto law: excentricity and scale lengths in pc 0.0140  7222.0  5777.0
      Local interstellar matter ( M/pc3) :.0500
      Absorption at s.n. : 0.00 mag./kpc   0 clouds                                                (unchanged: NO extinction)
```

If the form refuses an open V range, use `Band V : -99.00- 45.00` (any limit fainter than M_V = 30
at 12 kpc, i.e. V > 45.4, removes the cut). If it refuses AbsMag `[-10, 30[`, keep bos9's `[-7, 24[`
(it only removes the coolest brown dwarfs, which are dark in our pipeline anyway).

**After bos10 arrives, before it replaces bos9:** check (1) bulge giants' mass sd is ~0 (it was
0.101); (2) single-age populations have age sd 0; (3) the thin-disc fraction of M_V > 13 is flat in
distance; (4) N rows per component. Then rebuild the npy cache (`analysis/besancon_sample.py`), and
re-decide lens_ml.dat and Deviation 82 on noise-free masses (OPEN_ITEMS audit items B2, B3) before
re-running `CMD/lens_ml_table.py` and `CMD/BolometricCorrection.py`.

## Run B -- bos11a-e (validation cones). OPTIONAL, but answers "is one sightline enough?"

Same block as Run A, changing only the first line. Purpose: (1) does each component's luminosity
function change across the footprint and Rubin scan region (thin/thick age mixes depend on height
above the plane; the bulge should not change by construction)? (2) Besancon's star counts along
each line, for comparison with the C++ density model and with observed counts at two classical
fields (Baade's window: Holtzman 1998, Zoccali 2000, OGLE; SWEEPS: Calamida 2015).

```
bos11a  (l =    -3.0; b =   -1.4)  Distance  0.0 a  12.0 kpc;  Solid angle      0.010  square degree    west end of the scan region
bos11b  (l =     3.5; b =   -1.4)  Distance  0.0 a  12.0 kpc;  Solid angle      0.010  square degree    east end (near side of the bar)
bos11c  (l =     0.5; b =   -0.4)  Distance  0.0 a  12.0 kpc;  Solid angle      0.005  square degree    near the plane (densest)
bos11d  (l =     1.06; b =  -3.81) Distance  0.0 a  12.0 kpc;  Solid angle      0.010  square degree    Baade's window
bos11e  (l =     1.25; b =  -2.65) Distance  0.0 a  12.0 kpc;  Solid angle      0.010  square degree    SWEEPS field (HST)
```

Expected sizes: ~0.3-1.5 M rows each (~0.1-0.5 GB).
