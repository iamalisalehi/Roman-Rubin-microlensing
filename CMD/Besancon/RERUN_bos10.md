# Requesting the Besancon catalogue

The simulator draws its source and lens stars (kind, mass, age, metallicity, temperature, gravity,
bolometric magnitude) from a synthetic catalogue made with the **Besancon Galaxy Model**. The model
can only be run through its web form, so this file is the recipe: which options to choose and why.
The production catalogue is called `bos10`; any catalogue made the same way can replace it.

## Quick recipe

1. Open <https://model.obs-besancon.fr/> and choose **model 1612** (the version every production run
   used; another version works but gives different numbers). Ask for a **catalogue of simulated stars**,
   not star counts.
2. Enter **one field**: (l, b) = (0.5, -1.4) deg, distances 0 to 12 kpc, solid angle 0.05 deg^2.
   The simulator uses this one catalogue for every sightline: it takes the *kinds* of stars from it and
   computes how many there are, and where they sit, from its own density model.
3. Set the options in the table below. Leave everything else at the model's defaults.
4. Download the result: a text file of roughly 10 million stars, about 3.7 GB.
5. Put it anywhere and point `BESANCON_CATALOGUE` in your pipeline config at it (the default is
   `CMD/Besancon/bos10`). `pipeline/pipeline.sh <config> check` then verifies the file's header.

| Form option | Value | Why |
|---|---|---|
| Field and distance | (0.5, -1.4), 0-12 kpc, 0.05 deg^2 | 12 kpc is the simulator's maximum source distance |
| Extinction | **off**: 0 mag/kpc, 0 clouds | the simulator applies its own 3D dust at each star's distance |
| Kinematics | **off** | the simulator draws its own velocities |
| Apparent magnitudes, every band | -99 to 99 (no limit) | a V limit quietly removes the faint dwarfs and white dwarfs of the disc |
| Colours (B-V, U-B, V-I, V-K) | the widest the form allows (-10 to 99) | colour cuts remove the hottest stars |
| Absolute magnitude | [-10, 30[ | keeps the faintest stars |
| Every "error" law (Teff, logg, [M/H], [a/Fe], age, mass, and the photometric errors V B U I K) | **0** | see below |
| Age | [0, 14[ Gyr | the whole model age range |
| Mass | [0, 90[ Msun | |
| Spectral types | O0.0 to D5.0 | includes white dwarfs |

**Why the errors must be zero.** By default the model adds random "observational" noise to every output
value: mass by 0.10 Msun, age by 0.1 Gyr, Teff by 50 K, and so on. The simulator uses the catalogue mass
as the star's true mass (it bins lens masses by it and draws lenses from it), so noisy masses bias the
results. Noise also pushes some masses below zero (those stars are dropped) and some halo ages past the
age range (those stars are dropped too).

**If the form refuses a value:** for an open V range use -99 to 45 (anything fainter than V = 45 is
no limit at 12 kpc). For the absolute magnitude, [-7, 24[ is acceptable: it only loses the coolest
brown dwarfs, which the simulator treats as dark anyway.

## What the file must look like

The usual Besancon text output. The line before the data is a `#` header with the 38 column names, and
it must contain `Teff logg Pop Age Mass Mbol [M/H] [a/Fe] CL Typ`. The pipeline's `check` and `prep`
stages test this before doing any work. A few damaged rows (a wrong number of fields) are skipped by
the star-list builder (`bos10` has 4).

## The reference header

For comparison, this is the parameter block at the top of `bos10`, the production catalogue. Your
file's block should match it except for anything you changed on purpose. (The block also lists the
model's own evolution tracks, density laws and metallicities. Those belong to model 1612 and are not
form inputs.)

```
  (l =     0.5; b =   -1.4)  Distance  0.0 a  12.0 kpc;  Solid angle      0.050  square degree
  Adaptive step in distance  (20,30,50)
      Apparent magnitudes Band V    :-99.00- 99.00
                          Band B    :-99.00- 99.00
                          Band U    :-99.00- 99.00
                          Band R    :-99.00- 99.00
                          Band I    :-99.00- 99.00
                          Band J    :-99.00- 99.00
                          Band H    :-99.00- 99.00
                          Band K    :-99.00- 99.00
                          Band L    :-99.00- 99.00

      Colours:  B-V        -10.00 -  99.00
                U-B        -10.00 -  99.00
                V-I        -10.00 -  99.00
                V-K        -10.00 -  99.00

      AbsMag   range [      -10.00-       30.00[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)
      Teff     range [      300.00-   100000.00[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)
      Logg     range [       -1.00-        9.50[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)
      Met      range [       -4.20-        1.50[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)
      [Alp/Fe] range [       -0.10-        0.70[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)
      Age      range [        0.00-       14.00[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)
      Mass     range [        0.00-       90.00[. Errors as a fct of magnitude:  0.0000+  0.0000*mag +  0.0000*mag^2)
      Selection in age: from .00E9Gyr to OldTh
      Spectral types: O0.0 - D5.0  Luminosity class    :
      No kinematics
   Computation mode :           1
      Poisson drawings (ran1) for N <=  6
      Photometric error laws: V   : 0.000000E+00+mag*0.000000E+00+mag^2*0.000000E+00
                              B   : 0.000000E+00+mag*0.000000E+00+mag^2*0.000000E+00
                              U   : 0.000000E+00+mag*0.000000E+00+mag^2*0.000000E+00
                              I   : 0.000000E+00+mag*0.000000E+00+mag^2*0.000000E+00
                              K   : 0.000000E+00+mag*0.000000E+00+mag^2*0.000000E+00
   Atmosphere models: Basel 3.1
   ...
   INTERSTELLAR MATTER AND EXTINCTION
      ...
      Absorption at s.n. : 0.00 mag./kpc   0 clouds
```

## Optional: other fields, for checks

The simulator needs only the one catalogue. Smaller catalogues of other fields, made with exactly the
same options, are useful to check that one field is representative (does each component's mix of stars
change across the survey area?) and to compare the model's star counts with observed fields. The ones
made for this project, 0.005-0.01 deg^2 each (0.1-0.6 GB):

| Name | (l, b) deg | Where |
|---|---|---|
| bos11a | (-3.0, -1.4) | west end of the Rubin scan region |
| bos11b | (3.5, -1.4) | east end, the near side of the bar |
| bos11c | (0.5, -0.4) | close to the plane, the densest field |
| bos11d | (1.06, -3.81) | Baade's window |
| bos11e | (1.25, -2.65) | the HST SWEEPS field |
