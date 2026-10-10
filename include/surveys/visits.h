// The two survey visit-list containers: `lsst` (Rubin) and `roman`, holding the per-visit arrays read
// from the Baseline/ files, and for Roman the photometric error lookup table.
#ifndef ROMAN_SURVEYS_VISITS_H
#define ROMAN_SURVEYS_VISITS_H

#include "common.h"

struct lsst {
    std::vector<int> filter;    // Nl
    
    std::vector<int> ct;        // Nl -- one slot per Rubin visit; see matchVisibleEpochs
    std::vector<double> RA;     // Nl
    std::vector<double> DEC;    // Nl
    std::vector<double> l;      // Nl
    std::vector<double> b;      // Nl
    std::vector<double> tim;    // Nl
    std::vector<double> sig5;   // Nl
    std::vector<double> dist;   // Nl
    std::vector<double> rot;    // Nl -- OpSim rotSkyPos [deg]
    std::vector<double> fwhm;   // Nl -- geometric PSF FWHM of the visit [arcsec] (from seeingFwhmEff)

    //ID  RA  Dec  l  b  start  filter  airmass  seeing  skyBrightness visittime sigma5 targetname distance
    lsst()
        : filter(Nl),
          ct(Nl),
          RA(Nl), DEC(Nl), l(Nl), b(Nl), tim(Nl), sig5(Nl), dist(Nl), rot(Nl), fwhm(Nl)
    {}
};

struct roman {
    std::vector<double> mag;   // NaRoman: mag-vs-error lookup (sigma_roman.txt)
    std::vector<double> err;   // NaRoman

    // Per-visit epoch bookkeeping, mirrors lsst's fields.
    std::vector<int> ct;        // NlRoman -- visible-epoch indices for the current sightline.
                                // Must be the full visit count: a sightline inside a GBTDS field
                                // matches ~51,500 visits and truncation silently ends the mission early.
    std::vector<double> RA;     // NlRoman
    std::vector<double> DEC;    // NlRoman
    std::vector<double> l;      // NlRoman
    std::vector<double> b;      // NlRoman
    std::vector<double> tim;    // NlRoman
    std::vector<double> sig5;   // NlRoman -- only needed if the photometric error varies per visit
    std::vector<int> field;     // NlRoman -- GBTDS field index, 0-4 contiguous block, 5 GC
    std::vector<int> layout;    // NlRoman -- roll of this visit, 0 spring / 1 autumn; selects
                                // which detector layout is placed at (l, b)

    // No `filter` array: only F146 (filter index 6) is modeled for Roman.

    roman()
        : mag(NaRoman), err(NaRoman),
          ct(NlRoman),
          RA(NlRoman), DEC(NlRoman), l(NlRoman), b(NlRoman), tim(NlRoman), sig5(NlRoman),
          field(NlRoman), layout(NlRoman)
    {}
};

#endif // ROMAN_SURVEYS_VISITS_H
