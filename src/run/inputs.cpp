// Reading the input files. Moved out of main() verbatim; the order of the printed lines is kept.
#include "run/inputs.h"

GbtdsLayout loadGbtdsLayout() {
    // The GBTDS detector layout (Deviation 69), read first: the footprint grid is checked
    // against the real detector size, and the scan region is built from the field outline.
    const GbtdsLayout gl = readGbtdsLayout();
    std::cout << "**** GBTDS layout read: " << GBTDS_NSCA << " detectors x " << GBTDS_NLAYOUT
              << " rolls, detector side " << gl.scaSide << " deg, field reach " << gl.rField
              << " deg ****\n";
    return gl;
}

int readRubinVisits(lsst& ls) {
    // --------------------- Read BulgeBaseline.dat ------------------
    std::ifstream fil(PATH_BULGE_BASELINE);
    if (!fil) {
        std::cerr << "Cannot read BulgeBaseline.dat\n";
        return 1;
    }

    int ID;
    double TV, airm, seeingVal, skyB, texp;
//    int tmpFilter;
    std::string header;

    std::getline(fil, header);   // Skip the header line
    for (int i = 0; i < Nl; ++i) {
        fil >> ID >> ls.RA[i] >> ls.DEC[i] >> ls.l[i] >> ls.b[i]
            >> ls.tim[i] >> ls.filter[i] >> airm >> seeingVal >> skyB
            >> TV >> ls.sig5[i] >> texp >> ls.dist[i] >> ls.rot[i];
//if (i == 0) cout << ID << endl;
        // A failed extraction is a silent no-op that leaves this row zero-initialised,
        // and zeros pass every CHECK below: (l,b)=(0,0) is inside the bulge region,
        // tim=0 is inside [0,Tobs], filter=0 is a valid u-band index. So the stream
        // state is the only thing that can catch a short or malformed baseline.
        if (!fil) {
            std::cerr << "BulgeBaseline.dat: read failed at row " << i << " of " << Nl
                      << ". File has fewer rows than Nl, or contains a stray header "
                      << "(append-mode duplicate). Regenerate with readbaselineBulge.py.\n";
            return 1;
        }
        CHECK(airm >= 0.0);
        CHECK(ls.filter[i] >= 0);
        CHECK(ls.filter[i] < 6);
        CHECK(ls.tim[i] >= 0.0);
        CHECK(ls.tim[i] <= Tobs);
        CHECK(texp >= 0.0);
   }
    fil.close();
    std::cout << "**** File BulgeBaseline.dat was read ****\n";
    return 0;
}

int readRubinAstromTable(lsst& ls) {
    std::ifstream fil;
    // --------------------- Read sigmaA_LSST.txt -------------------
    // Atrometric error?
    fil.open(PATH_SIGMA_A_LSST);
    if (!fil) { std::cerr << "Cannot read sigmaA_LSST.txt\n"; return 1; }

    for (int i = 0; i < Na; ++i) {
        fil >> ls.mag[i] >> ls.err[i];

        CHECK(ls.mag[i] >= 16.0);
        CHECK(ls.mag[i] <= 25.0);
        CHECK(ls.err[i] >= 0.2);
        CHECK(ls.err[i] <= 5.0);
    }
    fil.close();
    std::cout << "**** File sigmaA_LSST.txt was read ****\n";
    return 0;
}

int readRomanErrorTable(roman& ro) {
    std::ifstream fil;
    // --------------------- Read sigma_Roman.txt ---------------------
    // Photometric error
    fil.open(PATH_SIGMA_ROMAN);
    if (!fil) { std::cerr << "Cannot read sigma_roman.txt\n"; return 1; }

    for (int i = 0; i < NaRoman; ++i) {
        fil >> ro.mag[i] >> ro.err[i];
        if (!fil) { std::cerr << "sigma_roman.txt: read failed at row " << i << "\n"; return 1; }
        CHECK(ro.mag[i] > 12.0);
        CHECK(ro.err[i] > 0.0);
        if (i > 0) CHECK(ro.mag[i] > ro.mag[i - 1]);
    }
    fil.close();
    {
        // Deviation 72: anchor the curve's 5-sigma point to ROMAN_DEPTH5_AB (see config/parameters.h).
        const double e5 = 1.0857 / 5.0;
        double m5 = -1.0;
        for (int i = 1; i < NaRoman; ++i)
            if (ro.err[i - 1] < e5 and ro.err[i] >= e5) {
                const double f = (std::log(e5) - std::log(ro.err[i - 1]))
                               / (std::log(ro.err[i]) - std::log(ro.err[i - 1]));
                m5 = ro.mag[i - 1] + f * (ro.mag[i] - ro.mag[i - 1]);
                break;
            }
        if (m5 < 0.0) { std::cerr << "sigma_roman.txt never reaches 5 sigma\n"; return 1; }
        const double shift = ROMAN_DEPTH5_AB - m5;
        for (int i = 0; i < NaRoman; ++i) {
            const double phot2 = std::max(ro.err[i] * ro.err[i] - ROMAN_PHOT_FLOOR * ROMAN_PHOT_FLOOR, 0.0);
            ro.mag[i] += shift;
            ro.err[i]  = std::sqrt(phot2 + ROMAN_PHOT_FLOOR * ROMAN_PHOT_FLOOR);
        }
        std::cout << "**** File sigma_roman.txt was read: Penny+2019 5-sigma point " << m5
                  << " AB shifted by " << shift << " mag to " << ROMAN_DEPTH5_AB << " ****\n";
    }
    return 0;
}

int readRomanVisits(roman& ro) {
    std::ifstream fil;
    int ID;
    std::string header;
    // --------------------- Read RomanBaseline.dat -------------------
    // Written by Baseline/generateRomanBaseline.py (adopted GBTDS layout, Deviation 69):
    // ID RA DEC l b time sig5 field layout.
    fil.open(PATH_ROMAN_BASELINE);
    if (!fil) {
        std::cerr << "Cannot read RomanBaseline.dat\n";
        return 1;
    }
    std::getline(fil, header); // skip header line
    for (int i = 0; i < NlRoman; ++i) {
        fil >> ID >> ro.RA[i] >> ro.DEC[i] >> ro.l[i] >> ro.b[i] >> ro.tim[i] >> ro.sig5[i]
            >> ro.field[i] >> ro.layout[i];
        // Same silent zero-fill failure mode as the Rubin read above.
        if (!fil) {
            std::cerr << "RomanBaseline.dat: read failed at row " << i << " of " << NlRoman
                      << ". Regenerate with generateRomanBaseline.py and update NlRoman.\n";
            return 1;
        }

        CHECK(ro.tim[i] >= 0.0);
        CHECK(ro.tim[i] <= Tobs);
        CHECK(ro.layout[i] >= 0 and ro.layout[i] < GBTDS_NLAYOUT);
        CHECK(ro.field[i] >= 0 and ro.field[i] < 6);
        // Deliberately no filter CHECK here — Roman is single-band (F146, index 6) for now.
    }
    fil.close();
    std::cout << "**** File RomanBaseline.dat was read ****\n";
    return 0;
}

int buildRomanSeasons(const roman& ro, RomanSchedule& sched) {
    // --------------------- Roman season geometry (Step D1) ---------------
    // Recovered from the epoch times just read, not restated from the generator.
    // Needed per event to place t0 relative to Roman's observing windows -- the
    // independent variable of the gap-filling result.
    sched = buildRomanSchedule(ro);

    // The clustering is only meaningful if the two spacing populations it separates
    // are genuinely separated: every within-season spacing below the threshold, every
    // between-season gap above it. If a future schedule ever samples a season more
    // sparsely than SEASON_GAP_MIN_DAYS, or packs seasons closer together than it, the
    // seasons come out wrong while dt_edge and t0zone still look perfectly reasonable.
    // Refuse to run rather than emit gap geometry that is quietly fiction.
    if (sched.seasons.size() < 2
        or sched.maxInSeasonSpacing >= SEASON_GAP_MIN_DAYS
        or sched.minSeasonGap       <= SEASON_GAP_MIN_DAYS
        or sched.minSeasonLength    <= 0.0) {
        std::cerr << "FATAL: cannot separate Roman observing seasons from inter-season gaps.\n"
                  << "       SEASON_GAP_MIN_DAYS = " << SEASON_GAP_MIN_DAYS << " d, but this "
                  << "schedule has\n"
                  << "       max in-season spacing " << sched.maxInSeasonSpacing << " d and "
                  << "min inter-season gap " << sched.minSeasonGap << " d,\n"
                  << "       shortest season " << sched.minSeasonLength << " d\n"
                  << "       (" << sched.seasons.size() << " season(s) found over days "
                  << sched.missionStart << " - " << sched.missionEnd << ").\n"
                  << "       Retune SEASON_GAP_MIN_DAYS in config/parameters.h against the cadence in\n"
                  << "       Baseline/generateRomanBaseline.py before trusting dt_edge/t0zone.\n";
        return 2;
    }
    std::cout << "**** Roman schedule: " << sched.seasons.size() << " seasons over days "
              << sched.missionStart << " - " << sched.missionEnd << " ****\n";
    return 0;
}

int readSkyTables(extin& ex, const lsst& ls) {
    // --------------------- Read extinction ------------------------
    readExtinction(ex, PATH_EXT_TABLES);
    readLensML(PATH_LENS_ML);   // luminous lenses (Deviation 74)
    readLsstCamMap(PATH_LSSTCAM_FOV);   // Rubin's footprint (Deviation 80)
    {
        // The galactic->ICRS conversion the coverage test relies on, checked against OpSim's own
        // RA/Dec for every Rubin visit (whose l, b readbaselineBulge.py derived with astropy).
        double worst = 0.0;
        for (int i = 0; i < Nl; ++i) {
            double ra, dec;
            galToIcrs(ls.l[i], ls.b[i], ra, dec);
            double dra = std::fabs(ra - ls.RA[i]); if (dra > 180.0) dra = 360.0 - dra;
            worst = std::max(worst, std::hypot(dra * std::cos(dec * M_PI / 180.0), dec - ls.DEC[i]));
        }
        std::cout << "galactic->ICRS check over " << Nl << " Rubin visits: worst " << worst * 3600.0 << " arcsec\n";
        if (worst > 1e-4) { std::cerr << "ERROR: galToIcrs disagrees with OpSim by " << worst << " deg\n"; return 1; }
    }
    return 0;
}

int readCmdTables(CMD& cm) {
    // --------------------- Call read_cmd --------------------------
    read_cmd(cm);
    std::cout << "******* read_cmd was done ************" << std::endl;
    return 0;
}

int loadInputs(lsst& ls, roman& ro, extin& ex, CMD& cm, RomanSchedule& sched) {
    if (int rc = readRubinVisits(ls))        return rc;
    if (int rc = readRubinAstromTable(ls))   return rc;
    if (int rc = readRomanErrorTable(ro))    return rc;
    if (int rc = readRomanVisits(ro))        return rc;
    if (int rc = buildRomanSeasons(ro, sched)) return rc;
    if (int rc = readSkyTables(ex, ls))      return rc;
    if (int rc = readCmdTables(cm))          return rc;
    return 0;
}
