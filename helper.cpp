#include "Bulge.h"


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double RandN(double sigm, double nnd) {
    std::normal_distribution<double> normal(0.0, 1.0);

    double x;
    do {
        x = sigm * normal(rng);
    } while (std::abs(x) > sigm * nnd); //[-N sigma:N sigma]

    return x;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double RandR(double down, double up) {
    std::uniform_real_distribution<double> dist(down, up);

    return dist(rng);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//               tE & ML & DL & U0 Functions                          //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
int FunctE(lens & l) {
   int gg = -1;

   if (l.tE <= tE_min)       gg = 0;

   else if (l.tE >= tE_max)  gg = GG;

   else {
       for (int i = 1; i <= GG; ++i) {
          if(double((l.tE - l.tEs[i-1]) * (l.tE - l.tEs[i])) < 0.0 or l.tE == l.tEs[i-1]) { gg = i - 1;  break; }
       }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(l.tE > 0.0);
   CHECK(l.tEs[1] > 0.0);
   CHECK(l.tEs[10] > 0.0);

   return(gg);
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
int FuncMl(lens & l) {
   int gg = -1;

   if (l.Ml <= Ml_min)      gg = 0;
   else if (l.Ml >= Ml_max) gg = GG;
   else {
      for(int i = 1; i <= GG; ++i) {
          if (double((l.Ml - l.Mls[i-1]) * (l.Ml - l.Mls[i])) < 0.0 or l.Ml == l.Mls[i-1]) { gg = i - 1;  break; }
      }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(l.Ml >= 3.0);

   return(gg);
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
int FuncPi(lens & l) {
   int gg = -1;
   double lpi = std::log10(l.pirel);

   if (lpi <= pi_min)      gg = 0;
   else if (lpi >= pi_max) gg = GG;
   else {
      for (int i = 1; i <= GG; ++i) {
          if (double((lpi - l.pis[i-1]) * (lpi - l.pis[i])) < 0.0 or lpi == l.pis[i-1]) { gg = i - 1;  break; }
      }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);

   return(gg);
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
int Funcu0(lens & l){
   int gg = -1;

   if (l.u0 <= u0_min)      gg = 0;
   else if (l.u0 >= u0_max) gg = GG;
   else {
      for (int i = 1; i <= GG; ++i) {
          if (double((l.u0 - l.u0s[i-1]) * (l.u0 - l.u0s[i])) < 0.0 or l.u0 == l.u0s[i-1]) { gg = i - 1;  break; }
      }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(l.u0 >= 0.0);
   CHECK(l.u0 <= u0m);

   return(gg);
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
int FuncMu(lens & l){
   int gg = -1;
   double mur = double(l.murel * year); //mas/years

   if (mur <= mu_min)      gg = 0;
   else if (mur >= mu_max) gg = GG;
   else {
       for (int i = 1; i <= GG; ++i) {
          if (double((mur - l.mus[i-1]) * (mur - l.mus[i])) < 0.0 or mur == l.mus[i-1]) { gg = i - 1;  break; }
       }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(mur >= 0.0);

   return(gg);
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
int FuncMb(lens & l, double gadr) {
   int gg = -1;

   if (gadr <= mb_min)       gg = 0;
   else if (gadr >= mb_max)  gg = GG;
   else {
       for (int i = 1; i <= GG; ++i) {
          if (double((gadr - l.mbs[i-1]) * (gadr - l.mbs[i])) < 0.0 or gadr == l.mbs[i-1]) { gg = i - 1; break; }
       }
   }
   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(gadr >= 0.0);

   return(gg);
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
int FuncFb(lens & l, double blen) {
   int gg = -1;

   if (blen <= fb_min)       gg = 0;
   else if (blen >= fb_max)  gg = GG;
   else {
       for (int i = 1; i <= GG; ++i) {
          if(double((blen - l.fbs[i-1]) * (blen - l.fbs[i])) < 0.0 or blen == l.fbs[i-1]) { gg = i - 1; break;}
       }
   }

   CHECK(gg >= 0);
   CHECK(gg <= GG);
   CHECK(blen >= 0.0);
   CHECK(blen <= 1.0);

   return(gg);
}


///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Error LSST calculations                    //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
double errlsstM(double mag, int fi, double sig5){ //LSST Photometric Error 

    double x, Delta1 = 0.0;
    x = std::pow(10.0, 0.4 * (mag - sig5));
    Delta1 = std::sqrt(std::fabs((0.04 - gama[fi]) * x + gama[fi] * x * x));
    
    if (Delta1 < 0.0001)   Delta1 = 0.0001;

    CHECK(sig5 >= 0.0);
    CHECK(sig5 <= 40.0);
    CHECK(Delta1 >= 0.00001);
    CHECK(x > 0.0);
   
    return std::sqrt(delta2 * delta2 + Delta1 * Delta1);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double errlsstA(lsst & ls, double ghadr){ //LSST Astrometric Error  //Change it!!

    double error = -1.0, shib = 0.0;

    if (ghadr < ls.mag[0] or  ghadr == ls.mag[0]) error = double(ls.err[0]);

    else if (ghadr > ls.mag[Na-1] or ghadr == ls.mag[Na-1]) {
        shib = (ls.err[Na-1] - ls.err[Na-2]) / (ls.mag[Na-1] - ls.mag[Na-2]);
        error = double(ls.err[Na-1] + shib * (ghadr - ls.mag[Na-1]));
    }

    else {
        for (int i = 1; i < Na; ++i) {
            if (double((ghadr - ls.mag[i]) * (ghadr - ls.mag[i-1])) < 0.0 or ghadr == ls.mag[i-1]) {
                shib = (ls.err[i] - ls.err[i-1]) / (ls.mag[i] - ls.mag[i-1]);
                error = double(ls.err[i-1] + shib * (ghadr - ls.mag[i-1]));
                break;
            }
        }
    }

    CHECK(error > 0.0);
    CHECK(error >= ls.err[0]);
    CHECK(ghadr >= 0.0);

    // Renormalise the shipped mission-averaged curve to a PER-VISIT error, which is what
    // l.erra[] means and what FisherM assumes. See the LSST_AST_* block in Bulge.h for the
    // two independent checks that fix the factor at 26.74. Applied here rather than by
    // editing files/sigmaA_LSST.txt so the input data stay as delivered and the correction
    // is visible in the code that depends on it.
    error *= LSST_AST_RENORM;

    CHECK(error >= LSST_AST_FLOOR * 0.999);
    CHECK(std::isfinite(error));

    return(error);
}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
// Roman WFI astrometric error for one F146 exposure, in milliarcseconds (Step H4).
//
// Replaces the errlsstA() placeholder that stood in for Roman -- Rubin's astrometric error
// curve evaluated at Roman's magnitude, which had no reason to be right and was flagged in
// OPEN_ITEMS.md. Constants, their sources and the per-exposure caveat are in Bulge.h.
//
// Three regimes:
//   m <= 20.62   1.1 mas       centroiding floor, 1% of the 110 mas pixel. A systematic,
//                              not photon noise, so it does NOT improve for brighter stars.
//   20.62 - 23.5 rises at 0.333/mag   interpolation between the two published anchors.
//   m >  23.5    rises at 0.4/mag     background dominated, SNR ~ counts.
//
// Unlike errlsstA this reads no data file and needs no instrument struct, so it takes the
// magnitude alone. Its photometric sibling errRomanM() lives in Bulge_LSST.cpp instead,
// because that one needs the roman struct.
double errRomanA(double magF146){

    double error = ROMAN_AST_FLOOR;

    if (magF146 > ROMAN_AST_MBKG) {
        error = ROMAN_AST_SBKG * std::pow(10.0, ROMAN_AST_SLOPE_BKG * (magF146 - ROMAN_AST_MBKG));
    }
    else if (magF146 > ROMAN_AST_MFLR) {
        error = ROMAN_AST_FLOOR * std::pow(10.0, ROMAN_AST_SLOPE_SRC * (magF146 - ROMAN_AST_MFLR));
    }

    // The floor is a floor: the source-dominated branch is continuous with it at
    // ROMAN_AST_MFLR by construction, but clamp anyway so no future edit to the constants can
    // silently return a precision better than Roman can centroid.
    if (error < ROMAN_AST_FLOOR) error = ROMAN_AST_FLOOR;

    CHECK(error >= ROMAN_AST_FLOOR);
    CHECK(std::isfinite(error));

    return(error);
}


///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Read CMD                                   //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void read_cmd(CMD & cm)
{
// mass  logT  Mbol  Age  Pop  Roman_F146  LSST_u  LSST_g  LSST_r  LSST_i  LSST_z  LSST_y  CL  Type //20/07/2026
// 0     1     2     3    4    5           6       7       8       9       10      11      12  13
    double Mbol, Pop;
    double dummy;
    std::string header;
//    double logL, gravity, metal, B, V, R, I, J, H;

    // ================================ THIN DISK =============================
    std::ifstream fp2("./CMD/components/thin_disk.dat");
    if (!fp2.is_open()) {
        throw std::runtime_error("cannot read ./CMD/components/thin_disk.dat");
    }

    std::getline(fp2, header);   // Skip the header line
    for (size_t j = 0; j < N1; ++j) {
        CHECK(cm.mass_thin[j]   >= 0.0);
        CHECK(cm.logT_thin[j]   >= 0.0);
        CHECK(cm.Mab_thin[j][2] <= 20.0);
        CHECK(cm.age_thin[j]    <= 10);
        CHECK(cm.cl_thin[j]     <= 7);
        CHECK(cm.typ_thin[j]    <= 9.0);

          if (!(fp2 >> cm.mass_thin[j]
                    >> cm.logT_thin[j]
                    >> Mbol
                    >> cm.age_thin[j]
                    >> Pop
                    >> cm.Mab_thin[j][6]
                    >> cm.Mab_thin[j][0]
                    >> cm.Mab_thin[j][1]
                    >> cm.Mab_thin[j][2]
                    >> cm.Mab_thin[j][3]
                    >> cm.Mab_thin[j][4]
                    >> cm.Mab_thin[j][5]
                    >> cm.cl_thin[j]
                    >> cm.typ_thin[j])) {
              throw std::runtime_error("Unexpected end of thin_disk.dat");
          }
      }

    // Make sure there's no extra data.
    CHECK(!(fp2 >> dummy));
    fp2.close();

    // ================================ BULGE ==================================
    fp2.open("./CMD/components/bulge.dat");
    if (!fp2.is_open()) {
        throw std::runtime_error("cannot read ./CMD/components/bulge.dat");
    }

    std::getline(fp2, header);   // Skip the header line
    for (size_t j = 0; j < N2; ++j) {
        CHECK(cm.mass_bulge[j]   >= 0.0);
        CHECK(cm.logT_bulge[j]   >= 0.0);
        CHECK(cm.Mab_bulge[j][2] <= 18.0);
        CHECK(cm.age_bulge[j]    <= 10);
        CHECK(cm.cl_bulge[j]     <= 7);
        CHECK(cm.typ_bulge[j]    <= 9.0);
          if (!(fp2 >> cm.mass_bulge[j]
                    >> cm.logT_bulge[j]
                    >> Mbol
                    >> cm.age_bulge[j]
                    >> Pop
                    >> cm.Mab_bulge[j][6]
                    >> cm.Mab_bulge[j][0]
                    >> cm.Mab_bulge[j][1]
                    >> cm.Mab_bulge[j][2]
                    >> cm.Mab_bulge[j][3]
                    >> cm.Mab_bulge[j][4]
                    >> cm.Mab_bulge[j][5]
                    >> cm.cl_bulge[j]
                    >> cm.typ_bulge[j])) {
              throw std::runtime_error("Unexpected end of bulge.dat");
          }
    }

    CHECK(!(fp2 >> dummy));
    fp2.close();

    // ================================ THICK DISK =============================
    fp2.open("./CMD/components/thick_disk.dat");
    if (!fp2.is_open()) {
        throw std::runtime_error("cannot read ./CMD/components/thick_disk.dat");
    }

    std::getline(fp2, header);   // Skip the header line
    for (size_t j = 0; j < N3; ++j) {
        CHECK(cm.mass_thick[j]   >= 0.0);
        CHECK(cm.logT_thick[j]   >= 0.0);
        CHECK(cm.Mab_thick[j][2] <= 20.0);
//        CHECK(cm.age_thick[j]    <= 8);
        CHECK(cm.age_thick[j]    <= 13);
        CHECK(cm.cl_thick[j]     <= 7);
        CHECK(cm.typ_thick[j]    <= 9.0);
          if (!(fp2 >> cm.mass_thick[j]
                    >> cm.logT_thick[j]
                    >> Mbol
                    >> cm.age_thick[j]
                    >> Pop
                    >> cm.Mab_thick[j][6]
                    >> cm.Mab_thick[j][0]
                    >> cm.Mab_thick[j][1]
                    >> cm.Mab_thick[j][2]
                    >> cm.Mab_thick[j][3]
                    >> cm.Mab_thick[j][4]
                    >> cm.Mab_thick[j][5]
                    >> cm.cl_thick[j]
                    >> cm.typ_thick[j])) {
              throw std::runtime_error("Unexpected end of thick_disk.dat");
          }
    }

    CHECK(!(fp2 >> dummy));
    fp2.close();

    // ================================ STELLAR HALO ===========================
    fp2.open("./CMD/components/halo.dat");
    if (!fp2.is_open()) {
        throw std::runtime_error("cannot read ./CMD/components/halo.dat");
    }

    std::getline(fp2, header);   // Skip the header line
    for (size_t j = 0; j < N4; ++j) {
        CHECK(cm.mass_halo[j]   >= 0.0);
        CHECK(cm.logT_halo[j]   >= 0.0);
        CHECK(cm.Mab_halo[j][2] <= 20.0);
//        CHECK(cm.age_halo[j]    <= 9);
        CHECK(cm.age_halo[j]    <= 14);
        CHECK(cm.cl_halo[j]     <= 7);
        CHECK(cm.typ_halo[j]   <= 9.0);
          if (!(fp2 >> cm.mass_halo[j]
                    >> cm.logT_halo[j]
                    >> Mbol
                    >> cm.age_halo[j]
                    >> Pop
                    >> cm.Mab_halo[j][6]
                    >> cm.Mab_halo[j][0]
                    >> cm.Mab_halo[j][1]
                    >> cm.Mab_halo[j][2]
                    >> cm.Mab_halo[j][3]
                    >> cm.Mab_halo[j][4]
                    >> cm.Mab_halo[j][5]
                    >> cm.cl_halo[j]
                    >> cm.typ_halo[j])) {
              throw std::runtime_error("Unexpected end of halo.dat");
          }
    }

    CHECK(!(fp2 >> dummy));
    fp2.close();

    std::cout << ">>>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
}

// Finds the sightline (one of `nfiles` unique l,b pointings) closest to
// (lon, lat). Call this ONCE per field pointing, not once per star.
int nearestSightline(const extin& ex, double lon, double lat) {
    int best = -1;

    double bestD2 = std::numeric_limits<double>::max();

    for (int k = 0; k < NFILES; k++) {
        double dl = ex.sightlines[k].l - lon;
        double db = ex.sightlines[k].b - lat;

        double d2 = dl * dl + db * db;

        if(d2 < bestD2){
            bestD2 = d2;
            best = k;
        }
    }

    CHECK(best >= 0);

    return best;
}

// Linearly interpolates extinction at distance `dist` (kpc) along sightline k.
// Assumes ex.dist[row..row+nlines) is sorted ascending within the block.
double interpExtinctionAlongSightline(const extin& ex, int k, double dist) {
    const auto& profile = ex.sightlines[k].profile;

    if(dist <= profile.dist[0])
        return profile.ext[0];

    if(dist >= profile.dist[NROWS-1])
        return profile.ext[NROWS-1];

    int lo = 0;
    int hi = NROWS-1;

    while(hi-lo > 1){
        int mid = lo + (hi-lo)/2;

        if(profile.dist[mid] <= dist)
            lo = mid;
        else
            hi = mid;
    }

    double t = (profile.dist[hi] > profile.dist[lo]) ? (dist - profile.dist[lo]) / (profile.dist[hi] - profile.dist[lo]) : 0.0;

    return profile.ext[lo] + t*(profile.ext[hi]-profile.ext[lo]);
}

//////////////////////////////////////////////////
void readBayestar(extin& ex, const std::string& folder) {
    int k = 0;

    for(const auto& entry : fs::directory_iterator(folder)) {
        if(!entry.is_regular_file() || entry.path().extension() != ".txt")
            continue;

        std::ifstream fin(entry.path());

        if(!fin) {
            std::cerr
                << "Cannot open "
                << entry.path()
                << "\n";

            continue;
        }

//        std::cout << "k = " << k << "\tNFILES = " << NFILES << endl;
        CHECK(k < NFILES);
        Sightline& s = ex.sightlines[k];

        fin >> s.l
            >> s.b
            >> s.profile.dist[0]
            >> s.profile.ext[0];
        
        double l, b;
        
        for (int i = 1; i < NROWS; ++i) {
            fin >> l
                >> b
                >> s.profile.dist[i]
                >> s.profile.ext[i];
/*
            if (l != s.l || b != s.b) {
                std::cerr << "Inconsistent l,b in file "
                          << entry.path() << '\n';
                std::exit(EXIT_FAILURE);
            }
*/
        }
    k++;
    }
//    std::cout << "k = " << k << "\tNFILES = " << NFILES << endl;
    CHECK(k == NFILES);

    std::cout
        << "Loaded "
        << k
        << " sightlines\n";
}

double CCM89_a(double lambda_um)
{
    double x = 1.0 / lambda_um;

    if (x < 1.1) {
        return 0.574 * std::pow(x, 1.61);
    }

    if (x < 3.3) {
        double y = x - 1.82;

        return 1.0
             + 0.17699 * y
             - 0.50447 * std::pow(y, 2)
             - 0.02427 * std::pow(y, 3)
             + 0.72085 * std::pow(y, 4)
             + 0.01979 * std::pow(y, 5)
             - 0.77530 * std::pow(y, 6)
             + 0.32999 * std::pow(y, 7);
    }

    throw std::runtime_error("CCM89 wavelength out of supported range");
}

double CCM89_b(double lambda_um)
{
    double x = 1.0 / lambda_um;

    if (x < 1.1) {
        return -0.527 * std::pow(x, 1.61);
    }

    if (x < 3.3) {
        double y = x - 1.82;

        return 1.41338 * y
             + 2.28305 * std::pow(y, 2)
             + 1.07233 * std::pow(y, 3)
             - 5.38434 * std::pow(y, 4)
             - 0.62251 * std::pow(y, 5)
             + 5.30260 * std::pow(y, 6)
             - 2.09002 * std::pow(y, 7);
    }

    throw std::runtime_error("CCM89 wavelength out of supported range");
}

double AlAv(double lambda_um, double Rv)
{
    double x = 1.0 / lambda_um;
    return CCM89_a(x) + CCM89_b(x) / Rv;
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///                 Roman observing-season geometry (Step D1)
///
/// Roman's visit list is a comb: ~70-day observing windows when the bulge is far enough
/// from the Sun, separated by ~110-day gaps when it is not. Everything downstream that
/// asks "did Roman have data near this event's peak?" needs those windows, and the only
/// authoritative statement of them is the epoch times themselves. So they are recovered
/// from the data rather than restated -- see the note above SEASON_GAP_MIN_DAYS in Bulge.h.
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
RomanSchedule buildRomanSchedule(const roman& ro)
{
    // A sorted, de-duplicated copy: the file lists every FIELD at every epoch, so the
    // same instant appears once per field and the raw column is neither unique nor
    // guaranteed monotonic. ~2.4 MB and one sort, once per run.
    std::vector<double> t(ro.tim.begin(), ro.tim.end());
    std::sort(t.begin(), t.end());
    t.erase(std::unique(t.begin(), t.end()), t.end());

    RomanSchedule sch;
    if (t.empty()) return sch;

    sch.missionStart = t.front();
    sch.missionEnd   = t.back();
    sch.minSeasonGap = std::numeric_limits<double>::infinity();

    double start = t.front();
    for (size_t i = 1; i < t.size(); ++i) {
        const double d = t[i] - t[i - 1];
        if (d > SEASON_GAP_MIN_DAYS) {
            sch.seasons.emplace_back(start, t[i - 1]);
            start = t[i];
            // Track the two spacing populations the threshold is meant to separate, so
            // main() can verify the separation is real instead of trusting it.
            sch.minSeasonGap = std::min(sch.minSeasonGap, d);
        } else {
            sch.maxInSeasonSpacing = std::max(sch.maxInSeasonSpacing, d);
        }
    }
    sch.seasons.emplace_back(start, t.back());
    if (sch.seasons.size() == 1) sch.minSeasonGap = 0.0; //nothing was ever classed as a gap

    // Shortest season. A season holding a single epoch has zero length, which is the signature
    // of a schedule sampled more coarsely than the threshold: every epoch becomes its own
    // "season", maxInSeasonSpacing never gets updated at all (so it stays a healthy-looking 0)
    // and minSeasonGap stays above the threshold. Both other margins pass; only this catches it.
    sch.minSeasonLength = std::numeric_limits<double>::infinity();
    for (const auto& sea : sch.seasons)
        sch.minSeasonLength = std::min(sch.minSeasonLength, sea.second - sea.first);

    return sch;
}

double RomanSchedule::dtToSeasonEdge(double t0) const
{
    if (seasons.empty()) return 0.0;

    double best   = std::numeric_limits<double>::infinity();
    bool   inside = false;
    for (const auto& s : seasons) {
        if (t0 >= s.first and t0 <= s.second) inside = true;
        best = std::min(best, std::min(std::fabs(t0 - s.first), std::fabs(t0 - s.second)));
    }
    // Negative inside a season, positive outside. See the header comment on the sign.
    return inside ? -best : best;
}

int RomanSchedule::zone(double t0) const
{
    if (seasons.empty())                      return T0_OFF_MISSION;
    // Off-mission is tested FIRST and wins: a peak before launch or after the mission ends
    // is not a gap in any useful sense, however close it happens to sit to an edge.
    if (t0 < missionStart or t0 > missionEnd) return T0_OFF_MISSION;
    for (const auto& s : seasons)
        if (t0 >= s.first and t0 <= s.second) return T0_IN_SEASON;
    return T0_IN_GAP;
}

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//            Kroupa initial mass function + stellar remnants         //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//

// Draw an INITIAL stellar mass from the Kroupa (2001) broken power law,
//
//     dN/dM  =  k_i * M^-alpha_i
//
// with alpha = 0.3, 1.3, 2.3 on [0.01, 0.08], [0.08, 0.5], [0.5, 120] Msun.
//
// Sampled by inverse CDF, not rejection. Rejection sampling is what the legacy IMnum 2-4
// branches use, and each of them carries hardcoded acceptance bounds (CHECK(f >= 1.4) and
// friends) that are only correct for the 3-5000 Msun range they were written against --
// change the range and the CHECKs fire. Inverse CDF has no such constants: the segment
// weights and the inverse are both computed from the breaks themselves, so the sampler
// stays correct if the breaks are ever retuned.
//
// The k_i are fixed by requiring the IMF be CONTINUOUS at each break, not chosen freely:
//     k1 * 0.08^-0.3 = k2 * 0.08^-1.3  =>  k2 = k1 * 0.08^(1.3-0.3)
//     k2 * 0.50^-1.3 = k3 * 0.50^-2.3  =>  k3 = k2 * 0.50^(2.3-1.3)
// A discontinuous IMF would put a step in the mass distribution at 0.08 and 0.5 Msun, and
// therefore a step in the tE distribution, which is the observable being compared to data.
double drawKroupaInitialMass()
{
    const double lo[3] = {KROUPA_MI_MIN, KROUPA_BREAK1, KROUPA_BREAK2};
    const double hi[3] = {KROUPA_BREAK1, KROUPA_BREAK2, KROUPA_MI_MAX};
    const double al[3] = {KROUPA_ALPHA1, KROUPA_ALPHA2, KROUPA_ALPHA3};

    // Continuity coefficients, k1 fixed to 1 (overall normalisation is irrelevant here --
    // we only ever sample from this, never evaluate an absolute number density).
    double k[3];
    k[0] = 1.0;
    k[1] = k[0] * std::pow(KROUPA_BREAK1, KROUPA_ALPHA2 - KROUPA_ALPHA1);
    k[2] = k[1] * std::pow(KROUPA_BREAK2, KROUPA_ALPHA3 - KROUPA_ALPHA2);

    // Number of stars per segment: k_i * integral of M^-alpha over the segment.
    // None of the Kroupa slopes equals 1, so the (1-alpha) form never divides by zero;
    // the CHECK keeps that assumption honest if the slopes are ever changed.
    double w[3], total = 0.0;
    for (int i = 0; i < 3; ++i) {
        CHECK(std::fabs(1.0 - al[i]) > 1e-9);
        const double p = 1.0 - al[i];
        w[i] = k[i] * (std::pow(hi[i], p) - std::pow(lo[i], p)) / p;
        CHECK(w[i] > 0.0);
        total += w[i];
    }

    // Pick a segment in proportion to its star count, then invert that segment's CDF.
    double u = RandR(0.0, total);
    int seg = 2;
    for (int i = 0; i < 3; ++i) {
        if (u <= w[i]) { seg = i; break; }
        u -= w[i];
    }

    const double p = 1.0 - al[seg];
    const double v = RandR(0.0, 1.0);
    const double lop = std::pow(lo[seg], p);
    const double Mi = std::pow(lop + v * (std::pow(hi[seg], p) - lop), 1.0 / p);

    CHECK(Mi >= KROUPA_MI_MIN * 0.999);
    CHECK(Mi <= KROUPA_MI_MAX * 1.001);
    return Mi;
}


// Map an initial mass to what is still there to act as a lens ~10 Gyr later.
//
// This is where the long-tE tail comes from, and it is not a detail. A star born at 25 Msun
// is long gone, but it left a ~6 Msun black hole -- and since the Einstein radius and the
// event timescale both scale as sqrt(Ml), that black hole lenses for roughly five times as
// long as the 0.3 Msun dwarf next to it. Long events are the ones that span Roman's season
// gaps, so the remnant prescription is what populates the regime this whole project is
// about. Dropping remnants would leave the short-tE yield science intact and quietly
// remove the long-tE precision science.
double remnantMass(double initialMass)
{
    const double Mi = initialMass;

    // Still burning hydrogen: a ~10 Gyr population has a turnoff near 1 Msun, so anything
    // lighter is unevolved and lenses at its birth mass.
    if (Mi < MS_TURNOFF) return Mi;

    // White dwarf. Kalirai et al. (2008) semi-empirical initial-final mass relation,
    // Mf = 0.109*Mi + 0.394, calibrated on open-cluster white dwarfs.
    if (Mi < WD_MI_MAX) return 0.109 * Mi + 0.394;

    // Neutron star. The observed mass distribution is narrow, so a single canonical value
    // is a better model than a spread invented to look sophisticated.
    if (Mi < NS_MI_MAX) return NS_MASS;

    // Black hole. Rough proportional fallback: Ml = 0.24*Mi gives ~4.8 Msun at Mi = 20 and
    // ~28.8 at Mi = 120, spanning the observed stellar-mass black hole range. This is the
    // crudest step in the chain -- black hole remnant masses depend on metallicity and
    // mass loss in ways no single slope captures -- and is flagged in OPEN_ITEMS.md.
    return BH_MI_SLOPE * Mi;
}
