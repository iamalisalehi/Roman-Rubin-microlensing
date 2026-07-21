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
          if(double((l.tE - l.tEs[i-1]) * (l.tE - l.tEs[i])) < 0.0 or l.tE == l.tEs[i-1]) { gg = int(i-1);  break; }
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
          if (double((l.Ml - l.Mls[i-1]) * (l.Ml - l.Mls[i])) < 0.0 or l.Ml == l.Mls[i-1]) { gg = int(i-1);  break; }
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
   double lpi = log10(l.pirel);

   if (lpi <= pi_min)      gg = 0;
   else if (lpi >= pi_max) gg = GG;
   else {
      for (int i = 1; i <= GG; ++i) {
          if (double((lpi - l.pis[i-1]) * (lpi - l.pis[i])) < 0.0 or lpi == l.pis[i-1]) { gg = int(i-1);  break; }
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
          if (double((l.u0 - l.u0s[i-1]) * (l.u0 - l.u0s[i])) < 0.0 or l.u0 == l.u0s[i-1]) { gg = int(i-1);  break; }
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
          if (double((mur - l.mus[i-1]) * (mur - l.mus[i])) < 0.0 or mur == l.mus[i-1]) { gg = int(i-1);  break; }
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
          if (double((gadr - l.mbs[i-1]) * (gadr - l.mbs[i])) < 0.0 or gadr == l.mbs[i-1]) { gg = int(i-1); break; }
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
          if(double((blen - l.fbs[i-1]) * (blen - l.fbs[i])) < 0.0 or blen == l.fbs[i-1]) { gg = int(i-1); break;}
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
    x = pow(10.0, float(0.4 * (mag- sig5 )));
    Delta1 = double(sqrt(fabs((0.04 - gama[fi]) * x + gama[fi] * x * x)));
    
    if (Delta1 < 0.0001)   Delta1 = 0.0001;

    CHECK(sig5 >= 0.0);
    CHECK(sig5 <= 40.0);
    CHECK(Delta1 >= 0.00001);
    CHECK(x > 0.0);
   
    return sqrt(delta2 * delta2 + Delta1 * Delta1);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double errlsstA(lsst & ls, double ghadr){ //LSST Astrometric Error 

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
    int j = 0;
    double Mbol, Pop;
//    double logL, gravity, metal, B, V, R, I, J, H;

    std::ifstream fp2("./BC/components/thin_disk.dat");
    if (!fp2.is_open()) {
        throw std::runtime_error("cannot read ./BC/components/thin_disk.dat");
    }

    while (fp2 >> cm.mass_thin[j] >> cm.logT_thin[j]    >> Mbol               >> cm.age_thin[j]     >> Pop
            >> cm.Mab_thin[j][6]  >> cm.Mab_thin[j][0]  >> cm.Mab_thin[j][1]  >> cm.Mab_thin[j][2]  >> cm.Mab_thin[j][3]  
            >> cm.Mab_thin[j][4]  >> cm.Mab_thin[j][5]  >> cm.cl_thin[j]      >> cm.typ_thin[j]) {

        CHECK(cm.mass_thin[j]   >= 0.0);
        CHECK(cm.logT_thin[j]   >= 0.0);
        CHECK(cm.Mab_thin[j][2] <= 20.0);
        CHECK(cm.age_thin[j]    <= 10);
        CHECK(cm.cl_thin[j]     <= 7);
        CHECK(cm.typ_thin[j]    <= 9.0);
       
        j++;
    }

    fp2.close();
    CHECK(j == N1);

    // ================================ BULGE ==================================
    j = 0;
    fp2.open("./BC/components/bulge.dat");
    if (!fp2.is_open()) {
        throw std::runtime_error("cannot read ./BC/components/bulge.dat");
    }

    while (fp2 >> cm.mass_bulge[j]  >> cm.logT_bulge[j]   >> Mbol                >> cm.age_bulge[j]     >> Pop
            >> cm.Mab_bulge[j][6]  >> cm.Mab_bulge[j][0]  >> cm.Mab_bulge[j][1]  >> cm.Mab_bulge[j][2]  >> cm.Mab_bulge[j][3]  
            >> cm.Mab_bulge[j][4]  >> cm.Mab_bulge[j][5]  >> cm.cl_bulge[j]      >> cm.typ_bulge[j]) {

        CHECK(cm.mass_bulge[j]   >= 0.0);
        CHECK(cm.logT_bulge[j]   >= 0.0);
        CHECK(cm.Mab_bulge[j][2] <= 18.0);
        CHECK(cm.age_bulge[j]    <= 10);
        CHECK(cm.cl_bulge[j]     <= 7);
        CHECK(cm.typ_bulge[j]    <= 9.0);
      
        j++;
    }

    fp2.close();
    CHECK(j == N2);

    // ================================ THICK DISK =============================
    j = 0;
    fp2.open("./BC/components/thick_disk.dat");
    if (!fp2.is_open()) {
        throw std::runtime_error("cannot read ./BC/components/thick_disk.dat");
    }
    while (fp2 >> cm.mass_thick[j] >> cm.logT_thick[j]    >> Mbol                >> cm.age_thick[j]     >> Pop
            >> cm.Mab_thick[j][6]  >> cm.Mab_thick[j][0]  >> cm.Mab_thick[j][1]  >> cm.Mab_thick[j][2]  >> cm.Mab_thick[j][3]  
            >> cm.Mab_thick[j][4]  >> cm.Mab_thick[j][5]  >> cm.cl_thick[j]      >> cm.typ_thick[j]) {

        CHECK(cm.mass_thick[j]   >= 0.0);
        CHECK(cm.logT_thick[j]   >= 0.0);
        CHECK(cm.Mab_thick[j][2] <= 20.0);
        CHECK(cm.age_thick[j]    <= 8);
        CHECK(cm.cl_thick[j]     <= 7);
        CHECK(cm.typ_thick[j]    <= 9.0);

        j++;
    }

    fp2.close();
    CHECK(j == N3);

    // ================================ STELLAR HALO ===========================
    j = 0;
    fp2.open("./BC/components/halo.dat");
    if (!fp2.is_open()) {
        throw std::runtime_error("cannot read ./BC/components/halo.dat");
    }
    while (fp2 >> cm.mass_halo[j] >> cm.logT_halo[j]    >> Mbol               >> cm.age_halo[j]     >> Pop
            >> cm.Mab_halo[j][6]  >> cm.Mab_halo[j][0]  >> cm.Mab_halo[j][1]  >> cm.Mab_halo[j][2]  >> cm.Mab_halo[j][3]  
            >> cm.Mab_halo[j][4]  >> cm.Mab_halo[j][5]  >> cm.cl_halo[j]      >> cm.typ_halo[j]) {

        CHECK(cm.mass_halo[j]   >= 0.0);
        CHECK(cm.logT_halo[j]   >= 0.0);
        CHECK(cm.Mab_halo[j][2] <= 20.0);
        CHECK(cm.age_halo[j]    <= 9);
        CHECK(cm.cl_halo[j]     <= 7);
        CHECK(cm.typ_halo[j]   <= 9.0);

        j++;
    }

    fp2.close();
    CHECK(j == N4);

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
        
            if (l != s.l || b != s.b)
            {
                std::cerr << "Inconsistent l,b in file "
                          << entry.path() << '\n';
                std::exit(EXIT_FAILURE);
            }
        }

    }
    CHECK(k == NFILES);

    std::cout
        << "Loaded "
        << k
        << " sightlines\n";
}
