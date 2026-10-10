// The core structs shared by most modules: one source star, one lens, the observer/Earth state,
// and one row of the per-event output table (EventRecord).
#ifndef ROMAN_TYPES_H
#define ROMAN_TYPES_H

#include "common.h"

///============================================================================
static std::vector<double> make_grid(double min, double max)
{
    std::vector<double> v(GG + 1);

    for (int i = 0; i <= GG; ++i)
    {
        v[i] = min + (max - min) * static_cast<double>(i) / GG;
    }

    return v;
}

// Geometric spacing, for a quantity whose range spans decades (a linear grid over 3-1000 Msun
// would put most of a log-uniform black-hole population in the first bin). Used only where the
// population asks for it.
static std::vector<double> make_grid_log(double min, double max)
{
    std::vector<double> v(GG + 1);
    const double lo = std::log(min), hi = std::log(max);

    for (int i = 0; i <= GG; ++i)
    {
        v[i] = std::exp(lo + (hi - lo) * static_cast<double>(i) / GG);
    }

    return v;
}

///============================================================================
struct source {
    int nums, cl;

    double mass, logT, typ, age;
    // Finite-source size: the physical radius [Rsun], the angular radius theta* = Rstar / Ds [mas, as
    // l.tetE], and rho = theta* / thetaE [], the source radius in Einstein radii (func_lens).
    double Rstar, thetaStar, rho;
    double mus1, mus2, mus;
    double xv, yv, zv;
    double Av; //, Avv;
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
    double pos1b, pos2b, pos1c, pos2c, def1a, def2a, def1c, def2c;
    double errM, errA, FWHM;
    double ut, ut0, Astar, xi, ux, uy;
    double Ds, TET, FI, lat, lon, vs;
    double Nstart, Rostart, Romaxs, Romins, nstart, nstarti;
    double od_thin, od_thick, od_bulge, od_halo, opt;

    std::array<double, 2> fb, mbs; // small fixed arrays
    // Share of each telescope's baseline flux that is the LENS's own light: 0 for dark lenses; for a
    // luminous lens it is already inside the blend (fb counts the source only).
    std::array<double, 2> fLens{};
    // Light centroid of the unresolved neighbours relative to the source's unlensed position, per
    // telescope [tt][x, y], in mas (func_source). Their light pulls the measured centroid toward this point.
    std::array<std::array<double, 2>, 2> blendOff{};

    std::vector<double> nssim; // size Num
    std::vector<double> nsdet; // size Num
    std::vector<double> rho_thin;
    std::vector<double> rho_thick;
    std::vector<double> rho_halo;
    std::vector<double> rho_stars;
    std::vector<double> rho_bulge;
    std::vector<double> Rostar0;
    std::vector<double> Rostari;
    std::vector<double> Nstari;

    std::vector<double> nsbl;  // size M+1
    std::vector<double> blend; // size M+1
    std::vector<double> Fluxb; // size M+1
    std::vector<double> magb;  // size M+1
    std::vector<double> Ai;    // size M+1
    std::vector<double> Mab;   // size M+1
    std::vector<double> Map;   // size M+1

    GalacticComponent struc;
    // Constructor
    source()
        : nssim(Num), nsdet(Num),
          rho_thin(Num), rho_thick(Num), rho_halo(Num), rho_stars(Num), rho_bulge(Num),
          Rostar0(Num), Rostari(Num), Nstari(Num),
          nsbl(M), blend(M), Fluxb(M), magb(M), Ai(M), Mab(M), Map(M)
    {}
};

struct lens {
    int numl;
    bool luminous = false;          // a living (main-sequence) star, whose light is blended

    double Ml, Dl, vl, Vt, xls, u0, A0, mi1, mi2;
    double rhomaxl, tE, RE, t0, murel, DeltaT;
    double piE, pirel, tetE, pos1, pos2, mul1, mul2, mul;
    double betal, betas, deltal, deltas, deltao;

    std::array<double, 2> Nhalo, Nself; // small fixed arrays

    std::vector<double> nstE; // size GG+1
    std::vector<double> ndtE; // size GG+1
    std::vector<double> NstE; // size GG+1
    std::vector<double> NdtE; // size GG+1
    std::vector<double> NsMl; // size GG+1
    std::vector<double> NdMl; // size GG+1
    std::vector<double> Nspi; // size GG+1
    std::vector<double> Ndpi; // size GG+1
    std::vector<double> Nsu0; // size GG+1
    std::vector<double> Ndu0; // size GG+1
    std::vector<double> Nsmb; // size GG+1
    std::vector<double> Ndmb; // size GG+1
    std::vector<double> Nsfb; // size GG+1
    std::vector<double> Ndfb; // size GG+1
    std::vector<double> Nsmu; // size GG+1
    std::vector<double> Ndmu; // size GG+1
    std::vector<double> timn; // size coun
    std::vector<double> magn; // size coun
    std::vector<double> soux; // size coun
    std::vector<double> souy; // size coun
    std::vector<double> errm; // size coun
    std::vector<double> erra; // size coun
    std::vector<int> tele;    // size coun
    std::vector<int> rseas;   // size coun: Roman season index of the epoch (-1 for Rubin)
    std::vector<int> rroll;   // size coun: Roman roll of the epoch, 0 spring / 1 autumn (-1 Rubin)
    
    std::vector<double> tEs;  // size GG+1
    std::vector<double> Mls;  // size GG+1
    std::vector<double> pis;  // size GG+1
    std::vector<double> u0s;  // size GG+1
    std::vector<double> mbs;  // size GG+1
    std::vector<double> fbs;  // size GG+1
    std::vector<double> mus;  // size GG+1

    GalacticComponent struc;
    // Constructor
    lens()
        : nstE(GG+1), ndtE(GG+1), NstE(GG+1), NdtE(GG+1),
          NsMl(GG+1), NdMl(GG+1),  Nspi(GG+1), Ndpi(GG+1),
          Nsu0(GG+1), Ndu0(GG+1),
          Nsmb(GG+1), Ndmb(GG+1),
          Nsfb(GG+1), Ndfb(GG+1),
          Nsmu(GG+1), Ndmu(GG+1),
          timn(coun), magn(coun), soux(coun), souy(coun), errm(coun), erra(coun),
          tele(coun), rseas(coun, -1), rroll(coun, -1),

          tEs(make_grid(tE_min, tE_max)),
          Mls(gPop->logGrid ? make_grid_log(mlMin(), mlMax()) : make_grid(mlMin(), mlMax())),
          pis(make_grid(pi_min, pi_max)),
          u0s(make_grid(u0_min, u0_max)),
          mbs(make_grid(mb_min, mb_max)),
          fbs(make_grid(fb_min, fb_max)),
          mus(make_grid(mu_min, mu_max))
    {}
};

struct astromet{
   double Ve_n1, Ve_n2;
   double ue_n1, ue_n2;

   // Multiplies L2_OFFSET_AU in lightcurve(): 1 = Roman at L2 (physical), 0 = Roman at the centre of
   // the Earth. A command-line switch for the satellite-parallax comparison.
   double satScale = 1.0;
};

struct EventRecord {
    int    counter, flagL;
    double tE, RE, piE, tetE, Vt, u0, Ml;
    double opt, Dl, Ds, vl, vs;
    double mbase, fblend;
    int    gg, struc;
    double FWHM, vsave, DeltaT, murel;
    double resu0, resu1, resu2, resu3, resu5, resu9, resu10, resu13, resu14;
    double Map2, nsbl2;
    int    flagi;
    double Ai2;

    // ---- per-survey bookkeeping ----
    // Lets downstream code tell a Rubin-only detection from a Roman-only one and apply the "both
    // surveys have data" cut. The positional aggregate initialiser depends on member order.
    int    ndwL, ndwR;          //epochs actually contributed by each survey
    int    detL, detR, detJ;    //independent per-survey and joint detection outcomes
    int    okJoint, okRubin, okRoman;   //1 = that Fisher partition is characterizable
    double sigtE_J,   sigtE_L,   sigtE_R;    //absolute 1-sigma on tE   [days]
    double sigpiE_J,  sigpiE_L,  sigpiE_R;   //absolute 1-sigma on piE  []
    double sigtetE_J, sigtetE_L, sigtetE_R;  //absolute 1-sigma on tetE [mas]
    int    detCls;              //DetClass -- which combination of surveys detected it
    int    synClass;            //SynergyClass -- see the note above the enum. Do NOT drop rows
                                //with a not-characterizable single-survey sigma; use this.
    double condA_J, condA_L, condA_R;   //photometric condition numbers (normalized matrix);
                                //-1 where that partition was not characterizable. Stored so the
                                //ill-conditioning cut can be chosen in analysis, not baked in.

    // ---- the rest of the per-event row ----
    // Order matters: the positional aggregate initialiser at the push_back site and every
    // downstream column index depend on it.
    double t0;                  //time of closest approach [d] on the simulation clock
                                //(day 0 = 2026-04-11, the first Rubin bulge visit)
    double xi;                  //source-trajectory angle [rad]; sets how the annual
                                //parallax ellipse projects onto the source track
    double lon, lat;            //Galactic sightline [deg]. Present so "per field" is a cut
                                //on this table rather than separate files.
    double mbs1, fb1;           //Roman F146 baseline magnitude [mag] and source-flux
                                //fraction []. mbs0/fb0 above are Rubin's. With RUBIN_BANDS = {r}
                                //these duplicate magb[6]/blend[6] and magb[2]/blend[2]; they are
                                //the quantities the Fisher matrix fits (params 6-8).
    std::array<double, M> magb; //blended baseline magnitude per FILTER [mag]
                                //(0-5 = LSST ugrizy, 6 = Roman F146) -- all stars in the
                                //seeing disc, not just the source
    std::array<double, M> blend;//fraction of aperture flux from the SOURCE, per filter []
    double relMl_J, relMl_L, relMl_R;  //fractional sigma(Ml); -1 = not measurable
    int    okB_J, okB_L, okB_R; //astrometric matrix inverted? NOT implied by okA: a row can
                                //have okRubin=1 while sigtetE_L is -1.
    double condB_J, condB_L, condB_R;  //normalized astrometric condition numbers; -1 as above
    double dtEdge;              //signed days from t0 to the nearest Roman season edge;
                                //NEGATIVE means t0 fell inside a season. See RomanSchedule.
    int    t0zone;              //T0Zone: 0 in-season, 1 mid-mission gap, 2 off-mission

    // ---- resolving the two images, counted per recorded epoch ----
    // Number of that survey's epochs at which both images were within the filter's
    // [saturation, single-visit depth] AND separated by more than the stated bar. Counts, since the
    // criterion is "at least three such epochs". nres5/nres20 use D*sigma_a at two anchors;
    // nresPSF uses the filter's PSF FWHM.
    long   nres5_L, nres20_L, nresPSF_L;
    long   nres5_R, nres20_R, nresPSF_R;
    // Largest separation reached at an epoch where both images were detectable [mas].
    // -1 means that never happened (sentinel); it must not enter a mean or a histogram.
    double dsepMax_L, dsepMax_R;
};

#endif // ROMAN_TYPES_H
