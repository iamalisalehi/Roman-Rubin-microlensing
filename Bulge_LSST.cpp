#include "Bulge.h"

///============================================================================

time_t _timeNow;
unsigned int _randVal;
unsigned int _dummyVal;
FILE * _randStream;

///==============================================================//
///     Instrument-agnostic epoch matching (LSST or Roman)       //
///                                                              //
///==============================================================//
// Extracted from the sky-position loop in main() so it can be called once per
// instrument instead of being hardwired to `ls`. Behavior is unchanged for LSST;
// calling it a second time with Roman's own l/b/tim arrays and its own FoV is what
// gives Roman its own epoch list instead of inheriting LSST's cadence.
//
// `label` is purely diagnostic ("LSST" / "Roman") — it identifies which call
// printed a given warning, since both calls share this one function.
//
// Two visits can legitimately share an identical recorded timestamp for a given
// sky position — most commonly two different fields/pointings whose FoV circles
// overlap and which happen to share the same observing schedule (this is exactly
// what can happen between adjacent Roman fields: FoVRoman=0.28 vs ~0.41 field
// spacing means neighboring fields' circles overlap). There is no meaningful
// cadence between two simultaneous visits, so rather than treat this as fatally
// corrupt data, we keep the first and skip the duplicate — but we log it, because
// if this fires constantly (not just occasionally near field boundaries) that's a
// sign of a genuine sorting/data problem in the baseline file that needs fixing at
// the source, not papering over here.
int matchVisibleEpochs(const char* label, double lon, double lat, double fov,
                       const std::vector<double>& l_arr,
                       const std::vector<double>& b_arr,
                       const std::vector<double>& tim_arr,
                       int nEpochs,
                       std::vector<int>& ct,
                       double& minCadence) {
    int ndd = 0;
    int nTies = 0;
    minCadence = 100000.0;

    for (int i = 0; i < 1000; ++i) ct[i] = -1; // same reset as before

    for (int i = 0; i < nEpochs; ++i) {
        double dl = lon - l_arr[i];
        double db = lat - b_arr[i];
        if (std::sqrt(dl * dl + db * db) <= fov) {

            if (ndd > 0) {
                double cade = tim_arr[i] - tim_arr[ct[ndd - 1]];

                if (cade <= 0.0) {
                    nTies++;
                    if (nTies <= 5) {
                        std::cerr << "[matchVisibleEpochs:" << label << "] skipping tied/out-of-order "
                                  << "epoch at index " << i << " (tim=" << tim_arr[i]
                                  << "), previous kept epoch tim=" << tim_arr[ct[ndd - 1]] << "\n";
                    }
                    continue; // keep the earlier one, don't record this as a new epoch
                }

                if (minCadence > cade) minCadence = cade;
                if (ndd >= 999) break;
                CHECK(tim_arr[i] >= 0.0);
                CHECK(tim_arr[i] <= Tobs);
                CHECK(minCadence > 0.0);
            }
            ct[ndd] = i;
            ndd += 1;
        }
    }

    if (nTies > 5) {
        std::cerr << "[matchVisibleEpochs:" << label << "] ... " << (nTies - 5)
                  << " more skipped (total " << nTies << " tied/out-of-order epochs for this sky position)\n";
    } else if (nTies > 0) {
        std::cerr << "[matchVisibleEpochs:" << label << "] " << nTies << " tied/out-of-order epoch(s) skipped\n";
    }

    return ndd;
}
// TODO(Ali): PLACEHOLDER. errlsstM interpolates against a per-visit `sig5` (5-sigma
// depth), which varies visit-to-visit for a ground-based survey (airmass, sky
// brightness, seeing). Roman is space-based with far more uniform per-visit depth,
// so a simple mag-vs-error lookup (no per-visit depth term) may be all you need —
// but confirm that's actually what sigma_roman.txt encodes before trusting this.
// This does nearest-neighbor lookup in ro.mag; swap in linear interpolation if
// sigma_roman.txt's mag sampling is coarse.
double errRomanM(const roman& ro, double mag)
{
    int    best = 0;
    double bestDiff = std::fabs(ro.mag[0] - mag);
    for (int i = 1; i < NaRoman; ++i) {
        double diff = std::fabs(ro.mag[i] - mag);
        if (diff < bestDiff) { bestDiff = diff; best = i; }
    }
    return ro.err[best];
}

///==============================================================//
///                                                              //                                                    /
///                  Main program                                //
///                                                              //
///==============================================================//

#ifndef FISHER_FIXTURE_BUILD
// main() is excluded when this translation unit is linked into the Fisher-matrix
// test fixture (tests/fisher_fixture.cpp), which supplies its own main(). Everything
// below main() -- FisherM, ErrorCal, lightcurve, invert_matrix -- is what the fixture
// links against, so it must stay outside this guard.
int main() {
    srand(static_cast<unsigned int>(time(0)));

    // --------------------- Allocate objects ------------------------
    auto s  = std::make_unique<source>();
    auto l  = std::make_unique<lens>();
    auto as = std::make_unique<astromet>();
    auto cm = std::make_unique<CMD>();
//    auto ga = std::make_unique<galactic>();
    auto ex = std::make_unique<extin>();
    auto ls = std::make_unique<lsst>();
    auto ro = std::make_unique<roman>();
    auto co = std::make_unique<covarian>();
    
    std::vector<EventRecord> records;
    records.reserve(1000); // rough upper bound on icon per field

    // --------------------- Read BulgeBaseline.dat ------------------
    std::ifstream fil("./Baseline/BulgeBaseline.dat");
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
        fil >> ID >> ls->RA[i] >> ls->DEC[i] >> ls->l[i] >> ls->b[i]
            >> ls->tim[i] >> ls->filter[i] >> airm >> seeingVal >> skyB
            >> TV >> ls->sig5[i] >> texp >> ls->dist[i];
//if (i == 0) cout << ID << endl;
        CHECK(airm >= 0.0);
        CHECK(ls->filter[i] >= 0);
        CHECK(ls->filter[i] < 6);
        CHECK(ls->l[i] >= l1 - FoV);
        CHECK(ls->l[i] <= l2 + FoV);
        CHECK(ls->b[i] >= b1 - FoV);
        CHECK(ls->b[i] <= b2 + FoV);
        CHECK(ls->tim[i] >= 0.0);
        CHECK(ls->tim[i] <= Tobs);
        CHECK(texp >= 0.0);
   }
    fil.close();
    std::cout << "**** File BulgeBaseline.dat was read ****\n";

    // --------------------- Read sigmaA_LSST.txt -------------------
    // Atrometric error?
    fil.open("./files/sigmaA_LSST.txt");
    if (!fil) { std::cerr << "Cannot read sigmaA_LSST.txt\n"; return 1; }

    for (int i = 0; i < Na; ++i) {
        fil >> ls->mag[i] >> ls->err[i];

        CHECK(ls->mag[i] >= 16.0);
        CHECK(ls->mag[i] <= 25.0);
        CHECK(ls->err[i] >= 0.2);
        CHECK(ls->err[i] <= 5.0);
    }
    fil.close();
    std::cout << "**** File sigmaA_LSST.txt was read ****\n";

    // --------------------- Read sigma_Roman.txt ---------------------
    // Photometric error
    fil.open("./files/sigma_roman.txt");
    if (!fil) { std::cerr << "Cannot read sigma_roman.txt\n"; return 1; }

    for (int i = 0; i < NaRoman; ++i) {
        fil >> ro->mag[i] >> ro->err[i];
        CHECK(ro->mag[i] > 12.0);
        CHECK(ro->err[i] > 0.0);
    }
    fil.close();
    std::cout << "**** File sigma_roman.txt was read ****\n";
 
    // --------------------- Read RomanBaseline.dat -------------------
    // TODO(Ali): generate this file from the ROTAC 2025 overguide season/cadence design
    // (analogous to readbaselineBulge.py, but sourced from Roman's own season structure
    // rather than an LSST OpSim). Format assumed to mirror BulgeBaseline.dat's columns
    // that matter for matching: ID RA DEC l b time [sig5] ... — adjust the >> list below
    // to whatever columns you actually emit.
    fil.open("./Baseline/RomanBaseline.dat");
    if (!fil) {
        std::cerr << "Cannot read RomanBaseline.dat\n";
        return 1;
    }
    std::getline(fil, header); // skip header line
    for (int i = 0; i < NlRoman; ++i) {
        fil >> ID >> ro->RA[i] >> ro->DEC[i] >> ro->l[i] >> ro->b[i] >> ro->tim[i] >> ro->sig5[i];

        CHECK(ro->tim[i] >= 0.0);
        CHECK(ro->tim[i] <= Tobs);
        // Deliberately no filter CHECK here — Roman is single-band (F146, index 6) for now.
    }
    fil.close();
    std::cout << "**** File RomanBaseline.dat was read ****\n";

    // --------------------- Read extinction ------------------------
    readBayestar(*ex,"./files/ext/");
    std::cout << "**** File extinctionf.txt was read ****\n";

    // --------------------- Call read_cmd --------------------------
    read_cmd(*cm);
    std::cout << "******* read_cmd was done ************" << std::endl;

    int    save = 0, flagm, flagL, gg = -1; //ss, qq, ww, vv, zz, pp, counter,
    int    nri = -1, nde = -1, icon;
    int    nlens;// hh; // nde1, nri1,
    int    nDetNeither, nDetRubinOnly, nDetRomanOnly, nDetJointOnly, nDetBoth; // per-field label counts
    int    gi,       ndw, sq, ndd;
    int    giR,      sqR, nddR; // Roman-side cursor/count, parallel to gi/sq/ndd
    int    flag_det; // nml = 0;
    int    ndw_L, ndw_R;           // per-instrument epoch counts (ndw stays the joint/shared total)
    int    flag_det_L, flag_det_R; // per-instrument run-test result (flag_det stays the joint one)
    int    detL, detR, detJ;       // per-instrument / joint detection booleans; FFG[0] = detL or detR or detJ
    int    flagf,  fi; // datf1, datf2;
    double errs,   errg, minc, cade; // fel, , mind
    double fdetRubin, testL, testR; // Step B2: per-survey pre-selection (fdet retired)
    bool   rubinDetectable, romanDetectable, acceptRubin, acceptRoman;
    double mincR, cadeR, errgR; // Roman-side cadence/error tracking, parallel to minc/cade/errg
    // PLACEHOLDER: no Roman astrometric error model exists yet — see the TODO(Ali) on
    // errsR's use in the Roman branch below, and JOINT_FIT_REFACTOR_PLAN.md's Deferred
    // section ("The Roman astrometric error model decision"). errsR/magnioR are scratch
    // for the placeholder computation, standing in until a real errRomanA()-style
    // function + F146 astrometric-error dataset exist.
    double errsR, magnioR; // Roman-side astrometric-error placeholder / noisy-magnitude scratch
    double magnio, test, deltaA; // dist,
    double Astar0, As1,  As0;
    double initial;
    double trajm, trajp; //  ddf;
    double chi1,  chi2,   chi3, chi1a, chi2a, chi3a, sil, sil2;
    // Per-instrument duplicates. chi1/chi2/chi3/chi1a/chi2a/chi3a above are the JOINT
    // accumulators (fed by both branches); _L/_R are Rubin-only/Roman-only respectively.
    double chi1_L, chi2_L, chi3_L, chi1a_L, chi2a_L, chi3a_L;
    double chi1_R, chi2_R, chi3_R, chi1a_R, chi2a_R, chi3a_R, silR, sil2R;
    double dchiL, dchiP,  dchiA;
    double dchiL_L, dchiP_L, dchiA_L, dchiL_R, dchiP_R, dchiA_R;
    double flag0, flag1,  flag2;
    double flag0_L, flag1_L, flag2_L, flag0_R, flag1_R, flag2_R;
    double vs1,   vs2,    def1p,  def2p, vsave, dt,    Mpeak;
    double ErtE,  ErpiE,  ErtetE, Erml,  Erdl,  Ermul, Ermus, Eru0, Erfb, nsim;
    double mbase, fblend, Gamma,  Neven, EFF,   EffiD, EffiL, nerr=0.0;
//    double shib,  Efi;
   
    std::array<int, 3> FFG;
    std::array<double, M> magni, magni0;
    std::array<double, 2> tE, RE, piE, tetE, Vt, u0, Ml, opd, Dl, Ds, vl, vs;
    std::array<double, 2> numd, Struc, murel, vsn, DelT, mbs, fb, fwhm, Map, nbl, Ext;
   
    //fil=fopen("./files/MONTLMC/files/EfLMC2B.dat","r");
    //for(int i=0; i<=GG; ++i){
    //fscanf(fil,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",//18
    //&l.NdtE[i],&l.NstE[i], &l.NdMl[i],&l.NsMl[i], &l.Ndpi[i],&l.Nspi[i],   &l.Ndu0[i], &l.Nsu0[i], &l.Ndmb[i],&l.Nsmb[i], 
    //&l.Ndfb[i],&l.Nsfb[i], &l.Ndmu[i],&l.Nsmu[i], &l.Nhalo[1],&l.Nhalo[0], &l.Nself[1],&l.Nself[0]);}
    //fclose(fil); 
    //cout<<"**** File extinctionf.txt was read ****"<<endl;    
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
// Create/clear files if IMnum == 1
    if (IMnum == 1) {
        std::string filnam0 = "./files/MONTLMC/files/BHLSSTMONTS.dat";
        std::ofstream(filnam0).close(); // create/clear file
    }

    // File names
//    std::string filnam0 = "./files/MONTLMC/files/BHLSSTMONTS.dat"; // now visible outside the if
    std::string filnam1 = "./files/MONTLMC/files/magC"   + std::to_string(save)  +  ".dat";
    std::string filnam2 = "./files/MONTLMC/files/datC"   + std::to_string(save)  +  ".dat";
    std::string fnLDt   = "./files/MONTLMC/files/LpLMC"  + std::to_string(IMnum) +  ".dat";
    std::string fnEff   = "./files/MONTLMC/files/EfLMC"  + std::to_string(IMnum) +  ".dat";
    std::string fnEffB  = "./files/MONTLMC/files/EfLMC"  + std::to_string(IMnum) + "B.dat";
    std::string fnGam   = "./files/MONTLMC/files/MapLMC" + std::to_string(IMnum) +  ".dat";
    std::string testf   = "./test"                       + std::to_string(IMnum) +  ".dat";

    // Open files
    std::ifstream fil0(fnLDt);
//    std::ifstream fil1(filnam0);
    std::ofstream fil2(fnEff);
    std::ofstream fil2b(fnEffB);
    std::ofstream fil4(filnam1);
    std::ofstream fil5(filnam2);
    std::ofstream fil3(fnGam, std::ios::app);
    std::ofstream filg_in(testf);

    // Check all
    if (!fil0 || !fil2 || !fil2b || !fil3 || !fil4 || !fil5 || !filg_in) {
        std::cerr << "Cannot open one or more files!" << std::endl;
        return 1;
    }

    save = 0;

///HHHHHHHHHHHHHHHHHHHHH Monte Carlo Simulation HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

//    for (s->lon = l1 - wid; s->lon <= l2 + wid; s->lon += dd) {
    for (s->lon = 0.5; s->lon <= 0.6; s->lon += dd) {
//    for (s->lon = 1.0; s->lon <= 1.1; s->lon += dd) {
        nri +=  1;
        nde  = -1;
//        for (s->lat = b1 - wid; s->lat <= b2 + wid; s->lat += dd) {
        for (s->lat = -1.0; s->lat <= -0.9; s->lat += dd) {
//        for (s->lat = -3.9; s->lat <= -3.8; s->lat += dd) {
            if (s->lon < lx and s->lat > bx) {
                continue;
            };
            nde += 1;
            cout << ">>>>>>>>>>> NEW STEP " << nde << " <<<<<<<<\t nri:  " << nri << endl;
            cout << "longtitude: " << s->lon << "\t latitude: " << s->lat << endl;

            cade = 0.0;

            ndd  = matchVisibleEpochs("LSST", s->lon, s->lat, FoV, ls->l, ls->b, ls->tim, Nl, ls->ct, minc);
            nddR = matchVisibleEpochs("Roman", s->lon, s->lat, FoVRoman, ro->l, ro->b, ro->tim, NlRoman, ro->ct, mincR);

            cout << "ndd (LSST): "  << ndd  << "\t minc (LSST): "  << minc  << endl;
            cout << "ndd (Roman): " << nddR << "\t minc (Roman): " << mincR << endl;
 
            icon  = 0;
            nlens = 0;
            nDetNeither = 0; nDetRubinOnly = 0; nDetRomanOnly = 0; nDetJointOnly = 0; nDetBoth = 0;
            nsim  = 0.0;
            nerr  = 0.0;
            for (int i = 0; i < Num; ++i) { s->nssim[i] = 0.0;  s->nsdet[i] = 0.0; }
            for (int i = 0; i <= GG; ++i) { l->nstE[i]  = 0.0;  l->ndtE[i]  = 0.0; }

            s->TET = (360.0 - s->lon) / RAa;///radian s.lon/RA;//
            s->FI  = s->lat / RAa;
            
            Disk_model(*s, 1);
            int sightlineIdx = nearestSightline(*ex, s->lon, s->lat);

            records.clear();
            do { //Start of visible star
                nsim += 1.0;
                func_source(*s, *cm, *ex, sightlineIdx);
                func_lens(*l, *s);
//                std::cerr << "nsim=" << nsim << "  Ds=" << s->Ds << "  mass=" << s->mass
//                          << "  nums=" << s->nums << "  Ml=" << l->Ml << "  u0=" << l->u0 << "\n";
                optical_depth(*s);

                s->nssim[s->nums] += 1.0;
                flagf   = 0;
                flagm   = 0;
                initial = 0.0;
                // (Step B2: the old single `test = RandR(0.0,1.0)` draw consumed here by
                // `test <= s->blend[2]` is gone — testL/testR are now drawn fresh right
                // before the per-survey pre-selection check, below.)

                for (int i = 0; i < coun; ++i) {
                    l->timn[i] = 0.0;  l->magn[i] = 0.0; l->soux[i] = 0.0;  l->souy[i] = 0.0;
                    l->errm[i] = 0.0;  l->erra[i] = 0.0; l->tele[i] = -1;
                }

                ndw     = 0;   flag_det = 0;
                ndw_L   = 0;   ndw_R    = 0;
                flag_det_L = 0; flag_det_R = 0;
                flag0   = 0.0; flag1    = 0.0; flag2 = 0.0;
                flag0_L = 0.0; flag1_L  = 0.0; flag2_L = 0.0;
                flag0_R = 0.0; flag1_R  = 0.0; flag2_R = 0.0;
                chi1    = 0.0; chi2     = 0.0; chi3  = 0.0;
                chi1_L  = 0.0; chi2_L   = 0.0; chi3_L = 0.0;
                chi1_R  = 0.0; chi2_R   = 0.0; chi3_R = 0.0;
                chi1a   = 0.0; chi2a    = 0.0; chi3a = 0.0;
                chi1a_L = 0.0; chi2a_L  = 0.0; chi3a_L = 0.0;
                chi1a_R = 0.0; chi2a_R  = 0.0; chi3a_R = 0.0;
                def1p   = 0.0; s->def1c = 0.0; vsave = 0.0;
                def2p   = 0.0; s->def2c = 0.0;
                s->errM = 0.0; s->errA  = 0.0;
                dchiL   = 0.0; dchiP    = 0.0; dchiA = 0.0;
                dchiL_L = 0.0; dchiP_L  = 0.0; dchiA_L = 0.0;
                dchiL_R = 0.0; dchiP_R  = 0.0; dchiA_R = 0.0;

                dt   = 60.0;///days
                fdetRubin = 0.0;
                romanDetectable = false;

                for (int i = 0; i < M; ++i) {
                    Mpeak = s->magb[i] - 2.5 * std::log10(l->A0 * s->blend[i] + 1.0 - s->blend[i]);
//                        cout << "i=" << i << "  Mab=" << s->Mab[i] << "  Map=" << s->Map[i]
//                             << "  blend=" << s->blend[i] << "  Mpeak=" << Mpeak << endl;
                    if (i < 6) { // LSST ugrizy
                        if (Mpeak <= thre[i] and s->magb[i] > satu[i])    fdetRubin += 1.0;
                    } else {     // i == 6, Roman F146 — single band, no ">=2 filters" bar applies
                        if (Mpeak <= thre[i] and s->magb[i] > satu[i])    romanDetectable = true;
                    }
                }
                rubinDetectable = (fdetRubin > 1.0); // at least 2 of the 6 LSST bands

                testL = RandR(0.0, 1.0);
                testR = RandR(0.0, 1.0);
                // Independent draws per survey — reusing one draw for both would correlate
                // the Rubin-accept and Roman-accept decisions for no physical reason. Each
                // draw is weighted by that survey's OWN blend fraction (Step B2), replacing
                // the old single test <= s->blend[2] (LSST r-band only, for both surveys).
                acceptRubin = rubinDetectable and (testL <= s->blend[2]);
                acceptRoman = romanDetectable and (testR <= s->blend[6]);

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                if (acceptRubin or acceptRoman) { //detectable by Rubin (>=2 LSST bands) or by Roman (F146)
                    cout << "************** DETECTABLE!!!!!! ********" << endl;
                    s->nsdet[s->nums] += 1.0;
                    flagf = 1;
                    test  = RandR(0.0, 100.0);
    
                    if (test < 1.0 && save < 0 && IMnum == 1) {
                        initial = 20.0 * year;
                        save += 1;
                        flagm = 1;
                    }
    
                    gi = 0;
                    gi = 0; giR = 0;
                    for (double tim = float(0.0 * year - 100.0 - initial);  tim < float(10.0 * year + 100.0 + initial); tim = tim + dt) {
                        lightcurve(*s, *l, *as, tim);
                        Astar0   = double(s->ut0 * s->ut0 + 2.0) / std::sqrt(s->ut0 * s->ut0 * (s->ut0 * s->ut0 + 4.0)); //MAgnification equation
                        s->Astar = double(s->ut  * s->ut  + 2.0) / std::sqrt(s->ut  * s->ut  * (s->ut  * s->ut  + 4.0)); //MAgnification equation
                        As0      = double(Astar0   * s->blend[2] + 1.0 - s->blend[2]);
                        As1      = double(s->Astar * s->blend[2] + 1.0 - s->blend[2]); //LSST r-band
                        vs1      = double(s->mus1 + (s->def1c - def1p) / dt); //[mas/days]
                        vs2      = double(s->mus2 + (s->def2c - def2p) / dt); //[mas/days]
                        def1p    = s->def1c; //pervious
                        def2p    = s->def2c;

                        // Computed once per timestep (not just inside the LSST branch) since both
                        // instruments' astrometric chi-square terms need it, and it only depends on
                        // s->pos1b/pos2b/pos1c/pos2c, which lightcurve() already refreshed above.
                        trajm = std::sqrt(s->pos1b * s->pos1b + s->pos2b * s->pos2b); //stright + parallax
                        trajp = std::sqrt(s->pos1c * s->pos1c + s->pos2c * s->pos2c); //stright + parallax+lensing

                        if (flagm > 0) {
                            fil4 << std::fixed << std::setprecision(4)
                                 << tim      << " " << s->def1c << " " << s->def2c << " " << As0      << " " << As1 << " "
                                 << s->pos1b << " " << s->pos2b << " " << s->pos1c << " " << s->pos2c << " "
                                 << s->def1a << " " << s->def2a << " " << l->pos1  << " " << l->pos2  << "\n";
                        }
    
                        for (int i = 0; i < M; ++i) {
                            magni0[i] = s->magb[i] - 2.5 * std::log10(Astar0   * s->blend[i] + 1.0 - s->blend[i]);
                            magni[i]  = s->magb[i] - 2.5 * std::log10(s->Astar * s->blend[i] + 1.0 - s->blend[i]);
                        }
                        // Rubin's representative-band model magnitude (RUBIN_REF_BANDS, Bulge.h) --
                        // replaces the old hardcoded magni[2] (r-band) at the two use sites below.
                        // Reduces to exactly magni[2] when RUBIN_REF_BANDS = {2} (the default), since
                        // s->mbs[0]/s->fb[0] were built from the same combination in func_source.
                        double magniRubinRef = s->mbs[0] - 2.5 * std::log10(s->Astar * s->fb[0] + 1.0 - s->fb[0]);
                        sq = int(ls->ct[gi]);
                        sq  = int(ls->ct[gi]);
                        sqR = (nddR > 0 and giR < nddR) ? int(ro->ct[giR]) : -1;

                        // ---------------- LSST (ugrizy) ----------------
    
                        if (tim >= 0.0 and tim <= Tobs and tim >= ls->tim[int(ls->ct[0])] and tim <= ls->tim[int(ls->ct[ndd - 1])] and
                            gi < ndd and sq >= 0 and sq <= static_cast<int>(Nl) and tim >= ls->tim[sq]) {
    
                            fi = int(ls->filter[sq]);
    
                            if (magni[fi] >= satu[fi] and magni[fi] <= thre[fi]) {
                                errg = errlsstM(magni[fi], int(fi), double(ls->sig5[sq])); //[mag]
                                errs = errlsstA(*ls, magniRubinRef); ///[mas]
    
                                deltaA = std::fabs(std::pow(10.0, -0.4 * errg) - 1.0) * (s->blend[fi] * s->Astar + 1.0 - s->blend[fi]);
                                magnio = magni[fi] + RandN(errg, 3.0);
    
                                chi1 += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg)); //real
                                chi2 += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg)); //real no parallax
                                chi3 += std::fabs((magnio - s->magb[fi]) * (magnio - s->magb[fi]) / (errg * errg)); //baseline
                                chi1_L += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg));
                                chi2_L += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg));
                                chi3_L += std::fabs((magnio - s->magb[fi]) * (magnio - s->magb[fi]) / (errg * errg));

                                sil  = RandN(errs * std::sqrt(2.0), 3.0);
                                sil2 = RandN(errs, 3.0);

                                chi1a += std::fabs((trajp + sil - trajp) * (trajp + sil - trajp) / (errs * errs * 2.0)); //real trajectory
                                chi2a += std::fabs((trajp + sil - trajm) * (trajp + sil - trajm) / (errs * errs * 2.0)); //real without lensing
                                chi3a += std::fabs( sil2  * sil2 / (errs * errs));
                                chi1a_L += std::fabs((trajp + sil - trajp) * (trajp + sil - trajp) / (errs * errs * 2.0));
                                chi2a_L += std::fabs((trajp + sil - trajm) * (trajp + sil - trajm) / (errs * errs * 2.0));
                                chi3a_L += std::fabs( sil2  * sil2 / (errs * errs));

                                CHECK(ndw < coun); // array-bounds guard — see note on `coun` sizing
                                l->timn[ndw] = tim;
                                l->magn[ndw] = magniRubinRef; // Rubin representative-band model value (RUBIN_REF_BANDS)
                                l->errm[ndw] = errg;
                                l->soux[ndw] = s->pos1c;
                                l->souy[ndw] = s->pos2c;
                                l->erra[ndw] = errs;
                                l->tele[ndw] = 0; // 0 = LSST
    
                                flag2 = 0.0;
                                if (std::fabs(magnio - s->magb[fi]) > std::fabs(3.0 * errg))    flag2 = 1.0;
                                if (ndw_L > 2 and float(flag0 + flag1 + flag2) > 2.0)   flag_det = 1;
                                flag0 = flag1;
                                flag1 = flag2;
                                flag2_L = flag2; // identical test; kept as an explicit Rubin-labeled copy
                                if (ndw_L > 2 and float(flag0_L + flag1_L + flag2_L) > 2.0)   flag_det_L = 1;
                                flag0_L = flag1_L;
                                flag1_L = flag2_L;
    
                                if (flagm > 0) {
                                    fil5 << std::fixed << std::setprecision(4)
                                         << tim     << "  "
                                         << std::setprecision(6) << As1 + RandN(deltaA, 3.0)    << "  "
                                         << deltaA  << "  "
                                         << std::setprecision(4) << s->pos1c + RandN(errs, 3.0) << "  "
                                         << s->pos2c + RandN(errs, 3.0) << "  "
                                         << s->def1c + RandN(errs, 3.0) << "  "
                                         << s->def2c + RandN(errs, 3.0) << "  "
                                         << errs    << "  "
                                         << 1.0     << "  "
                                         << int(fi) << "\n";
                                }
    
                                CHECK(sq >= 0);
                                CHECK(sq <= int(Nl - 1));
                                CHECK(fi >= 0);
                                CHECK(fi <= 5);
                                CHECK(errs >= 0.0);
                                CHECK(errg >= 0.0);
                                CHECK(deltaA >= 0.0);
                                CHECK(gi <= ndd);
    
                                s->errM += deltaA;
                                s->errA += errs;
                                vsave += std::sqrt(vs1 * vs1 + vs2 * vs2);

                                ndw += 1;
                                ndw_L += 1;
                            }//magnitude limit
                            gi += 1;
                        }
                        // ---------------- Roman (F146) — new ----------------
                        // Roman now feeds its own chi1_R/chi2_R/chi3_R (+ astrometric _R,
                        // placeholder error — see errsR below) and the joint chi1/chi2/chi3,
                        // alongside the Rubin-only _L versions built in the LSST branch above.
                        // See Step B1 of JOINT_FIT_REFACTOR_PLAN.md. s->errM/s->errA remain
                        // Rubin-only by design (ORIENTATION.md: they're specifically the
                        // LSST-only running-error accumulators, unrelated to detection).
                        if (nddR > 0 and tim >= 0.0 and tim <= Tobs and
                            tim >= ro->tim[int(ro->ct[0])] and tim <= ro->tim[int(ro->ct[nddR - 1])] and
                            giR < nddR and sqR >= 0 and sqR <= static_cast<int>(NlRoman) and tim >= ro->tim[sqR]) {

                            constexpr int fiR = 6; // F146 — the only Roman band modeled so far

                            if (magni[fiR] >= satu[fiR] and magni[fiR] <= thre[fiR]) {
                                // TODO(Ali): confirm this matches how sigma_roman.txt / ro->mag,ro->err
                                // are meant to be interpolated (see errRomanM stub near matchVisibleEpochs).
                                errgR = errRomanM(*ro, magni[fiR]); //[mag]

                                magnioR = magni[fiR] + RandN(errgR, 3.0);
                                chi1 += std::fabs((magnioR -   magni[fiR]) * (magnioR -   magni[fiR]) / (errgR * errgR));
                                chi2 += std::fabs((magnioR -  magni0[fiR]) * (magnioR -  magni0[fiR]) / (errgR * errgR));
                                chi3 += std::fabs((magnioR - s->magb[fiR]) * (magnioR - s->magb[fiR]) / (errgR * errgR));
                                chi1_R += std::fabs((magnioR -   magni[fiR]) * (magnioR -   magni[fiR]) / (errgR * errgR));
                                chi2_R += std::fabs((magnioR -  magni0[fiR]) * (magnioR -  magni0[fiR]) / (errgR * errgR));
                                chi3_R += std::fabs((magnioR - s->magb[fiR]) * (magnioR - s->magb[fiR]) / (errgR * errgR));

                                // TODO(Ali): PLACEHOLDER astrometric error — no Roman/F146 astrometric
                                // error model exists yet. Stands in with LSST's astrometric-error curve
                                // (errlsstA) evaluated at Roman's own magnitude/epoch, so it's at least
                                // deterministic and epoch-correct (not a stale `errs` left over from
                                // whichever epoch last set it). Per JOINT_FIT_REFACTOR_PLAN.md's
                                // Deferred section ("Roman astrometric error model decision" — the
                                // ~100 mas FWHM and gamma value for the F146~22 transition), replace
                                // this call with a real errRomanA()-style function reading a real F146
                                // astrometric-error dataset (mirroring errlsstA/sigmaA_LSST.txt) before
                                // relying on chi1a_R/chi2a_R/chi3a_R, dchiA_R, or the astrometric Fisher
                                // matrix (inputB/inverB, Ny params tetE/mus1/mus2/piE — Phase C, Step C4)
                                // for any real conclusion. Logged in OPEN_ITEMS.md.
                                errsR = errlsstA(*ls, magni[fiR]); //[mas]
                                silR  = RandN(errsR * std::sqrt(2.0), 3.0);
                                sil2R = RandN(errsR, 3.0);
                                chi1a += std::fabs((trajp + silR - trajp) * (trajp + silR - trajp) / (errsR * errsR * 2.0));
                                chi2a += std::fabs((trajp + silR - trajm) * (trajp + silR - trajm) / (errsR * errsR * 2.0));
                                chi3a += std::fabs( sil2R * sil2R / (errsR * errsR));
                                chi1a_R += std::fabs((trajp + silR - trajp) * (trajp + silR - trajp) / (errsR * errsR * 2.0));
                                chi2a_R += std::fabs((trajp + silR - trajm) * (trajp + silR - trajm) / (errsR * errsR * 2.0));
                                chi3a_R += std::fabs( sil2R * sil2R / (errsR * errsR));

                                flag2_R = 0.0;
                                if (std::fabs(magnioR - s->magb[fiR]) > std::fabs(3.0 * errgR))    flag2_R = 1.0;
                                if (ndw_R > 2 and float(flag0_R + flag1_R + flag2_R) > 2.0)   flag_det_R = 1;
                                flag0_R = flag1_R;
                                flag1_R = flag2_R;

                                CHECK(ndw < coun); // array-bounds guard — see note on `coun` sizing
                                l->timn[ndw] = tim;
                                l->magn[ndw] = magni[fiR]; // F146 magnitude — NOT magni[2]; FisherM's
                                                            // tt==1 branch compares against s.mbs[1]/s.fb[1],
                                                            // which are F146-based, so the reference point
                                                            // recorded here must be F146 too.
                                l->errm[ndw] = errgR;
                                l->soux[ndw] = s->pos1c;
                                l->souy[ndw] = s->pos2c;
                                // TODO(Ali): Roman has its own astrometric precision (it's the whole
                                // point of the astrometric-microlensing science case). Reusing LSST's
                                // `errs` here is a placeholder so the array has *something* consistent
                                // in it; replace with a Roman astrometric error model when ready.
                                l->erra[ndw] = errs;
                                l->tele[ndw] = 1; // 1 = Roman/F146

                                CHECK(sqR >= 0);
                                CHECK(sqR <= int(NlRoman - 1));
                                CHECK(errgR >= 0.0);
                                CHECK(giR <= nddR);

                                ndw += 1;
                                ndw_R += 1;
                            }//magnitude limit
                            giR += 1;
                        }

                        // ---------------- Adaptive step size ----------------
                        // dt must respect whichever instrument's *next* epoch is sooner.
                        // Sizing dt to LSST's cadence alone (the original behavior) would
                        // silently step right over Roman's much denser epochs.
                        if ( tim < -50.0 or tim > (10.0 * year + 50.0)) {
                            dt = 60.0; //days
                        }
    
                        else {
                            bool lsstInWindow  = (tim >= ls->tim[int(ls->ct[0])] and tim <= ls->tim[int(ls->ct[ndd - 1])]);
                            bool romanInWindow = (nddR > 0 and tim >= ro->tim[int(ro->ct[0])] and tim <= ro->tim[int(ro->ct[nddR - 1])]);

                            if (lsstInWindow) {
                                if (gi > 0 and gi < ndd) cade = float(ls->tim[int(ls->ct[gi])] - ls->tim[int(ls->ct[gi - 1])]); //days
                                else cade = minc;
                            } else {
                                cade = 3.0; //days — matches the original "outside LSST window" fallback
                            }

                            if (romanInWindow) {
                                if (giR > 0 and giR < nddR) cadeR = float(ro->tim[int(ro->ct[giR])] - ro->tim[int(ro->ct[giR - 1])]); //days
                                else cadeR = mincR;
                            } else {
                                cadeR = cade; // Roman not active right now — don't let it constrain dt
                            }

                            dt = double(std::min(cade, cadeR));
                        }
    
                    }//end of loop time
                    if (flagm > 0) {
                        fil4.close();
                        fil5.close();
                    }
                }// end of visible star

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH 
                for (int i = 0; i < nq; ++i) co->resu[i] = -1.0;
                // Same sentinel discipline for the per-survey results: an event that never
                // reaches FisherM must not inherit the previous event's sigmas.
                for (int q = 0; q < NSURV; ++q) {
                    co->okA[q] = 0; co->okB[q] = 0; co->nepochA[q] = 0;
                    co->condA[q] = -1.0; co->condB[q] = -1.0;
                    for (int k = 0; k < Nx; ++k) co->Era[q][k] = -1.0;
                    for (int k = 0; k < Ny; ++k) co->Erb[q][k] = -1.0;
                }
    
                FFG[0] = 0;
                FFG[1] = 0;
                FFG[2] = 0;
                detL = 0; detR = 0; detJ = 0;

                //cout << "flagf: " << flagf << "ndw: " << ndw << endl;

                if (flagf == 0 or ndw <= 2) {
                    errg    = errlsstM(s->magb[2], 2, double(24.43)); //r-band
                    s->errA = errlsstA(*ls, s->magb[2]); //r-band
                    s->errM = std::fabs(std::pow(10.0, - 0.4 * errg) - 1.0); //r-band
                    dchiL = 0.0;
                    dchiP = 0.0;
                    dchiA = 0.0;
                    dchiL_L = 0.0; dchiP_L = 0.0; dchiA_L = 0.0;
                    dchiL_R = 0.0; dchiP_R = 0.0; dchiA_R = 0.0;
                    vsave = s->mus;
                }

                if (flagf > 0 and ndw > 2) { //if star is visible
                    cout << "************** DETECTABLE!!!!!! ********" << endl;
                    icon +=1;
                    vsave   = double(vsave   / (ndw + 0.000065645));
                    s->errM = double(s->errM / (ndw + 0.000065645));
                    s->errA = double(s->errA / (ndw + 0.000065645));
                    dchiL   = std::fabs(chi3  - chi1);  //lensing_effect
                    dchiP   = std::fabs(chi2  - chi1);  //parallax_effect
                    dchiA   = std::fabs(chi2a - chi1a); //deflection_effect
                    dchiL_L = std::fabs(chi3_L  - chi1_L);
                    dchiP_L = std::fabs(chi2_L  - chi1_L);
                    dchiA_L = std::fabs(chi2a_L - chi1a_L);
                    dchiL_R = std::fabs(chi3_R  - chi1_R);
                    dchiP_R = std::fabs(chi2_R  - chi1_R);
                    dchiA_R = std::fabs(chi2a_R - chi1a_R);

                    // Three independent detection tests. detL/detR use each instrument's own
                    // epoch count (ndw_L/ndw_R) and run-test result — mixing in the joint ndw
                    // here would let Roman's dense epochs silently raise the bar for a purely
                    // Rubin-driven signal (and vice versa). detJ requires a persistent run in
                    // EITHER instrument's own cadence (flag_det_L or flag_det_R) rather than an
                    // interleaved run across two different cadences, which wouldn't mean anything
                    // (see Step B1 teaching brief). FFG[0] — the boolean that gates FisherM below,
                    // preserving the existing detection-then-Fisher ordering (plan §0.4) — is the
                    // union of all three, so an event Rubin alone clearly detects is never dropped
                    // from characterization just because the joint-scaled threshold happens to miss.
                    if (s->FWHM < Tobs and dchiL_L > float(2.0 * ndw_L) and flag_det_L > 0 and ndw_L > 10) detL = 1;
                    if (s->FWHM < Tobs and dchiL_R > float(2.0 * ndw_R) and flag_det_R > 0 and ndw_R > 10) detR = 1;
                    if (s->FWHM < Tobs and dchiL   > float(2.0 * ndw)   and (flag_det_L > 0 or flag_det_R > 0) and ndw > 10) detJ = 1;

                    if (detL or detR or detJ) { //lensing — detected by Rubin, Roman, or the joint test
                        FFG[0] = 1; //Lensing
                        nlens += FFG[0];
                        FisherM(*s, *l, *as, *co, ndw);

                        if (co->flagi > 0) {
                            nerr += 1.0;
                            ErrorCal(*co, *l, *s);

                            std::ofstream fil0_append(fnLDt, std::ios::app);
                            fil0_append << std::fixed << std::setprecision(5)
                                        << l->Ml       << " " << l->Dl       << " " << l->mul       << " " << s->fb[0]     << " " << s->mbs[0]   << " "
          <<std::setprecision(7)           << l->tE       << " " << l->murel    << " " << l->u0        << " " << s->lon       << " " << s->lat      << " "
                                        << l->piE      << " " << l->tetE     << " " << co->resu[1]  << " " << co->resu[2]  << " " << co->resu[3] << " "
                                        << co->resu[5] << " " << co->resu[9] << " " << co->resu[10] << " " << co->resu[13] << "\n";
                            fil0_append.close();
                        }
                    }

                    // Five-way detection label — diagnostics for Step B1's acceptance check.
                    // Not yet persisted to EventRecord; that's Step D1's job.
                    if      (detL and detR)  nDetBoth      += 1;
                    else if (detL)           nDetRubinOnly += 1;
                    else if (detR)           nDetRomanOnly += 1;
                    else if (detJ)           nDetJointOnly += 1;
                    else                     nDetNeither   += 1;

                    gg = FunctE(*l);
                    //ss = int(FuncMl(*l));
                    //qq = int(FuncPi(*l));
                    //ww = int(Funcu0(*l));
                    //vv = int(FuncMu(*l));
                    //zz = int(FuncMb(*l, s->Map[2]));
                    //pp = int(FuncFb(*l, s->blend[2]));
                }
//            CHECK(gg >= 0);

            records.push_back(EventRecord{
                icon, static_cast<int>(FFG[0]),
                l->tE, l->RE/AU, l->piE, l->tetE, l->Vt, l->u0, l->Ml,
                s->opt*1.0e6, l->Dl, s->Ds, l->vl, s->vs,
                s->mbs[0], s->fb[0],
                gg, static_cast<int>(l->struc),
                s->FWHM/year, vsave/s->mus, l->DeltaT/s->errA, l->murel*year,
                co->resu[0], co->resu[1], co->resu[2], co->resu[3], co->resu[5],
                co->resu[9], co->resu[10], co->resu[13], co->resu[14],
                s->Map[2], s->nsbl[2],
                co->flagi, s->Ai[2],
                ndw_L, ndw_R,
                detL, detR, detJ,
                co->okA[SJOINT], co->okA[SRUBIN], co->okA[SROMAN],
                co->Era[SJOINT][1], co->Era[SRUBIN][1], co->Era[SROMAN][1],   //sigma(tE)
                co->Era[SJOINT][3], co->Era[SRUBIN][3], co->Era[SROMAN][3],   //sigma(piE)
                co->Erb[SJOINT][0], co->Erb[SRUBIN][0], co->Erb[SROMAN][0],   //sigma(tetE)
                synergyClass(*co),
                co->condA[SJOINT], co->condA[SRUBIN], co->condA[SROMAN]
            });
   
            filg_in.open(testf, std::ios::app);
            filg_in << icon           << " " << FFG[0]         << " " << l->tE                      << " "
                    << l->RE / AU     << " " << l->piE         << " " << l->tetE                    << " "
                    << l->Vt          << " " << l->u0          << " " << l->Ml                      << " " 
                    << s->opt * 1.0e6 << " " << l->Dl          << " " << s->Ds                      << " "
                    << l->vl          << " " << s->vs          << " " << s->mbs[0]                  << " "
                    << s->fb[0]       << " " << gg             << " " << static_cast<int>(l->struc) << " "
                    << s->FWHM / year << " " << vsave / s->mus << " " << l->DeltaT / s->errA        << " " << l->murel * year << " "
                    << co->resu[0]    << " " << co->resu[1]    << " " << co->resu[2]                << " " 
                    << co->resu[3]    << " " << co->resu[5]    << " "
                    << co->resu[9]    << " " << co->resu[10]   << " " << co->resu[13]               << " " << co->resu[14]    << " "
                    << s->Map[2]      << " " << s->nsbl[2]     << " " << co->flagi                  << " " << s->Ai[2]        << " "
                    // per-survey bookkeeping (Step C5), appended so existing column indices hold
                    << ndw_L << " " << ndw_R << " "
                    << detL  << " " << detR  << " " << detJ << " "
                    << co->okA[SJOINT] << " " << co->okA[SRUBIN] << " " << co->okA[SROMAN] << " "
                    << co->Era[SJOINT][1] << " " << co->Era[SRUBIN][1] << " " << co->Era[SROMAN][1] << " "
                    << co->Era[SJOINT][3] << " " << co->Era[SRUBIN][3] << " " << co->Era[SROMAN][3] << " "
                    << co->Erb[SJOINT][0] << " " << co->Erb[SRUBIN][0] << " " << co->Erb[SROMAN][0] << " "
                    << synergyClass(*co) << " "
                    << co->condA[SJOINT] << " " << co->condA[SRUBIN] << " " << co->condA[SROMAN] << "\n";
            filg_in.close();
//          


            if (flagm > 0 && IMnum == 1) {
                std::ofstream fil1_append("./files/MONTLMC/files/BHLSSTMONTS.dat", std::ios::app);

                fil1_append << save            << " " << std::fixed << std::setprecision(5)
                            << s->lat          << " " << s->lon         << " "
                            << static_cast<int>(l->struc)               << " " << l->Ml           << " "
                            << l->Dl           << " " << l->vl          << " "
                            << static_cast<int>(s->struc)               << " " << s->cl           << " "
                            << s->Ds           << " " << s->logT        << " " << s->mass         << " "
                            << s->vs           << " " << s->Mab[2]      << " " << s->Mab[6]       << " "
                            << s->Map[2]       << " " << s->Map[6]      << " "
                            << s->mbs[0]       << " " << s->mbs[1]      << " " << s->fb[0]        << " "
                            << s->fb[1]        << " " << s->nsbl[2]     << " "
                            << s->nsbl[6]      << " " << s->Ai[2]       << " " << s->Ai[6]        << " "
                            << std::setprecision(7)
                            << l->tE           << " " << l->RE / AU     << " " << l->t0           << " "
                            << l->mul          << " " << l->Vt          << " " << l->u0           << " "
                            << s->opt*1.0e6    << " " << l->tetE        << " "
                            << s->mus1         << " " << s->mus2        << " " << s->xi           << " "
                            << l->mul1         << " " << l->mul2        << " " << l->piE          << " "
                            << s->errM         << " " << s->errA        << " "
                            << flagf           << " " << flag_det       << " " << dchiL           << " "
                            << dchiP           << " " << dchiA          << " " << ndw             << " "
                            << FFG[0]          << " " << FFG[1]         << " " << FFG[2]          << " "
                            << s->Rostart      << " " << s->Nstart      << " " << l->mi1          << " "
                            << l->mi2          << " " << vsave          << " " << s->FWHM         << " "
                            << chi3a           << " " << nri            << " " << nde             << " "
                            << l->betal  * RAa << " " << l->betas * RAa << " " << l->deltal * RAa << " "
                            << l->deltas * RAa << " "
                            << l->deltao * RAa << " "
                            << co->resu[0]     << " " << co->resu[1]    << " " << co->resu[2]     << " "
                            << co->resu[3]     << " " << co->resu[5]    << " "
                            << co->resu[9]     << " " << co->resu[10]   << " " << co->resu[13]    << " "
                            << co->resu[14]    << "\n";
                
                fil1_append.close();
            }
            //cout << "** End of saving in the file *********" << save << endl;
            //cout << "icon: " << icon << "\tnlens: " << nlens << "\tnerr: " << nerr << endl;
//            } while (icon < 850 or nlens < 150 or nerr < 2.0); //end of loop icon
            } while (icon < 20 or nlens < 5 or nerr < 1.0); //end of loop icon
   
            for (int i = 0; i <= GG; ++i) {
                l->NstE[i] += l->nstE[i];
                l->NdtE[i] += l->ndtE[i];
                l->ndtE[i]  = double(l->ndtE[i] / (l->nstE[i] + eps)); // [0.0, 1.0]

                fil2 << std::fixed << std::setprecision(4)
                     << l->tEs[i]  << " "
                     << double(l->NdtE[i]  * 100.0 / (l->NstE[i] + eps)) << " "
                     << l->Mls[i]  << " "
                     << double(l->NdMl[i]  * 100.0 / (l->NsMl[i] + eps)) << " "
                     << l->pis[i]  << " "
                     << double(l->Ndpi[i]  * 100.0 / (l->Nspi[i] + eps)) << " "
                     << l->u0s[i]  << " "
                     << double(l->Ndu0[i]  * 100.0 / (l->Nsu0[i] + eps)) << " "
                     << l->mbs[i]  << " "
                     << double(l->Ndmb[i]  * 100.0 / (l->Nsmb[i] + eps)) << " "
                     << l->fbs[i]  << " "
                     << double(l->Ndfb[i]  * 100.0 / (l->Nsfb[i] + eps)) << " "
                     << l->mus[i]  << " "
                     << double(l->Ndmu[i]  * 100.0 / (l->Nsmu[i] + eps)) << " "
                     << double(l->Nhalo[1] * 100.0 /  l->Nhalo[0])       << " "
                     << double(l->Nself[1] * 100.0 /  l->Nself[0])       << "\n";

                fil2b << std::fixed  << std::setprecision(1)
                      << l->NdtE[i]  << " " << l->NstE[i]  << " "
                      << l->NdMl[i]  << " " << l->NsMl[i]  << " "
                      << l->Ndpi[i]  << " " << l->Nspi[i]  << " "
                      << l->Ndu0[i]  << " " << l->Nsu0[i]  << " "
                      << l->Ndmb[i]  << " " << l->Nsmb[i]  << " "
                      << l->Ndfb[i]  << " " << l->Nsfb[i]  << " "
                      << l->Ndmu[i]  << " " << l->Nsmu[i]  << " "
                      << l->Nhalo[1] << " " << l->Nhalo[0] << " "
                      << l->Nself[1] << " " << l->Nself[0] << "\n";
            }
   
            s->nstart = 0.0;
            for (int i = 0; i < Num; ++i) {
                test = double(s->nsdet[i]   / (s->nssim[i] + eps));
                s->nstart  += s->Rostari[i] * (s->Nstart   / s->Rostart) * test;
            } //Number/deg^{2}

            for (int i = 0; i < 2; ++i){
                tE[i]   = 0.0, RE[i]  = 0.0, piE[i]  = 0.0, tetE[i]  = 0.0; Vt[i]    = 0.0;
                u0[i]   = 0.0, Ml[i]  = 0.0, opd[i]  = 0.0; Dl[i]    = 0.0, Ds[i]    = 0.0;
                vl[i]   = 0.0, vs[i]  = 0.0, mbs[i]  = 0.0, fb[i]    = 0.0; numd[i]  = 0.0;
                fwhm[i] = 0.0; vsn[i] = 0.0; DelT[i] = 0.0; Struc[i] = 0.0; murel[i] = 0.0;
                Map[i]  = 0.0; nbl[i] = 0.0; Ext[i]  = 0.0;
            }

            EffiL = 0.0; EFF   = 0.0; Gamma  = 0.0; Neven = 0.0; EffiD = 0.0;
            ErtE  = 0.0; ErpiE = 0.0; ErtetE = 0.0; Erml  = 0.0; Erfb  = 0.0;
            Erdl  = 0.0; Ermul = 0.0; Ermus  = 0.0; Eru0  = 0.0;
  
    int tempStruc;
    for (const auto& r : records) {
//        counter = r.counter;
        flagL = r.flagL;
        l->tE = r.tE;   l->RE = r.RE;   l->piE = r.piE;  l->tetE = r.tetE;
        l->Vt = r.Vt;   l->u0 = r.u0;   l->Ml  = r.Ml;
        s->opt = r.opt; l->Dl = r.Dl;   s->Ds  = r.Ds;   l->vl = r.vl; s->vs = r.vs;
        mbase = r.mbase; fblend = r.fblend; gg = r.gg; tempStruc = r.struc;
        s->FWHM = r.FWHM; vsave = r.vsave; l->DeltaT = r.DeltaT; l->murel = r.murel;
        co->resu[0]=r.resu0;  co->resu[1]=r.resu1;   co->resu[2]=r.resu2;
        co->resu[3]=r.resu3;  co->resu[5]=r.resu5;   co->resu[9]=r.resu9;
        co->resu[10]=r.resu10; co->resu[13]=r.resu13; co->resu[14]=r.resu14;
        s->Map[2]=r.Map2; s->nsbl[2]=r.nsbl2; co->flagi=r.flagi; s->Ai[2]=r.Ai2;

        l->struc = static_cast<GalacticComponent>(tempStruc);

        // ------------------ Update cumulative sums ------------------
        Struc[0] += 1.0;
        tE[0]    += l->tE;    RE[0]  += l->RE;     piE[0]  += l->piE;
        tetE[0]  += l->tetE;  Vt[0]  += l->Vt;     u0[0]   += l->u0;
        Ml[0]    += l->Ml;    opd[0] += s->opt;    Dl[0]   += l->Dl;
        Ds[0]    += s->Ds;    vl[0]  += l->vl;     vs[0]   += s->vs;
        mbs[0]   += mbase;    fb[0]  += fblend;    numd[0] += 1.0;
        fwhm[0]  += s->FWHM;  vsn[0] += vsave;     DelT[0] += l->DeltaT;
        murel[0] += l->murel; Map[0] += s->Map[2]; nbl[0]  += s->nsbl[2];
        Ext[0]   += s->Ai[2];

        // ------------------ Flagged case ---------------------------
        if (flagL > 0) {
            Struc[1] += 1.0;
            tE[1]    += l->tE;    RE[1]  += l->RE;     piE[1]  += l->piE;
            tetE[1]  += l->tetE;  Vt[1]  += l->Vt;     u0[1]   += l->u0;
            Ml[1]    += l->Ml;    opd[1] += s->opt;    Dl[1]   += l->Dl;
            Ds[1]    += s->Ds;    vl[1]  += l->vl;     vs[1]   += s->vs;
            mbs[1]   += mbase;    fb[1]  += fblend;    numd[1] += 1.0;
            fwhm[1]  += s->FWHM;  vsn[1] += vsave;     DelT[1] += l->DeltaT;
            murel[1] += l->murel; Map[1] += s->Map[2]; nbl[1]  += s->nsbl[2];
            Ext[1]   += s->Ai[2];

            EFF += static_cast<double>(l->ndtE[gg] / (l->tE / year)); // 1/years

            if (co->flagi > 0) {
                Eru0  += co->resu[0];  ErtE   += co->resu[1];  Erfb  += co->resu[2];
                ErpiE += co->resu[3];  ErtetE += co->resu[5];  Erml  += co->resu[9];
                Erdl  += co->resu[10]; Ermul  += co->resu[13]; Ermus += co->resu[14];
            }
        }
    }

    for (int i = 0; i < 2; ++i) { // What is the purpose of this block?
        tE[i]    = double(tE[i]            / (numd[i] + eps)) / year; //[years]  
        RE[i]    = double(RE[i]            / (numd[i] + eps));//[AU]  
        piE[i]   = double(piE[i]           / (numd[i] + eps));//[]     
        tetE[i]  = double(tetE[i]          / (numd[i] + eps));//[mas]  
        Vt[i]    = double(Vt[i]            / (numd[i] + eps));//[km/s]  
        u0[i]    = double(u0[i]            / (numd[i] + eps));//[]  
        Ml[i]    = double(Ml[i]            / (numd[i] + eps));//[Msun]  
        opd[i]   = double(opd[i]           / (numd[i] + eps)) * u0m * u0m;//[] x10^{6}
        Dl[i]    = double(Dl[i]            / (numd[i] + eps));//[kpc]  
        Ds[i]    = double(Ds[i]            / (numd[i] + eps));//[kpc]  
        vl[i]    = double(vl[i]            / (numd[i] + eps));//[km/s]  
        vs[i]    = double(vs[i]            / (numd[i] + eps));//[km/s]  
        mbs[i]   = double(mbs[i]           / (numd[i] + eps));//[mag]  
        fb[i]    = double(fb[i]            / (numd[i] + eps));//[]
        fwhm[i]  = double(fwhm[i]          / (numd[i] + eps));//[years]
        vsn[i]   = double(vsn[i]           / (numd[i] + eps));//[]
        DelT[i]  = double(DelT[i]          / (numd[i] + eps));//[]
        murel[i] = double(murel[i]         / (numd[i] + eps));//[mas/year]
        Map[i]   = double(Map[i]           / (numd[i] + eps));//[mag]
        nbl[i]   = double(nbl[i]           / (numd[i] + eps));//[]
        Ext[i]   = double(Ext[i]           / (numd[i] + eps));//[mag]
        Struc[i] = double(Struc[i] * 100.0 / (numd[i] + eps));
    }

    EFF = double(EFF / (numd[1] + eps));//1/[years]
    Gamma = double(2.0 / M_PI * opd[0] * 1.0e-6 * EFF) / u0m; // 1/[star*year]  
    EffiL = double(numd[1] * 100.0 / (numd[0] + eps)); // probability of detecting lensing  
    EffiD = double(numd[0] * 100.0 / (nsim    + eps)); // probability of detecting stars  
    Neven = double(s->nstart * Gamma * 10.0);//deg^{-2}
   
    Eru0   = double(Eru0   / (nerr + eps));  
    Erfb   = double(Erfb   / (nerr + eps));  
    ErtE   = double(ErtE   / (nerr + eps));  
    ErpiE  = double(ErpiE  / (nerr + eps));  
    ErtetE = double(ErtetE / (nerr + eps));  
    Erml   = double(Erml   / (nerr + eps));  
    Erdl   = double(Erdl   / (nerr + eps));  
    Ermul  = double(Ermul  / (nerr + eps));  
    Ermus  = double(Ermus  / (nerr + eps));  
   
    //Maps
    fil3 << std::fixed << std::setprecision(6)
         << tE[0]    << " " << tE[1]    << " "
         << RE[0]    << " " << RE[1]    << " "
         << piE[0]   << " " << piE[1]   << " "
         << tetE[0]  << " " << tetE[1]  << " "
         << Vt[0]    << " " << Vt[1]    << " "
         << u0[0]    << " " << u0[1]    << " "
         << Ml[0]    << " " << Ml[1]    << " "
         << opd[0]   << " " << opd[1]   << " "
         << Dl[0]    << " " << Dl[1]    << " "
         << Ds[0]    << " " << Ds[1]    << " "
         << vl[0]    << " " << vl[1]    << " "
         << vs[0]    << " " << vs[1]    << " "
         << mbs[0]   << " " << mbs[1]   << " "
         << fb[0]    << " " << fb[1]    << " "
         << fwhm[0]  << " " << fwhm[1]  << " "
         << vsn[0]   << " " << vsn[1]   << " "
         << DelT[0]  << " " << DelT[1]  << " "
         << Struc[0] << " " << Struc[1] << " "
         << murel[0] << " " << murel[1] << " "
         << Map[0]   << " " << Map[1]   << " "
         << nbl[0]   << " " << nbl[1]   << " "
         << Ext[0]   << " " << Ext[1]   << " ";

    fil3 << std::setprecision(8)
         << EffiD      << " " << EffiL        << " "
         << std::log10(EFF) << " " << std::log10(Gamma) << " " << std::log10(Neven) << " "
         << Eru0       << " " << ErtE         << " " << Erfb         << " " << ErpiE << " " << ErtetE << " "
         << Erml       << " " << Erdl         << " " << Ermul        << " " << Ermus << " ";

    fil3 << std::setprecision(1)
         << nsim              << " " << numd[0]          << " " << numd[1] << " "
         << nerr              << " " << nri              << " " << nde     << " "
         << std::log10(s->Rostart) << " " << std::log10(s->Nstart) << " " << std::log10(s->nstart)
         << "\n";
   
    cout << "nsim:  "  << nsim    << "\t Ndetected:  " << icon    << "\t Nlensing:  " << nlens << "\t NError:  " << nerr << endl;
    cout << "Detection label counts — neither: " << nDetNeither
         << "  Rubin-only: " << nDetRubinOnly
         << "  Roman-only: " << nDetRomanOnly
         << "  joint-only: " << nDetJointOnly
         << "  both: "       << nDetBoth << endl;
    cout << "numd0:  " << numd[0] << "\t numd1:  "     << numd[1] << endl;
    cout << "EFF:  "   << EFF     << "\t Gamma:  "     << Gamma   << "\t Neven:  "    << Neven << endl;
    cout << "l.tE:  "  << l->tE   << "\t l.Ml:  "      << l->Ml   << endl;
   
    CHECK(EFF > 0.0);
    CHECK(Gamma > 0.0);
    CHECK(EffiL > 0.0);
    CHECK(EffiD > 0.0);
    CHECK(Neven > 0.0);
    
    CHECK(DelT[0] > 0.0);
    CHECK(DelT[1] > 0.0);
    CHECK(vsn[0] > 0.0);
    CHECK(vsn[1] > 0.0);
    CHECK(fwhm[0] > 0.0);
    CHECK(fwhm[1] > 0.0);
    
    CHECK(fb[0] > 0.0);
    CHECK(fb[0] <= 1.0);
    CHECK(fb[1] > 0.0);
    CHECK(fb[1] <= 1.0);
    
    CHECK(mbs[0] > 0.0);
    CHECK(mbs[1] > 0.0);
    
    CHECK(vs[0] > 0.0);
    CHECK(vs[1] > 0.0);
    CHECK(vl[0] > 0.0);
    CHECK(vl[1] > 0.0);
    
    CHECK(Ds[0] > 0.0);
    CHECK(Ds[1] > 0.0);
    CHECK(Dl[0] > 0.0);
    CHECK(Dl[1] >= 0.0);
    
    CHECK(piE[1] >= 0.0);
    CHECK(tetE[1] > 0.0);
    
    CHECK(murel[0] > 0.0);
    CHECK(murel[1] > 0.0);
    
    CHECK(ErtE > 0.0);
    CHECK(ErpiE > 0.0);
    CHECK(ErtetE > 0.0);
    CHECK(Erml > 0.0);
    CHECK(Erdl > 0.0);
    CHECK(Ermul > 0.0);
    CHECK(Ermus > 0.0);
    CHECK(Eru0 > 0.0);
    CHECK(Erfb > 0.0);
    
    CHECK(nerr != 0.0);
    CHECK(numd[0] != 0.0);
    CHECK(numd[1] != 0.0);
    CHECK(nsim != 0.0);
    
    CHECK(icon == numd[0]);
    CHECK(numd[1] == nlens);
     
    cout << "==============================================================" << endl;  
  
   }}//right_accention and declinaton
   //fclose(filh);  
   return(0);
}
#endif // FISHER_FIXTURE_BUILD

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Fisher and Covariance matrixes             //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void FisherM(source & s, lens & l, astromet & as,  covarian & co, int ndw)
{
    co.flagi =+ 1;
    int tt; // tev;

// (K + H + J)/3 = F146
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///        Photometry                                  ///
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    co.Delta1[0] = 0.1507586576;//u0
    co.Delta1[1] = l.tE * 0.254674;//tE [days]
    co.Delta1[3] = 0.2509463534656 * l.piE;//piE
    co.Delta1[4] = 3.0 * M_PI / 180.0; //for xi [radian]
    co.Delta1[5] = l.tE * 0.254674; // t0 [days] -- same fractional-of-tE step as tE itself,
    // since t0's local curvature scale is set by the same characteristic timescale.
    // Placeholder pending Step C3's step-size convergence sweep, like the other four.
    co.Delta1[6] = 0.05; // mbs0 [mag] -- Rubin baseline magnitude
    co.Delta1[8] = 0.05; // mbs1 [mag] -- Roman baseline magnitude
    // The model magnitude depends on mbs linearly with unit slope, so the finite difference is
    // exact for any step and this value only has to avoid underflow. fb0 (index 2) and fb1
    // (index 7) reuse the telescope-keyed co.bb[] steps set inside the data loop below, which
    // are binned so that fb + step never leaves the physical range [0,1].

    // Three matrices, zeroed together: joint, Rubin-only, Roman-only (see SurveyIdx in Bulge.h).
    for (int q = 0; q < NSURV; ++q) {
        co.nepochA[q] = 0;
        co.okA[q] = 0;
        co.okB[q] = 0;
        // Must be reset here too: a partition rejected for having too few epochs never reaches
        // invert_matrix, so without this it would silently carry the previous event's value.
        co.condA[q] = -1.0;
        co.condB[q] = -1.0;
        for (int j = 0; j < Nx; ++j) {
            for (int k = 0; k < Nx; ++k) {
                gsl_matrix_set(co.inputA[q].get(), j, k, 0.0);
                gsl_matrix_set(co.inverA[q].get(), j, k, 0.0);
            }
        }
    }



    for (int i = 0; i < ndw; ++i) {//data
        tt = int(l.tele[i]);//telescope[0,1] LSST, ELT
        const int surv = surveyOfTele(tt); // which single-survey matrix this epoch also feeds
        co.nepochA[SJOINT] += 1;
        co.nepochA[surv]   += 1;

        // Step size for the fb (blend fraction) derivative must be binned against
        // *this epoch's own telescope's* blend fraction, s.fb[tt] -- not always
        // s.fb[0] (Rubin). fb is a flux ratio, physically bounded to [0,1]; these
        // bin edges (0.15, 0.85) and step sizes (0.07, 0.15) are chosen so that
        // fb[tt] + bb[1] never crosses 1.0 and fb[tt] + bb[0] never crosses 0.0 --
        // but only if bb[] is computed from the same fb it perturbs. Binning off
        // s.fb[0] while perturbing s.fb[tt] (tt=1, Roman) broke that guarantee:
        // Roman's F146 blend fraction is routinely much higher than Rubin's r-band
        // one (see Step B3 -- F146's PSF is ~10x smaller), so a low-fb[0] bin's
        // large "+0.15" step, applied to an already-high fb[1], pushed fb[tt] past
        // 1.0 and tripped CHECK(s.fb[tt] <= 1.0) -- an uncaught std::runtime_error
        // that hard-crashed the whole program. Discovered running Step C1's
        // acceptance test; unrelated to that step, fixed here first because it
        // blocks numerical verification of every Fisher-matrix step.
        if (s.fb[tt] < 0.15)      {co.bb[0] =+ 0.07; co.bb[1] =+ 0.15;}
        else if (s.fb[tt] < 0.85) {co.bb[0] =- 0.07; co.bb[1] =+ 0.07;}
        else                      {co.bb[0] =- 0.07; co.bb[1] =- 0.15;}

        for (int j = 0; j < Nx; ++j) {
            for (int h = 0; h < 2; ++h) {

                if (j == 0) {co.diff = double(+co.Delta1[j] * sig[h]) ;      l.u0 += co.diff;}
                if (j == 1) {co.diff = double(+co.Delta1[j] * sig2[h]);      l.tE += co.diff;}
                if (j == 2) {co.diff = double(co.bb[h]);               s.fb[tt] += co.diff;}
                // Per-telescope flux parameters (Step C2b). Each affects ONLY its own telescope's
                // epochs. On an epoch from the other telescope the model magnitude is unchanged, so
                // the derivative is identically zero and that epoch contributes nothing to this
                // parameter's row -- which is exactly right: Rubin's data says nothing about Roman's
                // blend fraction. co.diff is still set to a nonzero dummy so the division is safe.
                if (j == 6) {co.diff = (tt == 0) ? double(+co.Delta1[j] * sig[h]) : 1.0;
                                 if (tt == 0) s.mbs[0] += co.diff;}
                if (j == 7) {co.diff = (tt == 1) ? double(co.bb[h]) : 1.0;
                                 if (tt == 1) s.fb[1]  += co.diff;}
                if (j == 8) {co.diff = (tt == 1) ? double(+co.Delta1[j] * sig[h]) : 1.0;
                                 if (tt == 1) s.mbs[1] += co.diff;}
                if (j == 3) {co.diff = double(+co.Delta1[j] * sig2[h]);     l.piE += co.diff;}
                if (j == 4) {co.diff = double(+co.Delta1[j] * sig[h]) ;      s.xi += co.diff;}
                if (j == 5) {co.diff = double(+co.Delta1[j] * sig[h]) ;      l.t0 += co.diff;}

                lightcurve(s, l, as, l.timn[i]);
                s.Astar = (s.ut * s.ut + 2.0) / std::sqrt(s.ut * s.ut * (s.ut * s.ut + 4.0));
                co.magw = s.mbs[tt] - 2.5 * std::log10(s.Astar * s.fb[tt] + 1.0 - s.fb[tt]);
                co.derm1[h] = double(co.magw - l.magn[i]) / co.diff;

                CHECK(l.tE > 0.0);
                CHECK(s.fb[tt] > 0.0);
                CHECK(s.fb[tt] <= 1.0);
                CHECK(l.piE > 0.0);
                CHECK(l.u0 != 0.0);
                CHECK(co.diff != 0.0);
                CHECK(l.errm[i] > 0.0);

                if (j == 0)      l.u0 -= co.diff;
                if (j == 1)      l.tE -= co.diff;
                if (j == 2)  s.fb[tt] -= co.diff;
                if (j == 3)     l.piE -= co.diff;
                if (j == 4)      s.xi -= co.diff;
                if (j == 5)      l.t0 -= co.diff;
                if (j == 6 and tt == 0) s.mbs[0] -= co.diff;
                if (j == 7 and tt == 1)  s.fb[1] -= co.diff;
                if (j == 8 and tt == 1) s.mbs[1] -= co.diff;
            }

            co.derm1f = double(co.derm1[0] + co.derm1[1]) * 0.5;

            for (int k = 0; k <= j; ++k) {
                for (int h = 0; h < 2;  ++h) {

                    if (k == 0) {co.diff = double(+co.Delta1[k] * sig[h]) ;     l.u0 += co.diff;}
                    if (k == 1) {co.diff = double(+co.Delta1[k] * sig2[h]);     l.tE += co.diff;}
                    if (k == 2) {co.diff = double(co.bb[h]);                s.fb[tt] += co.diff;}
                    // Per-telescope flux parameters (Step C2b). Each affects ONLY its own telescope's
                    // epochs. On an epoch from the other telescope the model magnitude is unchanged, so
                    // the derivative is identically zero and that epoch contributes nothing to this
                    // parameter's row -- which is exactly right: Rubin's data says nothing about Roman's
                    // blend fraction. co.diff is still set to a nonzero dummy so the division is safe.
                    if (k == 6) {co.diff = (tt == 0) ? double(+co.Delta1[k] * sig[h]) : 1.0;
                                     if (tt == 0) s.mbs[0] += co.diff;}
                    if (k == 7) {co.diff = (tt == 1) ? double(co.bb[h]) : 1.0;
                                     if (tt == 1) s.fb[1]  += co.diff;}
                    if (k == 8) {co.diff = (tt == 1) ? double(+co.Delta1[k] * sig[h]) : 1.0;
                                     if (tt == 1) s.mbs[1] += co.diff;}
                    if (k == 3) {co.diff = double(+co.Delta1[k] * sig2[h]);    l.piE += co.diff;}
                    if (k == 4) {co.diff = double(+co.Delta1[k] * sig[h]) ;     s.xi += co.diff;}
                    if (k == 5) {co.diff = double(+co.Delta1[k] * sig[h]) ;     l.t0 += co.diff;}

                    lightcurve(s, l, as, l.timn[i]);
                    s.Astar = (s.ut * s.ut + 2.0) / std::sqrt(s.ut * s.ut * (s.ut * s.ut + 4.0));
                    co.magw = s.mbs[tt] - 2.5 * std::log10(s.Astar * s.fb[tt] + 1.0 - s.fb[tt]);
                    co.derm2[h] = double(co.magw - l.magn[i]) / co.diff;

                    CHECK(l.tE > 0.0);
                    CHECK(s.fb[tt] >= 0.0);
                    CHECK(s.fb[tt] <= 1.0);
                    CHECK(l.piE > 0.0);
                    CHECK(l.u0 != 0.0);
                    CHECK(co.diff != 0.0);
                    CHECK(l.errm[i] > 0.0);

                    if (k == 0)      l.u0 -= co.diff;
                    if (k == 1)      l.tE -= co.diff;
                    if (k == 2)  s.fb[tt] -= co.diff;
                    if (k == 3)     l.piE -= co.diff;
                    if (k == 4)      s.xi -= co.diff;
                    if (k == 5)      l.t0 -= co.diff;
                    if (k == 6 and tt == 0) s.mbs[0] -= co.diff;
                    if (k == 7 and tt == 1)  s.fb[1] -= co.diff;
                    if (k == 8 and tt == 1) s.mbs[1] -= co.diff;
                }

                co.derm2f = double(co.derm2[0] + co.derm2[1]) * 0.5;

                // ACCUMULATE, do not overwrite. A Fisher matrix is by definition a sum of
                // information over independent measurements:
                //     F_jk = sum_i (dm_i/dtheta_j)(dm_i/dtheta_k) / sigma_i^2
                // Each light-curve epoch i contributes one term. Using gsl_matrix_set with the
                // bare term (as this line previously did) discarded every epoch but the last,
                // leaving a rank-1 matrix: singular for any Nx > 1, so invert_matrix fell into
                // its deter==0 branch, nudged the diagonal by 1e-10, and returned an inverse of
                // order 1e10 -- which is why every reported sigma was astronomically large
                // (relative sigma ~1e7 on tE for real detected events). The pre-loop zeroing of
                // inputA is what makes this running sum well-defined.
                // Partitioned accumulation (Step C5). The derivatives above are identical for
                // all three matrices -- same event, same model. What differs is only which
                // epochs are allowed to contribute. Each epoch feeds the joint matrix and
                // exactly one single-survey matrix, which is why
                //     F[SJOINT] == F[SRUBIN] + F[SROMAN]
                // holds exactly, element by element.
                {
                    const double contrib = co.derm1f * co.derm2f / (l.errm[i] * l.errm[i]);
                    gsl_matrix_set(co.inputA[SJOINT].get(), j, k,
                                   gsl_matrix_get(co.inputA[SJOINT].get(), j, k) + contrib);
                    gsl_matrix_set(co.inputA[surv].get(), j, k,
                                   gsl_matrix_get(co.inputA[surv].get(), j, k) + contrib);
                }
            }
        }//end of for J
    }//end of data for

    for (int q = 0; q < NSURV; ++q) {
        for (int j = 0; j < Nx; ++j) {
            for (int k = (j + 1); k < Nx; ++k) {
                double element = gsl_matrix_get(co.inputA[q].get(), k, j);
                gsl_matrix_set(co.inputA[q].get(), j, k, element);
            }
        }
        // A partition with fewer epochs than free parameters cannot possibly constrain them:
        // the information matrix is rank-deficient by construction. Report that honestly
        // instead of inverting it into a meaningless ~1e10 sigma. This is a physical statement,
        // not a numerical workaround -- a short event peaking in a Roman gap genuinely has no
        // Roman data, and "Roman cannot characterize this event" is the correct answer.
        // Guard against a partition with fewer epochs than the parameters it must constrain.
        // The relevant count is the ACTIVE subset for that survey, not the full Nx.
        co.okA[q] = (co.nepochA[q] >= static_cast<int>(
                         activePhotParams(q, co.nepochA[SRUBIN], co.nepochA[SROMAN]).size()))
                        ? invert_matrix(co, 0, q) : 0;
    }

/*cout<<"Covariance matrix Photometry"<<endl;
   for(int i=0; i<Nx; ++i){ for(int j=0; j<Nx; ++j){cout<<co.inputA[i][j]<<"\t   ";}  cout<<"\n"; }
   cout<<"*************Input*******************"<<endl;
   for(int i=0; i<Nx; ++i){ for(int j=0; j<Nx; ++j){cout<<co.inverA[i][j]<<"\t   ";}  cout<<"\n"; }
   cout<<"************Inverse*****************"<<endl;  */
    //double summ;

    /*for(int i = 0; i < Nx; ++i) {
        for(int j = 0; j < Nx; ++j) {
            summ = 0.0;
            for(int k = 0; k < Nx;  ++k){
                summ += gsl_matrix_get(co.inputA, i, k) * gsl_matrix_get(co.inverA, k, j);
            }
            gsl_matrix_set(co.summ, i, j, summ);
        }
    }*/
    // Source - https://stackoverflow.com/q/40687568
    // Posted by Sara Fuerst, modified by community. See post 'Timeline' for change history
    // Retrieved 2026-02-09, License - CC BY-SA 3.0
/*
    gsl_matrix_set_zero(co.summA);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, co.inverA, co.inputA, 0.0, co.summA);

    gsl_matrix *I = gsl_matrix_alloc(Nx, Nx);
    gsl_matrix *D = gsl_matrix_alloc(Nx, Nx);

    gsl_matrix_memcpy(D, co.summA);
    gsl_matrix_sub(D, I);

    //if ((gsl_matrix_max(D) > 0.2) or (gsl_matrix_min(D) < -0.2)) co.flagi = -1;
    double minA, maxA;
    gsl_matrix_minmax(D, &minA, &maxA);
    if (maxA > 0.2) co.flagi -= 2;
    cout << co.flagi;
    gsl_matrix_free(D);
    gsl_matrix_free(I);
*/
/*
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Nx; ++j) {
            if ((i == j and std::fabs(gsl_matrix_get(co.summA, i, j) - 1.0) > 0.2) or (i != j and std::fabs(gsl_matrix_get(co.summA, i, j) - 0.0) > 0.2)) {
                //cout<<"Error_Inverse sum:  "<<co.summ[i][j]<<"\t i:  "<<i<<"\t j:  "<<j<<endl;
                //cout<<""<<l.t0<<"\t"<<l.u0<<"\t"<<l.tE<<"\t"<<s.xi*M_PI/180.0<<"\t"<<s.fb[tt]<<"\t"<<s.mbs[tt]<<"\t"<<l.piE<<endl;
                co.flagi = -1;
            }
        }
    }
*/

/*
    cout << "Covariance matrix Photometry" << endl;
    print_mat_contents(co.inputA, Nx);
    cout << "*************Input*******************" << endl;
    print_mat_contents(co.inverA, Nx);
    cout << "************Inverse*****************" << endl;
*/
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///                   Astrometry                                            ///
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    co.Delta2[0] = l.tetE * 0.2500467345692345;
    co.Delta2[1] = s.mus1 * 0.2500465923465493;
    co.Delta2[2] = s.mus2 * 0.25000465936443;
    co.Delta2[3] = 0.25000375423543264 * l.piE;
    for (int q = 0; q < NSURV; ++q)
    for(int j = 0; j < Ny; ++j){
        for(int k = 0; k < Ny; ++k){
            gsl_matrix_set(co.inputB[q].get(), j, k, 0.0);
            gsl_matrix_set(co.inverB[q].get(), j, k, 0.0);
            //gsl_matrix_set(co.adjunB, j, k, 0.0);
            //gsl_matrix_set(co.tempoB, j, k, 0.0);
            //gsl_matrix_set(co.summ, j, k, 0.0);
        }
    }

    for (int i = 0; i < ndw; ++i) {
        const int survB = surveyOfTele(int(l.tele[i])); // same partition as the photometric side
        for (int j = 0; j < Ny;  ++j) {

            for (int h = 0; h < 2; ++h) {
                if (j == 0) { co.diff = double(co.Delta2[j] * sig2[h]); l.tetE += co.diff; }
                if (j == 1) { co.diff = double(co.Delta2[j] * sig[h]);  s.mus1 += co.diff; }
                if (j == 2) { co.diff = double(co.Delta2[j] * sig[h]);  s.mus2 += co.diff; }
                if (j == 3) { co.diff = double(co.Delta2[j] * sig2[h]); l.piE  += co.diff; }


                lightcurve(s, l, as, l.timn[i]);
                co.dera1[h] = double(s.pos1c - l.soux[i]) / co.diff;
                co.derb1[h] = double(s.pos2c - l.souy[i]) / co.diff;

                CHECK(l.tetE > 0.0);
                CHECK(l.piE > 0.0);
                CHECK(co.diff != 0.0);
                CHECK(!(s.pos1c == l.soux[i] && s.pos2c == l.souy[i]));
                CHECK(l.erra[i] > 0.0);

                if (j==0) l.tetE -= co.diff;
                if (j==1) s.mus1 -= co.diff;
                if (j==2) s.mus2 -= co.diff;
                if (j==3)  l.piE -= co.diff;
            }

            co.dera1f = (co.dera1[0] + co.dera1[1]) * 0.5;
            co.derb1f = (co.derb1[0] + co.derb1[1]) * 0.5;


            for (int k = 0; k <= j; ++k) {
                for (int h = 0; h < 2; ++h) {
                    if (k==0) {co.diff = double(co.Delta2[k] * sig2[h]);  l.tetE += co.diff;}
                    if (k==1) {co.diff = double(co.Delta2[k] * sig[h] );  s.mus1 += co.diff;}
                    if (k==2) {co.diff = double(co.Delta2[k] * sig[h] );  s.mus2 += co.diff;}
                    if (k==3) {co.diff = double(co.Delta2[k] * sig2[h]);  l.piE  += co.diff;}

                    lightcurve(s, l, as, l.timn[i]);
                    co.dera2[h] = double(s.pos1c - l.soux[i]) / co.diff;
                    co.derb2[h] = double(s.pos2c - l.souy[i]) / co.diff;

                    CHECK(l.tetE > 0.0);
                    CHECK(l.piE > 0.0);
                    CHECK(co.diff != 0.0);
                    CHECK(!(s.pos1c == l.soux[i] && s.pos2c == l.souy[i]));

                    if (k==0) l.tetE -= co.diff;
                    if (k==1) s.mus1 -= co.diff;
                    if (k==2) s.mus2 -= co.diff;
                    if (k==3)  l.piE -= co.diff;
                }

                co.dera2f = (co.dera2[0] + co.dera2[1]) * 0.5;
                co.derb2f = (co.derb2[0] + co.derb2[1]) * 0.5;

                // ACCUMULATE, do not overwrite -- same defect and same reasoning as the
                // photometric matrix above, applied to the astrometric one. Here each epoch
                // contributes both sky coordinates (pos1c/pos2c), hence the two derivative
                // products summed before dividing by the per-epoch astrometric variance.
                {
                    const double contribB = (co.dera1f * co.dera2f + co.derb1f * co.derb2f)
                                              / (l.erra[i] * l.erra[i] * 2.0);
                    gsl_matrix_set(co.inputB[SJOINT].get(), j, k,
                                   gsl_matrix_get(co.inputB[SJOINT].get(), j, k) + contribB);
                    gsl_matrix_set(co.inputB[survB].get(), j, k,
                                   gsl_matrix_get(co.inputB[survB].get(), j, k) + contribB);
                }
            }
        }//end of loop J
    }//end of loop data

    for (int q = 0; q < NSURV; ++q) {
        for(int j = 0; j < Ny; ++j) {
            for(int k = (j+1); k < Ny; ++k) {
                double element = gsl_matrix_get(co.inputB[q].get(), k, j);
                gsl_matrix_set(co.inputB[q].get(), j, k, element);
            }
        }
        co.okB[q] = (co.nepochA[q] >= Ny) ? invert_matrix(co, 1, q) : 0;
    }

    /*for(int i = 0; i < Ny; ++i) {
        for(int j = 0; j < Ny; ++j) {
            summ = 0.0;
            for(int k = 0; k < Ny;  ++k){
                summ += gsl_matrix_get(co.inputB, i, k) * gsl_matrix_get(co.inverB, k, j);
            }
            gsl_matrix_set(co.summ, i, j, summ);
        }
    }*/

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, co.inverB[SJOINT].get(), co.inputB[SJOINT].get(), 0.0, co.summB.get());
/*
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Ny; ++j) {
            if ((i == j and std::fabs(gsl_matrix_get(co.summB, i, j) - 1.0) > 0.2) or (i != j and std::fabs(gsl_matrix_get(co.summB, i, j) - 0.0) > 0.2)) {
                //cout<<"Error_Inverse sum:  "<<co.summ[i][j]<<"\t i:  "<<i<<"\t j:  "<<j<<endl;
                //cout<<""<<l.t0<<"\t"<<l.u0<<"\t"<<l.tE<<"\t"<<s.xi*M_PI/180.0<<"\t"<<s.fb[tt]<<"\t"<<s.mbs[tt]<<"\t"<<l.piE<<endl;
                co.flagi = -1;
            }
        }
    }
*/
/*
    cout << "Covariance matrix Astrometry" << endl;
    print_mat_contents(co.inputB, Ny);
    cout << "*************Input*******************" << endl;
    print_mat_contents(co.inverB, Ny);
    cout << "************Inverse*****************" << endl;
*/
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Error parameters Calculation               //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void ErrorCal(covarian & co, lens & l , source & s){
  double corr1;

  // Per-survey 1-sigma forecasts (Step C5). An invalid partition -- too few epochs, or a
  // singular information matrix -- gets sigma = -1.0 as an explicit "not characterizable"
  // sentinel, rather than a plausible-looking number derived from a nudged singular matrix.
  // Downstream analysis MUST test okA/okB before using these.
  for (int q = 0; q < NSURV; ++q) {
      // -1.0 marks both "partition not characterizable" and "parameter not in this partition's
      // active subset" -- Rubin's matrix says nothing about Roman's flux parameters and vice
      // versa. Downstream code must test okA[] and skip negative sigmas rather than treat them
      // as measurements.
      for (int k = 0; k < Nx; ++k) co.Era[q][k] = -1.0;
      if (co.okA[q]) {
          for (int k : activePhotParams(q, co.nepochA[SRUBIN], co.nepochA[SROMAN])) {
              co.Era[q][k] = std::sqrt(std::fabs(gsl_matrix_get(co.inverA[q].get(), k, k)));
          }
      }
      for (int k = 0; k < Ny; ++k) {
          co.Erb[q][k] = co.okB[q] ? std::sqrt(std::fabs(gsl_matrix_get(co.inverB[q].get(), k, k)))
                                   : -1.0;
      }
  }

  // Everything below is the JOINT event summary, unchanged in meaning: resu[] keeps its existing
  // semantics and index layout so the output columns and per-field aggregation still work.
  // The per-survey numbers live in Era[]/Erb[] and are reported separately.
  const auto& Era = co.Era[SJOINT];
  const auto& Erb = co.Erb[SJOINT];

  corr1 = double(gsl_matrix_get(co.inverA[SJOINT].get(), 1, 4) / std::fabs(Era[1] * Era[4])); //correlation tE  && xi

  co.resu[0]  = double(Era[0] / (std::fabs(l.u0)    + eps));
  co.resu[1]  = double(Era[1] / (std::fabs(l.tE)    + eps));
  co.resu[2]  = double(Era[2] / (std::fabs(s.fb[0]) + eps));
  co.resu[3]  = double(Era[3] / (std::fabs(l.piE)   + eps));
  co.resu[4]  = double(Era[4] / (std::fabs(s.xi)    + eps));
  co.resu[5]  = double(Erb[0] / (std::fabs(l.tetE)  + eps));
  co.resu[6]  = double(Erb[1] / (std::fabs(s.mus1)  + eps));
  co.resu[7]  = double(Erb[2] / (std::fabs(s.mus2)  + eps));
  co.resu[8]  = double(Erb[3] / (std::fabs(l.piE)   + eps));

  co.resu[3]  = MIN(co.resu[3], co.resu[8]);

  co.resu[9]  = std::sqrt(co.resu[3] * co.resu[3] + co.resu[5] * co.resu[5]); //Lens Mass
  co.resu[10] = std::fabs(co.resu[9] * (s.Ds - l.Dl) / s.Ds); //Dl

  co.f1 = co.resu[5] * co.resu[5] + co.resu[1] * co.resu[1] + (Era[4] * std::tan(s.xi)) * (Era[4] * std::tan(s.xi))
    - 2.0 * co.resu[1] * Era[4] * std::tan(s.xi) * corr1;
  co.f2 = co.resu[5] * co.resu[5] + co.resu[1] * co.resu[1] + (Era[4] / std::tan(s.xi)) * (Era[4] / std::tan(s.xi))
    - 2.0 * co.resu[1] * Era[4] / std::tan(s.xi) * corr1;

  co.sigmul1 = std::sqrt(Erb[1] * Erb[1] + l.murel * l.murel * std::cos(s.xi) * std::cos(s.xi) * std::fabs(co.f1));
  co.sigmul2 = std::sqrt(Erb[2] * Erb[2] + l.murel * l.murel * std::sin(s.xi) * std::sin(s.xi) * std::fabs(co.f2));

  co.resu[11] = double(co.sigmul1 / (std::fabs(l.mul1) + eps)); //mu_lens_1
  co.resu[12] = double(co.sigmul2 / (std::fabs(l.mul2) + eps)); //mu_lens_2
  co.resu[13] = std::sqrt(co.sigmul1 * co.sigmul1 + co.sigmul2 * co.sigmul2) / (std::fabs(l.mul) + eps); //mu_lens
  co.resu[14] = std::sqrt(Erb[1]  * Erb[1]  +  Erb[2] * Erb[2])  / (std::fabs(s.mus) + eps); // mu_source
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Light curves and Astrometry                //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void lightcurve(source & s, lens & l, astromet & as, double timh)
{
    double dvex, dvey, Ve_x, tt;
    double pis = double(1.0/s.Ds); //[mas] source parallax
    double pil = double(1.0/l.Dl); //[mas] lens parallax
    double int1 = 0.0, int2 = 0.0;


    for (int ig = 0; ig < 2; ++ig) {
        if(ig == 0) tt = timh;
        if(ig == 1) tt = 0.0;
        
        dvex = +vearth * std::sin(omegae * tt + M_PI / 2.0) / omegae;//km/s
        dvey = -vearth * std::cos(omegae * tt + M_PI / 2.0) / omegae;//km/s
        as.Ve_n1 =  std::cos(tetp) * dvex * std::sin(l.deltao) - dvey * std::cos(l.deltao);
        Ve_x     = -std::cos(tetp) * dvex * std::cos(l.deltao) - dvey * std::sin(l.deltao);
        as.Ve_n2 = -std::sin(s.FI) * Ve_x + std::cos(s.FI) * std::sin(tetp) * dvex;

        if(ig == 0) { int1  = as.Ve_n1; int2  = as.Ve_n2; }
        if(ig == 1) { int1 -= as.Ve_n1; int2 -= as.Ve_n2; }
    }
    as.ue_n1 = int1; //[radian] without dimention
    as.ue_n2 = int2; //[radian] without dimention

    s.ux=-l.u0*std::sin(s.xi) + (timh-l.t0)*std::cos(s.xi)/l.tE + l.piE*as.ue_n1;//[] +parallax
    s.uy= l.u0*std::cos(s.xi) + (timh-l.t0)*std::sin(s.xi)/l.tE + l.piE*as.ue_n2;//[] +parallax

    s.ut0=std::sqrt((s.ux-l.piE*as.ue_n1)*(s.ux-l.piE*as.ue_n1) + (s.uy-l.piE*as.ue_n2)*(s.uy-l.piE*as.ue_n2));
    s.ut= std::sqrt(s.ux*s.ux+ s.uy*s.uy);

    if(s.ut==0.0)   s.ut=1.0e-50;
    if(s.ut0==0.0)  s.ut0=1.0e-50;


    s.def1a = (s.ux - l.piE * as.ue_n1) * l.tetE / (s.ut0 * s.ut0 + 2.0); //x-deflection, without parallax[mas]
    s.def2a = (s.uy - l.piE * as.ue_n2) * l.tetE / (s.ut0 * s.ut0 + 2.0); //y-deflection, without parallax[mas]

    s.def1c = s.ux * l.tetE / (s.ut * s.ut + 2.0); //x-deflection[mas]
    s.def2c = s.uy * l.tetE / (s.ut * s.ut + 2.0); //y-deflection[mas]

    s.pos1b = -l.u0 * l.tetE * std::sin(s.xi) + s.mus1 * (timh - l.t0) - as.ue_n1 * pis; //x-source trajectory+parallax[mas]
    s.pos2b = +l.u0 * l.tetE * std::cos(s.xi) + s.mus2 * (timh - l.t0) - as.ue_n2 * pis; //y-source trajectory+parallax[mas]

    s.pos1c = -l.u0 * l.tetE * std::sin(s.xi) + s.mus1 * (timh - l.t0) - as.ue_n1 * pis + s.def1c; //x-source trajectory+parallax[mas]+lensing
    s.pos2c = +l.u0 * l.tetE * std::cos(s.xi) + s.mus2 * (timh - l.t0) - as.ue_n2 * pis + s.def2c; //y-source trajectory+parallax[mas]+lensing

    l.pos1  = l.mul1 * (timh - l.t0) - as.ue_n1 * pil ;//x-lens trajectory && parallax[mas]
    l.pos2  = l.mul2 * (timh - l.t0) - as.ue_n2 * pil ;//y-lens trajectory && parallax[mas]

   //if(s.ut<double(5.0*s.ros) or s.ut0<double(5.0*s.ros)){
   //cout<<"Error finite ros: "<<s.ros<<"\t ut:  "<<s.ut<<"\t ut0:  "<<s.ut0<<endl;
   //int uue; cin>>uue;
  //}


}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Optical Depth calculations                 //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void optical_depth(source& s)
{
    double ds = (double)s.nums * step; //kpc
    double CC = 4.0 * G * M_PI * ds * ds * std::pow(10.0, 9.0) * Msun / (velocity * velocity * KP);
    double dl, x, dx;

    s.od_thin = s.od_thick = s.od_bulge = s.od_halo = s.opt = 0.0;
    
    for (int k = 1; k < s.nums; ++k) {
        dl = (double)k * step;///kpc
        x  = dl / ds;
        dx = (double)step / ds / 1.0;
        s.od_thin  += s.rho_thin[k]  * x * (1.0 - x) * dx * CC;
        s.od_thick += s.rho_thick[k] * x * (1.0 - x) * dx * CC;
        s.od_bulge += s.rho_bulge[k] * x * (1.0 - x) * dx * CC;
        s.od_halo  += s.rho_halo[k]  * x * (1.0 - x) * dx * CC;
    }
    s.opt = std::fabs(s.od_thin + s.od_thick + s.od_bulge + s.od_halo);
//    cout << "total_opticalD: " << s.opt << "\t od_thin: " << s.od_thin << endl;
//    cout << "od_thick: " << s.od_thick << "\t od_bulge: " << s.od_bulge << "\t od_halo: " << s.od_halo << endl;
}

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Galactic model calculations                //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void Disk_model(source& s, int numt)
{
    double x, xb, yb, zb, r4, r2, rdi, Rb;
    double nnf = 0.4 / 0.8;
    double mBarre; //stars/pc^3
    double Rx0, Ry0, Rz0, Rc, cp, cn, rhoE, rhoS;
    double alfa = 12.89 / RAa;
    double xf, yf, zf, rho;

    s.Romaxs = s.Nstart = s.Rostart = 0.0;
    s.Romins = 10000000000.0;

    double fd    = 1.0;  //see the program mass_averaged.cpp.  we do not apply any limitation
    double fb    = 1.0;  //just stars brighter than V=11.5, but we change to consider all stars
    double fh    = 1.0;  //No limitation
    double Rdd   = 2.17; //2.53; ///2.17;
    double Rhh   = 1.33; //1.32; //1.33;

//    double frac = 0.05; //fraction of halo in the form of compact objects
//    frac was replaced with 1.0 and some rescaling was applyed manually

    char filename[40];
    FILE *filj;
    FILE *fill;

    int flagf = 0;
    if (numt < 10) {
       flagf = 1;
       filj  = fopen("./files/density/lb.txt", "a+");
       sprintf(filename, "./files/density/%c%d.dat", 'D', numt);
       fill = fopen(filename, "w");
    }

    for (int i = 1; i < Num; ++i) {
        s.rho_thin[i]  = 0;
        s.rho_bulge[i] = 0;
        s.rho_halo[i]  = 0;
        s.rho_thick[i] = 0;
        x  = double(i  * step); //kpc
        zb = std::sin(s.FI) * x;
        yb = std::cos(s.FI) * std::sin(s.TET) * x;
        xb = x * std::cos(s.FI) * std::cos(s.TET) - Dsun;
        Rb = std::sqrt(xb * xb + yb * yb);

///========== Galactic Thin Disk =====================
        for (int ii = 0; ii < 8; ++ii) {
            rdi = Rb * Rb + zb * zb / (epci[ii] * epci[ii]);
            if (ii == 0) {
                rho = std::exp(-rdi / 25.0) - std::exp(-rdi / 9.0);
            }
            else if (ii > 0) {
                rho = std::exp(-std::sqrt(0.25 + rdi / (Rdd * Rdd))) - std::exp(-std::sqrt(0.25 + rdi / (Rhh * Rhh)));
            }
            rho *= 1.2; // totalmass= 4.25e10
            s.rho_thin[i] = std::fabs(s.rho_thin[i] + rho0[ii] * corr[ii] * 0.001 * rho/d0[ii]);
        } //Msun/pc^3

///========== Galactic Thick Disk =====================
        double rho00 = 1.34 * 0.001 + 3.04 * 0.0001;
        if (std::fabs(zb) < 0.4) {
            s.rho_thick[i] = std::fabs((rho00 / 0.999719) * std::exp(-(Rb - Dsun) / 2.5) * (1.0 - zb * zb / (0.4 * 0.8 * (2.0 + nnf))));
        }
        else {
            s.rho_thick[i] = std::fabs((rho00 / 0.999719) * std::exp(-(Rb - Dsun) / 2.5) * std::exp(nnf) * std::exp(-std::fabs(zb) / 0.8) / (1.0 + 0.5 * nnf)); //Msun/pc^3
        }
        s.rho_thick[i] *= 2.67; //total_mass=0.8e10

///========== Galactic Stellar Halo=================
        rdi = std::sqrt(Rb * Rb + zb * zb / (0.76 * 0.76));
        if (rdi <= 0.5) {
            s.rho_halo[i] = std::fabs(1.0 * (0.932 * 0.00001 / 867.067) * std::pow(0.5 / Dsun, - 2.44));
        }
        else {
            s.rho_halo[i] = std::fabs(1.0 * (0.932 * 0.00001 / 867.067) * std::pow(rdi / Dsun, - 2.44)); //Msun/pc^3
        }
        s.rho_halo[i] *= 5281.0; //Total_mass=1.2e9

///========== Galactic bulge =====================
        constexpr double barMassRescale = 0.24529; // calibrated so bulge column density toward
                                                   // Baade's Window (l=1, b=-3.9) matches the
                                                   // Han & Gould (2003) HST benchmark: 2086 Msun/pc^2
        xf =  xb * std::cos(alfa) + yb * std::sin(alfa);
        yf = -xb * std::sin(alfa) + yb * std::cos(alfa);
        zf =  zb;

        Rx0    = 1.46;
        Ry0    = 0.49;
        Rz0    = 0.39;
        Rc     = 3.43;
        cp     = 3.007;
        cn     = 3.329;
        mBarre = 35.45 / 3.84723 * barMassRescale;

        r4  = std::pow(std::pow(std::fabs(xf / Rx0), cn)
                     + std::pow(std::fabs(yf / Ry0), cn), cp / cn)
            + std::pow(std::fabs(zf / Rz0), cp);
        r4  = std::pow(std::fabs(r4), 1.0 / cp);
        r2  = std::sqrt(std::fabs(xf * xf + yf * yf));

        if (r2 <= Rc) {
            rhoS = mBarre * 1.0 / (std::cosh(-r4) * std::cosh(-r4));
        }
        else {
            rhoS = mBarre * 1.0 / (std::cosh(-r4) * std::cosh(-r4)) * std::exp(-4.0 * (r2 - Rc) * (r2 - Rc));
        }

        Rx0    = 4.44;
        Ry0    = 1.31;
        Rz0    = 0.80;
        Rc     = 6.83;
        cp     = 2.786;
        cn     = 3.917;
        mBarre = 2.27 / 87.0 * barMassRescale; // normalized

        r4  = std::pow(std::fabs(std::pow(std::fabs(xf / Rx0), cn)
                               + std::pow(std::fabs(yf / Ry0), cn)), cp / cn)
            + std::pow(std::fabs(zf / Rz0), cp);
        r4  = std::pow(r4, 1.0 / cp);
        r2  = std::sqrt(std::fabs(xf * xf + yf * yf));

        if (r2 <= Rc) {
            rhoE = mBarre * std::exp(-r4);
        }
        else {
            rhoE = mBarre * std::exp(-r4) * std::exp(-4.0 * (r2 - Rc) * (r2 - Rc));
        }

        s.rho_bulge[i]  = std::fabs(rhoS) + std::fabs(rhoE); ///Msun/pc^3
        s.rho_bulge[i] *= 0.45; //total mass= 1.7e10
///==================================================================

        s.Rostar0[i] = std::fabs(s.rho_thin[i] + s.rho_thick[i] + s.rho_bulge[i] + s.rho_halo[i]); //[Msun/pc^3]
        s.Rostari[i] = s.Rostar0[i] * x * x * step * 1.0e9 * (M_PI / 180.0) * (M_PI / 180.0); //[Msun/deg^2]
        s.Nstari[i]  = binary_fraction
                     * (s.rho_thin[i] * fd / 0.403445
                      + s.rho_thick[i] * fh / 0.4542 
                      + s.rho_halo[i] * fh / 0.4542 
                      + s.rho_bulge[i] * fb / 0.308571); //[Nt/pc^3]
        s.Nstari[i]  = s.Nstari[i] * x * x * step * 1.0e9 * (M_PI / 180.0) * (M_PI / 180.0); //[Ni/deg^2]

        s.Nstart  += s.Nstari[i];  //[Nt/deg^2]
        s.Rostart += s.Rostari[i]; //[Mt/deg^2]

        if (s.Rostari[i] > s.Romaxs) {
            s.Romaxs = s.Rostari[i]; //source selection
        }
        if (s.Rostari[i] < s.Romins) {
            s.Romins = s.Rostari[i]; //source selection
        }


        if (flagf > 0) {
            fprintf(fill, "%e   %e   %e   %e   %e  %e  %e\n",
                x, s.rho_thin[i], s.rho_bulge[i], s.rho_thick[i], s.rho_halo[i], s.Rostar0[i], s.Nstari[i]);
        }
    }
    if (flagf > 0) {
        fprintf(filj, "%.5lf  %.5lf  %.5lf  %.5lf\n", s.lon, s.lat, s.TET * RAa, s.FI * RAa);
        fclose(fill);
        fclose(filj);
    }

    cout << "Nstart [Nt/deg^2]: "    << s.Nstart << "\t Ro_star [Mass/deg^2]: " << s.Rostart << endl;
    cout << "Romaxs:  "              << s.Romaxs << "\t Romins:  "              << s.Romins  << endl;
    cout << ">>>>>>>>>>>>>>>>>>>>>> END OF DISK MODLE <<<<<<<<<<<<<<<<<<<<" << endl;
    //exit(0);
}

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

////////////////////// Matrix
// Largest condition number we will still trust an inverse from.
//
// A double carries ~16 significant decimal digits; inverting a matrix with condition number 10^k
// costs roughly k of them. The precision gains this project sets out to measure are at the
// few-percent level, so we need at least 3-4 trustworthy digits in the result, leaving room for
// about 10^12. Past that the reported sigma is dominated by round-off rather than by data.
//
// This threshold is applied to the NORMALIZED matrix, where a large value can no longer be
// blamed on parameters being measured in different units and instead means the data genuinely
// cannot constrain some combination of them. The per-event condition number is stored regardless,
// so a stricter (or looser) cut can be applied during analysis without rerunning the simulation.
constexpr double kMaxCondition = 1.0e12;

int invert_matrix(covarian & co, int flag, int surv)
{
    // Returns 1 on a usable inverse, 0 if the information matrix cannot be trusted -- singular,
    // or too ill-conditioned. A rejected partition means "this data cannot characterize this
    // event", a physical result that must be reported as such rather than smoothed over.
    //
    // For the photometric matrix only a subset of parameters is inverted, since a single-survey
    // partition carries no information about the other telescope's flux parameters (see
    // activePhotParams). The reduced inverse is scattered back into the full-size matrix, leaving
    // inactive rows and columns at zero; ErrorCal reports sigma = -1 for those.
    const bool photometric = (flag == 0);
    gsl_matrix*  in   = photometric ? co.inputA[surv].get() : co.inputB[surv].get();
    gsl_matrix*  out  = photometric ? co.inverA[surv].get() : co.inverB[surv].get();
    double&      cond = photometric ? co.condA[surv] : co.condB[surv];

    static const std::vector<int> kAllAst = {0, 1, 2, 3};
    const std::vector<int> actPhot = photometric
        ? activePhotParams(surv, co.nepochA[SRUBIN], co.nepochA[SROMAN])
        : std::vector<int>();
    const std::vector<int>& act = photometric ? actPhot : kAllAst;
    const int dim = static_cast<int>(act.size());

    cond = -1.0;
    gsl_matrix_set_zero(out);

    // ---- 1. Normalizing scale D = diag(1/sqrt(F_ii)) over the active parameters ----
    //
    // F_jk carries units of 1/(theta_j theta_k), and with tE in days (~30), u0 dimensionless
    // (~0.3), xi in radians and mbs in magnitudes, the entries span many orders of magnitude
    // before any physics enters. That ill-conditioning is an artifact of our choice of units and
    // is removable. A non-positive diagonal means the parameter has no information at all.
    std::vector<double> scale(dim);
    for (int i = 0; i < dim; ++i) {
        const double d = gsl_matrix_get(in, act[i], act[i]);
        if (!std::isfinite(d) || d <= 0.0) return 0;
        scale[i] = 1.0 / std::sqrt(d);
    }

    // ---- 2. Reduced, normalized matrix Ftilde = D F D, unit diagonal ----
    gsl_matrix *Ft = gsl_matrix_alloc(dim, dim);
    if (!Ft) throw std::runtime_error("GSL matrix allocation failed");
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            gsl_matrix_set(Ft, i, j,
                           gsl_matrix_get(in, act[i], act[j]) * scale[i] * scale[j]);

    // ---- 3. Condition number of the normalized matrix ----
    //
    // Symmetric positive semi-definite, so the 2-norm condition number is lambda_max/lambda_min.
    // gsl_eigen_symm destroys its input, hence the copy.
    {
        gsl_matrix *ev = gsl_matrix_alloc(dim, dim);
        gsl_vector *ew = gsl_vector_alloc(dim);
        gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(dim);
        gsl_matrix_memcpy(ev, Ft);
        gsl_eigen_symm(ev, ew, w);

        double lmin = gsl_vector_get(ew, 0), lmax = lmin;
        for (int i = 1; i < dim; ++i) {
            const double v = gsl_vector_get(ew, i);
            if (v < lmin) lmin = v;
            if (v > lmax) lmax = v;
        }
        gsl_eigen_symm_free(w);
        gsl_vector_free(ew);
        gsl_matrix_free(ev);

        if (!std::isfinite(lmin) || !std::isfinite(lmax) || lmin <= 0.0) {
            gsl_matrix_free(Ft);
            return 0;   // some parameter combination carries no information at all
        }
        cond = lmax / lmin;
        if (cond > kMaxCondition) {
            gsl_matrix_free(Ft);
            return 0;   // genuinely degenerate: report not-characterizable
        }
    }

    // ---- 4. Invert the reduced normalized matrix ----
    gsl_matrix      *lu = gsl_matrix_alloc(dim, dim);
    gsl_matrix      *Fi = gsl_matrix_alloc(dim, dim);
    gsl_permutation *p  = gsl_permutation_alloc(dim);
    gsl_matrix_memcpy(lu, Ft);

    int s;
    gsl_linalg_LU_decomp(lu, p, &s);
    const double det = gsl_linalg_LU_det(lu, s);
    co.deter = det;

    if (det == 0.0 || !std::isfinite(det)) {
        gsl_permutation_free(p); gsl_matrix_free(Fi);
        gsl_matrix_free(lu); gsl_matrix_free(Ft);
        return 0;
    }
    gsl_linalg_LU_invert(lu, p, Fi);

    // ---- 5. Undo the scaling and scatter back to full size: F^-1 = D Ftilde^-1 D ----
    //
    // Exact, not approximate: (D F D)^-1 = D^-1 F^-1 D^-1. The whole manoeuvre is algebraically
    // a no-op; its purpose is that the matrix handed to LU has unit diagonal and a far smaller
    // condition number.
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            gsl_matrix_set(out, act[i], act[j],
                           gsl_matrix_get(Fi, i, j) * scale[i] * scale[j]);

    gsl_permutation_free(p);
    gsl_matrix_free(Fi);
    gsl_matrix_free(lu);
    gsl_matrix_free(Ft);

    // A valid covariance matrix has non-negative variances on the diagonal.
    for (int i = 0; i < dim; ++i) {
        const double v = gsl_matrix_get(out, act[i], act[i]);
        if (!std::isfinite(v) || v < 0.0) {
            gsl_matrix_set_zero(out);
            return 0;
        }
    }
    return 1;
}

void print_mat_contents(gsl_matrix *matrix, int size)
{
    int i, j;
    double element;

    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            element = gsl_matrix_get(matrix, i, j);
            printf("%f ", element);
        }
        printf("\n");
    }
}

