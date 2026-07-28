#include "Bulge.h"

///============================================================================

time_t _timeNow;
unsigned int _randVal;
unsigned int _dummyVal;
FILE * _randStream;

///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//

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
 
    // --------------------- Read extinction ------------------------
    readBayestar(*ex,"./files/ext/");
    std::cout << "**** File extinctionf.txt was read ****\n";

    // --------------------- Call read_cmd --------------------------
    read_cmd(*cm);
    std::cout << "******* read_cmd was done ************" << std::endl;

    int    save = 0, flagm, flagL, gg = -1; //ss, qq, ww, vv, zz, pp, counter, 
    int    nri = -1, nde = -1, icon;
    int    nlens;// hh; // nde1, nri1, 
    int    gi,       ndw, sq, ndd;
    int    flag_det; // nml = 0;
    int    flagf,  fi; // datf1, datf2;
    double errs,   errg, fdet, minc, cade; // fel, , mind
    double magnio, test, deltaA; // dist,  
    double Astar0, As1,  As0;
    double initial;
    double trajm, trajp; //  ddf;
    double chi1,  chi2,   chi3, chi1a, chi2a, chi3a, sil, sil2;
    double dchiL, dchiP,  dchiA;
    double flag0, flag1,  flag2;
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
        for (s->lat = -1.4; s->lat <= -1.3; s->lat += dd) {
//        for (s->lat = -3.9; s->lat <= -3.8; s->lat += dd) {
            if (s->lon < lx and s->lat > bx) continue;
            nde += 1;
            cout << ">>>>>>>>>>> NEW STEP " << nde << " <<<<<<<<\t nri:  " << nri << endl;
            cout << "longtitude: " << s->lon << "\t latitude: " << s->lat << endl;

                ndd = 0; minc = 100000.0; cade = 0.0;
                
                for (int i = 0; i < 1000; ++i) ls->ct[i] = -1; // Initialize ct
                
                for (int i = 0; i < Nl; ++i) {
                    if (std::sqrt((s->lon - ls->l[i]) * (s->lon - ls->l[i]) + (s->lat - ls->b[i]) * (s->lat - ls->b[i])) <= FoV) {
                        // If LSST sees the source
                        ls->ct[ndd] = i;
                        if (ndd > 0) {
                            cade = ls->tim[ls->ct[ndd]] - ls->tim[ls->ct[ndd]-1];

                            if (minc > cade) minc = cade;

                            if (ndd >= 999) break; // Break out of the loop when ndd reached 999
                            CHECK(cade > 0.0);
                            CHECK(ls->tim[i] >= 0.0);
                            CHECK(ls->tim[i] <= Tobs);
                            CHECK(ls->filter[i] >= 0);
                            CHECK(ls->filter[i] < 6);
                            CHECK(ls->sig5[i] >= 0.0);
                            CHECK(ls->sig5[i] <= 40.0);
                            CHECK(minc > 0.0);
                        }
                        ndd += 1;
                    }
                }
                cout << "ndd: " << ndd << "\t minc: " << minc << endl;
 
                icon  = 0;
                nlens = 0;
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
//                    std::cerr << "nsim=" << nsim << "  Ds=" << s->Ds << "  mass=" << s->mass
//                              << "  nums=" << s->nums << "  Ml=" << l->Ml << "  u0=" << l->u0 << "\n";
                    optical_depth(*s);

                    s->nssim[s->nums] += 1.0;
                    flagf   = 0;
                    flagm   = 0;
                    initial = 0.0;
                    test = RandR(0.0, 1.0);

                    for (int i = 0; i < coun; ++i) {
                        l->timn[i] = 0.0;  l->magn[i] = 0.0; l->soux[i] = 0.0;  l->souy[i] = 0.0;
                        l->errm[i] = 0.0;  l->erra[i] = 0.0; l->tele[i] = -1;
                    }

                    ndw     = 0;   flag_det = 0;
                    flag0   = 0.0; flag1    = 0.0; flag2 = 0.0;
                    chi1    = 0.0; chi2     = 0.0; chi3  = 0.0;
                    chi1a   = 0.0; chi2a    = 0.0; chi3a = 0.0;
                    def1p   = 0.0; s->def1c = 0.0; vsave = 0.0;
                    def2p   = 0.0; s->def2c = 0.0;
                    s->errM = 0.0; s->errA  = 0.0;
                    dchiL   = 0.0; dchiP    = 0.0; dchiA = 0.0;

                    dt   = 60.0;///days
                    fdet = 0.0;

                    for (int i = 0; i < M; ++i) {
                        Mpeak = s->magb[i] - 2.5 * std::log10(l->A0 * s->blend[i] + 1.0 - s->blend[i]);
                            cout << "i=" << i << "  Mab=" << s->Mab[i] << "  Map=" << s->Map[i]
                                 << "  blend=" << s->blend[i] << "  Mpeak=" << Mpeak << endl;
                        if (Mpeak <= thre[i] and s->magb[i] > satu[i])    fdet += 1.0;
                    }

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
                    if (test <= s->blend[2] and fdet > 1.0) { //at least detectable in two filters
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
                            sq = int(ls->ct[gi]);
    
                            if (tim >= 0.0 and tim <= Tobs and tim >= ls->tim[int(ls->ct[0])] and tim <= ls->tim[int(ls->ct[ndd - 1])] and
                                gi < ndd and sq >= 0 and sq <= static_cast<int>(Nl) and tim >= ls->tim[sq]) {
    
                                fi = int(ls->filter[sq]);
    
                                if (magni[fi] >= satu[fi] and magni[fi] <= thre[fi]) {
                                    errg = errlsstM(magni[fi], int(fi), double(ls->sig5[sq])); //[mag]
                                    errs = errlsstA(*ls, magni[2]); ///[mas]
    
                                    deltaA = std::fabs(std::pow(10.0, -0.4 * errg) - 1.0) * (s->blend[fi] * s->Astar + 1.0 - s->blend[fi]);
                                    magnio = magni[fi] + RandN(errg, 3.0);
    
                                    chi1 += std::fabs((magnio -   magni[fi]) * (magnio -   magni[fi]) / (errg * errg)); //real
                                    chi2 += std::fabs((magnio -  magni0[fi]) * (magnio -  magni0[fi]) / (errg * errg)); //real no parallax
                                    chi3 += std::fabs((magnio - s->magb[fi]) * (magnio - s->magb[fi]) / (errg * errg)); //baseline
    
                                    sil  = RandN(errs * std::sqrt(2.0), 3.0);
                                    sil2 = RandN(errs, 3.0);
    
                                    trajm = std::sqrt(s->pos1b * s->pos1b + s->pos2b * s->pos2b); //stright + parallax
                                    trajp = std::sqrt(s->pos1c * s->pos1c + s->pos2c * s->pos2c); //stright + parallax+lensing
    
                                    chi1a += std::fabs((trajp + sil - trajp) * (trajp + sil - trajp) / (errs * errs * 2.0)); //real trajectory
                                    chi2a += std::fabs((trajp + sil - trajm) * (trajp + sil - trajm) / (errs * errs * 2.0)); //real without lensing
                                    chi3a += std::fabs( sil2  * sil2 / (errs * errs));
     
                                    l->timn[ndw] = tim;
                                    l->magn[ndw] = magni[2]; //scale to r-LSST band, but the errors are different
                                    l->errm[ndw] = errg;
                                    l->soux[ndw] = s->pos1c;
                                    l->souy[ndw] = s->pos2c;
                                    l->erra[ndw] = errs;
                                    l->tele[ndw] = 0;
    
                                    flag2 = 0.0;
                                    if (std::fabs(magnio - s->magb[fi]) > std::fabs(3.0 * errg))    flag2 = 1.0;
                                    if (ndw > 2 and float(flag0 + flag1 + flag2) > 2.0)   flag_det = 1;
                                    flag0 = flag1;
                                    flag1 = flag2;
    
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
                                    CHECK(ndw <= ndd);
                                    CHECK(errs >= 0.0);
                                    CHECK(errg >= 0.0);
                                    CHECK(deltaA >= 0.0);
                                    CHECK(gi <= ndd);
    
                                    s->errM += deltaA;
                                    s->errA += errs;
                                    vsave += std::sqrt(vs1 * vs1 + vs2 * vs2);
    
                                    ndw += 1;
                                }//magnitude limit
                                gi += 1;
                            }
    
                            if ( tim < -50.0 or tim > (10.0 * year + 50.0)) dt = 60.0; //days
    
                            else if (tim < ls->tim[int(ls->ct[0])] or tim > ls->tim[int(ls->ct[ndd-1])]) dt=3.0; //days
    
                            else {
                                if (gi > 0 and gi < ndd) cade = float(ls->tim[int(ls->ct[gi])] - ls->tim[int(ls->ct[gi-1])]); //days
                                else  cade = minc;
                                dt = double(cade);
                            }
    
                        }//end of loop time
                        if (flagm > 0) {
                            fil4.close();
                            fil5.close();
                        }
                    }// end of visible star

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH 
                    for (int i = 0; i < nq; ++i) co->resu[i] = -1.0;
    
                    FFG[0] = 0;
                    FFG[1] = 0;
                    FFG[2] = 0;
    
                    //cout << "flagf: " << flagf << "ndw: " << ndw << endl;
    
                    if (flagf == 0 or ndw <= 2) {
                        errg    = errlsstM(s->magb[2], 2, double(24.43)); //r-band
                        s->errA = errlsstA(*ls, s->magb[2]); //r-band
                        s->errM = std::fabs(std::pow(10.0, - 0.4 * errg) - 1.0); //r-band
                        dchiL = 0.0;
                        dchiP = 0.0;
                        dchiA = 0.0;
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
       
                        if (s->FWHM < Tobs and dchiL > float(2.0 * ndw) and flag_det > 0 and ndw > 10) { //lensing
                            FFG[0] = 1; //Lensing
                            nlens += FFG[0];
                            FisherM(*s, *l, *as, *co, ndw);
    
                            if (co->flagi > 0) {
                                nerr += 1.0;
                                ErrorCal(*co, *l, *s);
    
                                std::ofstream fil0_append(fnLDt, std::ios::app);
                                fil0_append << std::fixed << std::setprecision(5)
                                            << l->Ml       << " " << l->Dl       << " " << l->mul       << " " << s->fb[0]     << " " << s->mbs[0]   << " "
          << std::setprecision(7)           << l->tE       << " " << l->murel    << " " << l->u0        << " " << s->lon       << " " << s->lat      << " "
                                            << l->piE      << " " << l->tetE     << " " << co->resu[1]  << " " << co->resu[2]  << " " << co->resu[3] << " "
                                            << co->resu[5] << " " << co->resu[9] << " " << co->resu[10] << " " << co->resu[13] << "\n";
                                fil0_append.close();
                            }
                        }
                        gg = FunctE(*l);
                        //ss = int(FuncMl(*l));
                        //qq = int(FuncPi(*l));
                        //ww = int(Funcu0(*l));
                        //vv = int(FuncMu(*l));
                        //zz = int(FuncMb(*l, s->Map[2]));
                        //pp = int(FuncFb(*l, s->blend[2]));
                    }
//                CHECK(gg >= 0);

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
                    co->flagi, s->Ai[2]
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
                        << s->Map[2]      << " " << s->nsbl[2]     << " " << co->flagi                  << " " << s->Ai[2]        << "\n";
                filg_in.close();
//            }


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

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Fisher and Covariance matrixes             //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void FisherM(source & s, lens & l, astromet & as,  covarian & co, int ndw)
{
    co.flagi =+ 1;
    int tt; // tev;

    if (s.fb[0] < 0.15)      {co.bb[0] =+ 0.07; co.bb[1] =+ 0.15;}
    else if (s.fb[0] < 0.85) {co.bb[0] =- 0.07; co.bb[1] =+ 0.07;}
    else                     {co.bb[0] =- 0.07; co.bb[1] =- 0.15;}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///        Photometry                                  ///
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    co.Delta1[0] = 0.1507586576;//u0
    co.Delta1[1] = l.tE * 0.254674;//tE [days]
    co.Delta1[3] = 0.2509463534656 * l.piE;//piE
    co.Delta1[4] = 3.0 * M_PI / 180.0; //for xi [radian]

    for (int j = 0; j < Nx; ++j) {
        for (int k = 0; k < Nx; ++k) {
            gsl_matrix_set(co.inputA.get(), j, k, 0.0);
            gsl_matrix_set(co.inverA.get(), j, k, 0.0);
        }
    }



    for (int i = 0; i < ndw; ++i) {//data
        tt = int(l.tele[i]);//telescope[0,1] LSST, ELT

        for (int j = 0; j < Nx; ++j) {
            for (int h = 0; h < 2; ++h) {

                if (j == 0) {co.diff = double(+co.Delta1[j] * sig[h]) ;      l.u0 += co.diff;}
                if (j == 1) {co.diff = double(+co.Delta1[j] * sig2[h]);      l.tE += co.diff;}
                if (j == 2) {co.diff = double(co.bb[h]);               s.fb[tt] += co.diff;}
                if (j == 3) {co.diff = double(+co.Delta1[j] * sig2[h]);     l.piE += co.diff;}
                if (j == 4) {co.diff = double(+co.Delta1[j] * sig[h]) ;      s.xi += co.diff;}

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
            }

            co.derm1f = double(co.derm1[0] + co.derm1[1]) * 0.5;

            for (int k = 0; k <= j; ++k) {
                for (int h = 0; h < 2;  ++h) {

                    if (k == 0) {co.diff = double(+co.Delta1[k] * sig[h]) ;     l.u0 += co.diff;}
                    if (k == 1) {co.diff = double(+co.Delta1[k] * sig2[h]);     l.tE += co.diff;}
                    if (k == 2) {co.diff = double(co.bb[h]);                s.fb[tt] += co.diff;}
                    if (k == 3) {co.diff = double(+co.Delta1[k] * sig2[h]);    l.piE += co.diff;}
                    if (k == 4) {co.diff = double(+co.Delta1[k] * sig[h]) ;     s.xi += co.diff;}

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
                }

                co.derm2f = double(co.derm2[0] + co.derm2[1]) * 0.5;
                gsl_matrix_set(co.inputA.get(), j, k, co.derm1f * co.derm2f / (l.errm[i] * l.errm[i]));
            }
        }//end of for J
    }//end of data for

    for (int j = 0; j < Nx; ++j) {
        for (int k = (j + 1); k < Nx; ++k) {
            double element = gsl_matrix_get(co.inputA.get(), k, j);
            gsl_matrix_set(co.inputA.get(), j, k, element);
        }
    }
    //gsl_matrix_transpose_memcpy()
    invert_matrix(co, 0);

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
    for(int j = 0; j < Ny; ++j){
        for(int k = 0; k < Ny; ++k){
            gsl_matrix_set(co.inputB.get(), j, k, 0.0);
            gsl_matrix_set(co.inverB.get(), j, k, 0.0);
            //gsl_matrix_set(co.adjunB, j, k, 0.0);
            //gsl_matrix_set(co.tempoB, j, k, 0.0);
            //gsl_matrix_set(co.summ, j, k, 0.0);
        }
    }

    for (int i = 0; i < ndw; ++i) {
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
                gsl_matrix_set(co.inputB.get(), j, k, (co.dera1f * co.dera2f + co.derb1f * co.derb2f) / (l.erra[i] * l.erra[i] * 2.0));
            }
        }//end of loop J
    }//end of loop data

    for(int j = 0; j < Ny; ++j) {
        for(int k = (j+1); k < Ny; ++k) {
            double element = gsl_matrix_get(co.inputB.get(), k, j);
            gsl_matrix_set(co.inputB.get(), j, k, element);
        }
    }
    invert_matrix(co, 1);

    /*for(int i = 0; i < Ny; ++i) {
        for(int j = 0; j < Ny; ++j) {
            summ = 0.0;
            for(int k = 0; k < Ny;  ++k){
                summ += gsl_matrix_get(co.inputB, i, k) * gsl_matrix_get(co.inverB, k, j);
            }
            gsl_matrix_set(co.summ, i, j, summ);
        }
    }*/

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, co.inverB.get(), co.inputB.get(), 0.0, co.summB.get());
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

  for(int k=0; k < Nx; ++k)  co.Era[k] = std::sqrt(std::fabs(gsl_matrix_get(co.inverA.get(), k, k)));
  for(int k=0; k < Ny; ++k)  co.Erb[k] = std::sqrt(std::fabs(gsl_matrix_get(co.inverB.get(), k, k)));

  corr1 = double(gsl_matrix_get(co.inverA.get(), 1, 4) / std::fabs(co.Era[1] * co.Era[4])); //correlation tE  && xi

  co.resu[0]  = double(co.Era[0] / (std::fabs(l.u0)    + eps));
  co.resu[1]  = double(co.Era[1] / (std::fabs(l.tE)    + eps));
  co.resu[2]  = double(co.Era[2] / (std::fabs(s.fb[0]) + eps));
  co.resu[3]  = double(co.Era[3] / (std::fabs(l.piE)   + eps));
  co.resu[4]  = double(co.Era[4] / (std::fabs(s.xi)    + eps));
  co.resu[5]  = double(co.Erb[0] / (std::fabs(l.tetE)  + eps));
  co.resu[6]  = double(co.Erb[1] / (std::fabs(s.mus1)  + eps));
  co.resu[7]  = double(co.Erb[2] / (std::fabs(s.mus2)  + eps));
  co.resu[8]  = double(co.Erb[3] / (std::fabs(l.piE)   + eps));

  co.resu[3]  = MIN(co.resu[3], co.resu[8]);

  co.resu[9]  = std::sqrt(co.resu[3] * co.resu[3] + co.resu[5] * co.resu[5]); //Lens Mass
  co.resu[10] = std::fabs(co.resu[9] * (s.Ds - l.Dl) / s.Ds); //Dl

  co.f1 = co.resu[5] * co.resu[5] + co.resu[1] * co.resu[1] + (co.Era[4] * std::tan(s.xi)) * (co.Era[4] * std::tan(s.xi))
    - 2.0 * co.resu[1] * co.Era[4] * std::tan(s.xi) * corr1;
  co.f2 = co.resu[5] * co.resu[5] + co.resu[1] * co.resu[1] + (co.Era[4] / std::tan(s.xi)) * (co.Era[4] / std::tan(s.xi))
    - 2.0 * co.resu[1] * co.Era[4] / std::tan(s.xi) * corr1;

  co.sigmul1 = std::sqrt(co.Erb[1] * co.Erb[1] + l.murel * l.murel * std::cos(s.xi) * std::cos(s.xi) * std::fabs(co.f1));
  co.sigmul2 = std::sqrt(co.Erb[2] * co.Erb[2] + l.murel * l.murel * std::sin(s.xi) * std::sin(s.xi) * std::fabs(co.f2));

  co.resu[11] = double(co.sigmul1 / (std::fabs(l.mul1) + eps)); //mu_lens_1
  co.resu[12] = double(co.sigmul2 / (std::fabs(l.mul2) + eps)); //mu_lens_2
  co.resu[13] = std::sqrt(co.sigmul1 * co.sigmul1 + co.sigmul2 * co.sigmul2) / (std::fabs(l.mul) + eps); //mu_lens
  co.resu[14] = std::sqrt(co.Erb[1]  * co.Erb[1]  +  co.Erb[2] * co.Erb[2])  / (std::fabs(s.mus) + eps); // mu_source
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
    double CC = 4.0 * G * M_PI * ds * ds * std::pow(10.0,9.0) * Msun / (velocity * velocity * KP);
    double dl, x, dx;

    s.od_disk = s.od_ThD = s.od_bulge = s.od_halo = s.opt = 0.0;
    
    for (int k = 1; k < s.nums; ++k) {
        dl = (double)k * step;///kpc
        x  = dl / ds;
        dx = (double)step / ds / 1.0;
        s.od_disk  += s.rho_thin[k]  * x * (1.0 - x) * dx * CC;
        s.od_ThD   += s.rho_thick[k]   * x * (1.0 - x) * dx * CC;
        s.od_bulge += s.rho_bulge[k] * x * (1.0 - x) * dx * CC;
        s.od_halo  += s.rho_halo[k]  * x * (1.0 - x) * dx * CC;
    }
    s.opt = std::fabs(s.od_disk + s.od_ThD + s.od_bulge + s.od_halo);// + s.od_dlmc + s.od_blmc + s.od_hlmc);///total
    //cout<<"total_opticalD: "<<s.opt<<"\t od_disk: "<<s.od_disk<<endl;
    //cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
    //cout<<"od_dlmc:  "<<s.od_dlmc<<"\t od_blmc: "<<s.od_blmc<<"\t od_hlmc:  "<<s.od_hlmc<<endl;
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

    double frac = 0.05; //fraction of halo in the form of compact objects

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
            s.rho_thin[i] = std::fabs(s.rho_thin[i] + rho0[ii] * corr[ii] * 0.001 * rho/d0[ii]);
        } //Msun/pc^3

///========== Galactic Thick Disk =====================
        double rho00 = 1.34 * 0.001 + 3.04 * 0.0001;
        if (std::fabs(zb) < 0.4) {
            s.rho_thick[i] = std::fabs((rho00 / 0.999719) * std::exp(-(Rb - Dsun) / 2.5) * (1.0 - zb * zb / (0.4 * 0.8 * (2.0 + nnf))));
        }
        else {
            s.rho_thick[i] = std::fabs((rho00 / 0.999719) * std::exp(-(Rb - Dsun) / 2.5) * std::exp(nnf) * std::exp(-std::fabs(zb) / 0.8)/(1.0 + 0.5 * nnf)); //Msun/pc^3
        }

///========== Galactic Stellar Halo=================
        rdi = std::sqrt(Rb * Rb + zb * zb / (0.76 * 0.76));
        if (rdi <= 0.5) {
            s.rho_halo[i] = std::fabs(frac * (0.932 * 0.00001 / 867.067) * std::pow(0.5 / Dsun, - 2.44));
        }
        else {
            s.rho_halo[i] = std::fabs(frac * (0.932 * 0.00001 / 867.067) * std::pow(rdi / Dsun, - 2.44)); //Msun/pc^3
        }

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
        mBarre = 2.27 / 87.0 * barMassRescale; // needs to be normalized

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

        s.rho_bulge[i] = std::fabs(rhoS) + std::fabs(rhoE); ///Msun/pc^3
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
void invert_matrix(covarian & co, int flag)
{
    if (flag == 0){//Nx
        gsl_permutation *p = gsl_permutation_alloc(Nx);
        int s;

        gsl_linalg_LU_decomp(co.inputA.get(), p, &s);
        co.deter = gsl_linalg_LU_det(co.inputA.get(), s);

        if (co.deter == 0.0) {
            for (int i = 0; i < Nx; i++) {
                double val = gsl_matrix_get(co.inputA.get(), i, i);
                gsl_matrix_set(co.inputA.get(), i, i, val + 1e-10);
            }
        }

        gsl_linalg_LU_invert(co.inputA.get(), p, co.inverA.get());

        gsl_permutation_free(p);
    }

    if (flag == 1){//Ny
        gsl_permutation *p = gsl_permutation_alloc(Ny);
        int s;

        gsl_linalg_LU_decomp(co.inputB.get(), p, &s);
        co.deter = gsl_linalg_LU_det(co.inputB.get(), s);

        if (co.deter == 0.0) {
            for (int i = 0; i < Ny; i++) {
                double val = gsl_matrix_get(co.inputB.get(), i, i);
                gsl_matrix_set(co.inputB.get(), i, i, val + 1e-10);
            }
        }

        gsl_linalg_LU_invert(co.inputB.get(), p, co.inverB.get());

        gsl_permutation_free(p);
    }
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

