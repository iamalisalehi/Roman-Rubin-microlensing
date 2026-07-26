#include "Bulge.h"

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Func source calculations                   //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void func_source(source& s, CMD& cm, const extin& ex, int sightlineIdx) {
    int nums, num;
    double rho, rf;
    double Ds, Av = -1, Alv;
    double maxnb = 0.0;
    std::array<double, M> Map{}, Mab{}, Ai{}; 
    GalacticComponent struc;

    for (int i = 0; i < M; ++i) {
        s.nsbl[i]  = std::fabs(s.Nstart * std::pow(FWHM[i] * 0.5, 2) * M_PI / (3600.0 * 3600.0));
        s.nsbl[i] += RandN(std::sqrt(s.nsbl[i]), 2.0);

        if (s.nsbl[i] <= 1.0)  {
            s.nsbl[i] = 1.0;
        }
        if (s.nsbl[i] > maxnb) {
            maxnb = s.nsbl[i];
        }
        //cout << "nsbl[i]: " << s.nsbl[i] << "\t maxnb: " << maxnb << endl;
    }

    for (int k = 1; k <= int(maxnb + 0.000000034756346); ++k) {
        do {
            nums = int(RandR(5.0, Num - 2.0));
            rho  = RandR(s.Romins, s.Romaxs);
            Ds   = double(nums * step);
        } while (rho > s.Rostari[nums] or Ds < 0.0 or Ds > MaxD); //distance larger than 20.0

//        cout<<"k:  "<<k<<"\t Ds:  "<<Ds<<"\t nums:  "<<nums<<endl;

        rf = RandR(0.0, s.Rostar0[nums]);
        if      (rf <=  s.rho_thin[nums]) {
            struc = GalacticComponent::THIN_DISK;
        }
        else if (rf <= (s.rho_thin[nums] + s.rho_bulge[nums])) {
            struc = GalacticComponent::BULGE;
        }
        else if (rf <= (s.rho_thin[nums] + s.rho_bulge[nums] + s.rho_thick[nums])) {
            struc = GalacticComponent::THICK_DISK;
        }
        else if (rf <= (s.rho_thin[nums] + s.rho_bulge[nums] + s.rho_thick[nums] + s.rho_halo[nums])) {
            struc = GalacticComponent::HALO;
        }
        else {
           std::cerr << "Selected Galactic Component can not be initialized with rf =" << rf << ".\n";
           std::exit(EXIT_FAILURE);
       }  // cout<<"struc :  "<<struc<<endl;

        if (struc == GalacticComponent::THIN_DISK) { //thin disk
            num = int(RandR(0.0, N1 - 1.0));
            for (int i = 0; i < M; ++i) {
                Mab[i] = cm.Mab_thin[num][i];
            }
            if (k == 1) {
                s.mass = cm.mass_thin[num];
                s.age  = cm.age_thin[num];
                s.logT = cm.logT_thin[num];
                s.cl   = cm.cl_thin[num];
                s.typ  = cm.typ_thin[num];
            }
        }

        if (struc == GalacticComponent::BULGE) {// bulge
            num = int(RandR(0.0, N2 - 1.0));
            for (int i = 0; i < M; ++i) {
                Mab[i] = cm.Mab_bulge[num][i];
            }
            if (k == 1) {
                s.mass = cm.mass_bulge[num];
                s.age  = cm.age_bulge[num];
                s.logT = cm.logT_bulge[num];
                s.cl   = cm.cl_bulge[num];
                s.typ  = cm.typ_bulge[num];
            }
        }

        if (struc == GalacticComponent::THICK_DISK) { //thick disk
            num = int(RandR(0.0, N3 - 1.0));
            for (int i = 0; i < M; ++i) {
                Mab[i] = cm.Mab_thick[num][i];
            }
            if (k == 1) {
                s.mass = cm.mass_thick[num];
                s.age  = cm.age_thick[num];
                s.logT = cm.logT_thick[num];
                s.cl   = cm.cl_thick[num];
                s.typ  = cm.typ_thick[num];
            }
        }

        if (struc == GalacticComponent::HALO) {// stellar halo
            num = int(RandR(0.0, N4 - 1.0));
            for (int i = 0; i < M; ++i) {
                Mab[i] = cm.Mab_halo[num][i];
            }
            if (k == 1) {
                s.mass = cm.mass_halo[num];
                s.age  = cm.age_halo[num];
                s.logT = cm.logT_halo[num];
                s.cl   = cm.cl_halo[num];
                s.typ  = cm.typ_halo[num];
            }
        }

        Av = interpExtinctionAlongSightline(ex, sightlineIdx, Ds);
//        cout << "Av: " << Av << endl;
        for (int i = 0; i < M; ++i) {
            Alv = AlAv(lambda_um[i], Rv[static_cast<int>(struc)]); // A_lambda / A_V
            Ai[i] = Av * Alv + RandN(sigma[i], 1.0); //extinction in other bands

            if (Ai[i] < 0.0) {
                Ai[i] = 0.0;
            }
//        cout << "Alv: " << Alv << "\t filter: " << i << "\tRv[struc]: " << Rv[static_cast<int>(struc)]
//             << "\tstruc: " << static_cast<int>(struc) << endl;
//        cout << "Av: " << Av << "\t Ai[i]: " << Ai[i] << endl;
            Map[i] = Mab[i] + 5.0 * std::log10(Ds * 100.0) + Ai[i];

            if(s.nsbl[i] >= k) {
                s.Fluxb[i] += std::pow(10.0, -0.4 * Map[i]);
            }
//        cout << "filter:  " << i << "\tAlv: " << Alv << "\tExt: " << Ai[i]
//             << "\tMab: " << Mab[i] << "\tMap: " << Map[i] << endl;
        }
//        cout << "*****************************" << endl;

        if (k == 1) {
//            cout << "mass:  " << s.mass << "\t age:  " << s.age << "\t type:  " << s.typ  << endl;
//            cout << "cl:  "   << s.cl   << "\t age:  " << s.age << "\t tef:  "  << s.logT << endl;
            s.struc = struc;
            s.Ds = Ds;
            s.nums = nums;
            s.Av = Av;

            for (int i = 0; i < (M); ++i) {
                s.Ai[i]  = Ai[i];
                s.Map[i] = Map[i];
                s.Mab[i] = Mab[i];
            }
        }
        //cout<<"=========================="<<endl;
    } //loop over the stars

    for (int i = 0; i < M; ++i) {
        s.magb[i]  = -2.5 * std::log10(s.Fluxb[i]);
        s.blend[i] = std::pow(10.0, -0.4 * s.Map[i]) / s.Fluxb[i];

        CHECK(s.Fluxb[i] > 0.0);
        CHECK(s.nsbl[i]  >= 1.0);
        CHECK(s.blend[i] <= 1.00001);
        CHECK(s.blend[i] > 0.0);
//        cout << "i = " << i << "\tnsbl[i]: " << s.nsbl[i] << "\tblend[i]: " << s.blend[i] << endl;
//        if (s.nsbl[i] == 1) {
//            CHECK(s.blend[i] >= 1.0);
//        }
        CHECK(std::fabs(s.Map[i] - s.Mab[i] - s.Ai[i] - 5.0 * std::log10(s.Ds * 100.0)) <= 0.1);
        CHECK(s.Ds > 0.0);
        CHECK(Av >= 0.0);
    }

    // For purposes of Fisher Matrix calculations
    s.fb[0]  = s.blend[2]; // r-LSST
    s.fb[1]  = s.blend[6]; // F146
    s.mbs[0] = s.magb[2];  // r-LSST
    s.mbs[1] = s.magb[6];  // F146
//    cout << "Ds:   " << s.Ds << "\t nums: " << s.nums << endl;
//    debug note: enabling the line above showed that func_source runs many many times,
//    but no star is detectable by criteria stated in main function.
//cout<<">>>>>>>>>>>>  End of Func_source <<<<<<<<<<<<<<<<<<<<<"<<endl;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Func lens  calculations                    //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void func_lens(lens & l, source & s){

    double f, test, tt, Am, DD;
    double mmin = Ml_min;
    double mmax = Ml_max;
    l.rhomaxl = 0.0;

    for (int k = 1; k < int(s.nums - 1) ;++k) {
        l.Dl = double(k * step);
        tt = std::sqrt((s.Ds - l.Dl) * l.Dl / s.Ds) * s.Rostar0[k];
        if (tt > l.rhomaxl)  l.rhomaxl = tt;
    }

    do {
        l.numl = int(RandR(1.0, s.nums - 1.0));
        test   =     RandR(0.0, l.rhomaxl);
        l.Dl   = double(l.numl * step);
        tt     = std::sqrt((s.Ds - l.Dl) * l.Dl / s.Ds) * s.Rostar0[l.numl];

        CHECK(tt <= l.rhomaxl);
        CHECK(l.Dl > 0.0);
        CHECK(s.Ds > 0.0);
        CHECK(l.Dl < s.Ds);

    } while(test > tt);
    l.Dl = double(l.numl * step);

    double   randflag = RandR(0.0, s.Rostar0[l.numl]);
    if      (randflag <= (s.rho_thin[l.numl])) {
        l.struc = GalacticComponent::THIN_DISK;
    }
    else if (randflag <= (s.rho_thin[l.numl] + s.rho_bulge[l.numl])) {
        l.struc = GalacticComponent::BULGE;
    }
    else if (randflag <= (s.rho_thin[l.numl] + s.rho_bulge[l.numl] + s.rho_thick[l.numl])) {
        l.struc = GalacticComponent::THICK_DISK;
    }
    else if (randflag <= (s.rho_thin[l.numl] + s.rho_bulge[l.numl] + s.rho_thick[l.numl] + s.rho_halo[l.numl])) {
        l.struc = GalacticComponent::HALO;
    }
    else {
        throw std::runtime_error("Invalid randflag: " + std::to_string(randflag) + ", rho_star0: " + std::to_string(s.Rostar0[l.numl]));
    }
     //cout<<"Dl:  "<<l.Dl<<"\t struc_lens :  "<<l.struc<<endl;

    if (IMnum == 1) {
        l.Ml = RandR(mmin , mmax);
    } //[in solar mass]

    if (IMnum == 2) {
        l.Ml = -1.0;
        do {
            l.Ml = RandR(mmin, mmax); //[3,5000]in solar mass]
            test = RandR(std::pow(mmax, -0.5) * 100.0, std::pow(mmin, -0.5) * 100.0); //[1.4, 58]
            f = std::pow(l.Ml, -0.5) * 100.0;

            CHECK(f >= 1.4);
            CHECK(f <= 58.0);
            CHECK(test >= 1.4);
            CHECK(test <= 58.0);
            CHECK(l.Ml >= 3.0);
            CHECK(l.Ml <= 5000.0);
       
        } while(test > f);
    }

    if (IMnum == 3) {
        l.Ml = -1.0;
        do {
           l.Ml = RandR(mmin, mmax); //[3,5000]in solar mass]
           test = RandR(std::pow(mmax, -1.0) * 1000.0, std::pow(mmin, -1.0) * 1000.0); //[0.2, 334]
           f = std::pow(l.Ml, -1.0) * 1000.0;
           CHECK(f >= 0.19);
           CHECK(f <= 334.0);
           CHECK(test >= 0.19);
           CHECK(test <= 334.0);
           CHECK(l.Ml >= 3.0);
           CHECK(l.Ml <= 5000.0);
        
        } while(test > f);
    }

    if (IMnum == 4) {
       l.Ml = -1.0;
       do {
           l.Ml = RandR(mmin, mmax);//[3,5000]in solar mass]
           test = RandR(std::pow(mmax, -2.0) * 10000.0, std::pow(mmin, -2.0) * 10000.0); //[0.0004, 1112]
           f = std::pow(l.Ml, -2.0) * 10000.0;
         
           CHECK(f >= 0.00039);
           CHECK(f <= 1112.0);
           CHECK(test >= 0.00039);
           CHECK(test <= 1112.0);
           CHECK(l.Ml >= 3.0);
           CHECK(l.Ml <= 5000.0);

       } while(test > f);
    }

    l.xls    = l.Dl / s.Ds;
    l.RE     = std::sqrt(4.0 * G * l.Ml * Msun * s.Ds * KP) / velocity;
    l.RE     = l.RE * std::sqrt(l.xls * (1.0 - l.xls)); //meter
    l.tetE   = l.RE / AU / l.Dl; //[mas]
    s.ros    = 1.0 * Rsun * l.xls / l.RE;
    l.pirel  = 1.0 / l.Dl - 1.0 / s.Ds; //[mas]
    l.piE    = l.pirel / l.tetE; //[]
    l.u0     = RandR(0.001, u0m);
    l.t0     = RandR(2.0, Tobs - 2.0);
    l.DeltaT = std::sqrt(4.0 + l.u0 * l.u0) * l.tetE; //[mas]

    vrel(s,l);
    l.tE  = std::fabs(l.RE / (l.Vt * 1000.0 * 3600.0 * 24.0));//[days]
    l.A0  = (l.u0 * l.u0 + 2.0) / std::sqrt(l.u0 * l.u0 * (l.u0 * l.u0 + 4.0));
    l.mi1 = s.magb[2] - 2.5 * std::log10(std::fabs(l.A0 + 1.0) * 0.5 * s.blend[2] + 1.0 - s.blend[2]);
    l.mi2 = s.magb[2] - 2.5 * std::log10(std::fabs(l.A0 - 1.0) * 0.5 * s.blend[2] + 1.0 - s.blend[2]);
    Am = l.A0 + 1.0;
    DD = double(4.0 - Am * Am + Am * std::sqrt(Am * Am - 4.0)) / (Am * Am * 0.5 - 2.0);
    s.FWHM = 2.0 * l.tE * std::sqrt(std::fabs(DD - l.u0 * l.u0));//days
    // if(DD<double(l.u0*l.u0))  s.FWHM=l.tE;

    //cout<<"Ml:  "<<l.Ml<<"\t RE:  "<<l.RE/AU<<"\t tE:  "<<l.tE<<endl;
    //cout<<"ros:  "<<s.ros<<"\t murel:  "<<l.murel<<"\t tetE:  "<<l.tetE<<endl;
    //cout<<"u0:  "<<l.u0<<"\t pirel:  "<<l.pirel<<"\t piE:  "<<l.piE<<endl;
    //cout<<"mus1:  "<<s.mus1<<"\t mus2:  "<<s.mus2<<"\t mul1:  "<<l.mul1<<"\t mul2:  "<<l.mul2<<endl;
    //cout<<"DD:  "<<DD<<"\t FWHM:  "<<s.FWHM<<endl;

    CHECK(l.tE > 0.0);
    CHECK(l.Dl <= s.Ds);
    CHECK(l.Vt > 0.0);
    CHECK(l.t0 <= Tobs);
    CHECK(l.t0 >= 0.0);
    CHECK(l.Ml <= mmax);
    CHECK(l.Ml >= mmin);
    CHECK(l.pirel > 0.0);
    CHECK(l.RE >= 0.0);
    CHECK(l.xls < 1.0);
    CHECK(l.Dl >= 0.0);
    CHECK(s.Ds >= 0.0);
    CHECK(s.Ds <= MaxD);
    CHECK(l.u0 >= 0.0);
    CHECK(DD >= double(l.u0 * l.u0));
    CHECK(s.FWHM >= 0.0);
    CHECK(Am >= 2.0);

//    cout<<">>>>>>>>>>>> End of func_lens <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
}

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Relative Velocity calculations             //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//

void vrel(source & s, lens & l){
    if (l.Dl == 0.0) l.Dl = 0.00034735;
    double Rlc = std::sqrt(l.Dl * l.Dl * std::cos(s.FI) * std::cos(s.FI) + Dsun * Dsun - 2. * Dsun * l.Dl * std::cos(s.TET) * std::cos(s.FI));
    double Rsc = std::sqrt(s.Ds * s.Ds * std::cos(s.FI) * std::cos(s.FI) + Dsun * Dsun - 2. * Dsun * s.Ds * std::cos(s.TET) * std::cos(s.FI));
    if (Rlc == 0.0) Rlc = 0.0000000000034346123;
    if (Rsc == 0.0) Rsc = 0.000000000004762654134;

    double LVx, SVx, vt2;
    double SVT, SVR, SVZ, LVT, LVR, LVZ;
//    double fv, testfv, age;
    double test, age = -1;
    double  VSunx, vls2, vls1;
    double  tetd ;

    double NN = 3.5;
    double sigma_R_Disk,          sigma_T_Disk,          sigma_Z_Disk;
    double sigma_R_DiskL,         sigma_T_DiskL,         sigma_Z_DiskL;
    double sigma_R_DiskS,         sigma_T_DiskS,         sigma_Z_DiskS;
    double sigma_R_TDisk = 67.0,  sigma_T_TDisk = 51.0,  sigma_Z_TDisk = 42.0;
    double sigma_R_halo  = 131.0, sigma_T_halo  = 106.0, sigma_Z_halo  = 85.0;
    double sigma_R_Bulge = 113.0, sigma_T_Bulge = 115.0, sigma_Z_Bulge = 100.0;

    double Rho[4] = {0.0};
    double maxr = 0.0;
    for (int i = 0; i < 4; ++i) {
        Rho[i] = rho0[i] * corr[i] / d0[i];
        maxr += Rho[i];
    }


    for (int i = 0; i < 2; ++i) {
        test = RandR(0.0, maxr); ///total ages

        if      (test <=  Rho[0])                       {sigma_R_Disk = 16.7; sigma_T_Disk = 10.8; sigma_Z_Disk = 6.0;  age = 0.075;}
        else if (test <= (Rho[0]+Rho[1]))               {sigma_R_Disk = 19.8; sigma_T_Disk = 12.8; sigma_Z_Disk = 8.0;  age = 0.575;}
        else if (test <= (Rho[0]+Rho[1]+Rho[2]))        {sigma_R_Disk = 27.2; sigma_T_Disk = 17.6; sigma_Z_Disk = 10.0; age = 1.5;  }
        else if (test <= (Rho[0]+Rho[1]+Rho[2]+Rho[3])) {sigma_R_Disk = 30.2; sigma_T_Disk = 19.5; sigma_Z_Disk = 13.2; age = 2.5;  }

        if (i == 0) {
            sigma_R_DiskS = sigma_R_Disk;
            sigma_T_DiskS = sigma_T_Disk;
            sigma_Z_DiskS = sigma_Z_Disk;
        }
        if (i == 1) {
            sigma_R_DiskL = sigma_R_Disk;
            sigma_T_DiskL = sigma_T_Disk;
            sigma_Z_DiskL = sigma_Z_Disk;
        }
    }
    CHECK(age >= 0);

//    double v_R_lmc   = -57.0;
//    double v_T_lmc   = -226.0;
//    double v_Z_lmc   =  221.0;
//    double sigma_LMC =  20.2;
//    double err_rlmc  =  13.0; ///error of global velocity
//    double err_tlmc  =  15.0;
//    double err_zlmc  =  19.0;


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if (s.struc == GalacticComponent::THIN_DISK) {///Galactic disk
        SVR = RandN(sigma_R_DiskS, NN);
        SVT = RandN(sigma_T_DiskS, NN);
        SVZ = RandN(sigma_Z_DiskS, NN);
    }

    else if (s.struc == GalacticComponent::BULGE) {///Galactic bulge
        SVR = RandN(sigma_R_Bulge, NN);
        SVT = RandN(sigma_T_Bulge, NN);
        SVZ = RandN(sigma_Z_Bulge, NN);
    }

    else if (s.struc == GalacticComponent::THICK_DISK) {///thick disk
        SVR = RandN(sigma_R_TDisk, NN);
        SVT = RandN(sigma_T_TDisk, NN);
        SVZ = RandN(sigma_Z_TDisk, NN);
    }

    else if (s.struc == GalacticComponent::HALO){///stellar halo
        SVR = RandN(sigma_R_halo, NN);
        SVT = RandN(sigma_T_halo, NN);
        SVZ = RandN(sigma_Z_halo, NN);
    }

    else {
        std::cerr << "Selected Galactic Component " << static_cast<int>(s.struc) << " does not exist.\n";
        std::exit(EXIT_FAILURE);
    }

    if (s.struc == GalacticComponent::THIN_DISK or s.struc == GalacticComponent::THICK_DISK) {
        SVT = SVT + vro_sun * (1.00762 * std::pow(Rsc / Dsun, 0.0394) + 0.00712);
    }

    s.vs = std::sqrt(SVR * SVR + SVT * SVT + SVZ * SVZ);
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if (l.struc == GalacticComponent::THIN_DISK) {///Galactic disk
        LVR = RandN(sigma_R_DiskL, NN);
        LVT = RandN(sigma_T_DiskL, NN);
        LVZ = RandN(sigma_Z_DiskL, NN);
    }

    else if (l.struc == GalacticComponent::BULGE) {///Galactic bulge
        LVR = RandN(sigma_R_Bulge, NN);
        LVT = RandN(sigma_T_Bulge, NN);
        LVZ = RandN(sigma_Z_Bulge, NN);
    }

    else if (l.struc == GalacticComponent::THICK_DISK){///thick disk
        LVR = RandN(sigma_R_TDisk, NN);
        LVT = RandN(sigma_T_TDisk, NN);
        LVZ = RandN(sigma_Z_TDisk, NN);
    }

    else if (l.struc == GalacticComponent::HALO) {///stellar halo
        LVR = RandN(sigma_R_halo, NN);
        LVT = RandN(sigma_T_halo, NN);
        LVZ = RandN(sigma_Z_halo, NN);
    }

    else {
        std::cerr << "Selected Galactic Component" << static_cast<int>(s.struc) << "does not exist.\n";
        std::exit(EXIT_FAILURE);
    }

    if (l.struc == GalacticComponent::THIN_DISK or l.struc == GalacticComponent::THICK_DISK) {
        LVT = LVT + vro_sun * (1.00762 * std::pow(Rlc / Dsun, 0.0394) + 0.00712);
    }

    l.vl = std::sqrt(LVT * LVT + LVZ * LVZ + LVR * LVR);

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  BETA  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    l.betal = 0.0; l.betas = 0.0;
    tetd = s.TET;
    test = double(l.Dl * std::cos(s.FI) * std::sin(tetd) / Rlc);

    double tol = 1e-11;

    if (std::fabs(test - 1.0) < tol)      l.betal = pi/2.0;
    else if (std::fabs(test + 1.0) < tol)  l.betal = -pi/2.0;
    else                          l.betal = std::asin(test);

    test = double(s.Ds * std::cos(s.FI) * std::sin(tetd) / Rsc);

    if (std::fabs(test - 1.0) < 0.01)      l.betas =  pi/2.0;
    else if (std::fabs(test + 1.0) < 0.01) l.betas = -pi/2.0;
    else                              l.betas = std::asin(test);

    if (Dsun < std::fabs(l.Dl * std::cos(s.FI) * std::cos(tetd)) and l.betal >= 0.0) l.betal =   pi - l.betal;
    if (Dsun < std::fabs(s.Ds * std::cos(s.FI) * std::cos(tetd)) and l.betas >= 0.0) l.betas =   pi - l.betas;
    if (Dsun < std::fabs(l.Dl * std::cos(s.FI) * std::cos(tetd)) and l.betal < 0.0)  l.betal = -(pi - std::fabs(l.betal));
    if (Dsun < std::fabs(s.Ds * std::cos(s.FI) * std::cos(tetd)) and l.betas < 0.0)  l.betas = -(pi - std::fabs(l.betas));


    if (!(std::fabs(l.Dl * std::cos(s.FI) * std::sin(tetd)) <= Rlc)) {
        std::cerr << "Check failed: |Dl*cos(FI)*sin(tetd)| <= Rlc\n"
                  << "Dl: " << l.Dl << "\t Ds: " << s.Ds << '\n'
                  << "FI: " << s.FI << "\t TET: " << tetd << "\t betal: " << l.betal << '\n'
                  << "Rlc: " << Rlc << "\t Rsc: " << Rsc << "\t betas: " << l.betas << '\n'
                  << "value: " << std::fabs(l.Dl * std::cos(s.FI) * std::sin(tetd)) << '\n';
        throw std::runtime_error("Lens position exceeds Rlc");
    }
    
    if (!(std::fabs(test) <= 1.0)) {
        std::cerr << "Check failed: |test| <= 1.0\n"
                  << "test: " << test << '\n'
                  << "Dl: " << l.Dl << "\t Ds: " << s.Ds << '\n'
                  << "FI: " << s.FI << "\t TET: " << tetd << "\t betas: " << l.betas << '\n';
        throw std::runtime_error("Invalid test value");
    }
//HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    if (s.TET > pi)  tetd = s.TET - 2.0 * pi;

    l.deltal = pi - std::fabs(tetd) - std::fabs(l.betal);
    l.deltas = pi - std::fabs(tetd) - std::fabs(l.betas);

    if (l.betal < 0.0) l.deltal = -1.0 * l.deltal;
    if (l.betas < 0.0) l.deltas = -1.0 * l.deltas;

    l.deltao = pi - std::fabs(tetd);
    if (tetd < 0.0) l.deltao = -1.0 * l.deltao;

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    s.SV_n1 =+ SVR * std::sin(l.deltas) - SVT * std::cos(l.deltas);
    s.LV_n1 =+ LVR * std::sin(l.deltal) - LVT * std::cos(l.deltal);
    s.VSun_n1 =+ VSunR * std::sin(l.deltao) - VSunT * std::cos(l.deltao);

    SVx = -SVR * std::cos(l.deltas) - SVT * std::sin(l.deltas);
    LVx = -LVR * std::cos(l.deltal) - LVT * std::sin(l.deltal);
    VSunx= -VSunR * std::cos(l.deltao) - VSunT * std::sin(l.deltao);

    s.SV_n2 = -std::sin(s.FI) * (SVx) + std::cos(s.FI) * SVZ;
    s.LV_n2 = -std::sin(s.FI) * (LVx) + std::cos(s.FI) * LVZ;
    s.VSun_n2 = -std::sin(s.FI) * (VSunx) + std::cos(s.FI) * (VSunZ);


    vls1 = l.xls * s.SV_n1 - s.LV_n1 + (1.0 - l.xls) * s.VSun_n1;  ///Source - lens
    vls2 = l.xls * s.SV_n2 - s.LV_n2 + (1.0 - l.xls) * s.VSun_n2;  /// Source -lens
    l.Vt = std::sqrt(std::fabs(vls1 * vls1 + vls2 * vls2 ));


    s.mus1 = double(s.SV_n1 - s.VSun_n1) * 1000.0 * 3600.0 * 24.0 / (s.Ds * AU);//[mas/days]
    s.mus2 = double(s.SV_n2 - s.VSun_n2) * 1000.0 * 3600.0 * 24.0 / (s.Ds * AU);//[mas/days]
    l.mul1 = double(s.LV_n1 - s.VSun_n1) * 1000.0 * 3600.0 * 24.0 / (l.Dl * AU);//[mas/days]
    l.mul2 = double(s.LV_n2 - s.VSun_n2) * 1000.0 * 3600.0 * 24.0 / (l.Dl * AU);//[mas/days]
    s.mus  = std::sqrt(s.mus1 * s.mus1 + s.mus2 * s.mus2);//constant
    l.mul  = std::sqrt(l.mul1 * l.mul1 + l.mul2 * l.mul2);

    l.murel = std::sqrt((s.mus1 - l.mul1) * (s.mus1 - l.mul1) + (s.mus2 - l.mul2) * (s.mus2 - l.mul2));//[mas/days]
    vt2 = l.murel * l.Dl * AU / (1000.0 * 24.0 * 3600.0);//[km/s]



    if (vls1 > 0.0 and vls2 >= 0.0)          s.xi  = std::atan(std::fabs(vls2) / std::fabs(vls1));//in radian
    else if (vls1 <= 0.0 and vls2 >  0.0)    s.xi  = std::atan(std::fabs(vls1) / std::fabs(vls2)) + M_PI / 2.0;
    else if (vls1 <  0.0 and vls2 <= 0.0)    s.xi  = std::atan(std::fabs(vls2) / std::fabs(vls1)) + M_PI;
    else if (vls1 >= 0.0 and vls2 <  0.0)    s.xi  = std::atan(std::fabs(vls1) / std::fabs(vls2)) + 3.0 * M_PI / 2.0;
    else if (s.xi >  2.0 * M_PI)             s.xi -= 2.0 * M_PI;

    CHECK(!(vls1 == 0.0 && vls2 == 0.0 && std::fabs(vt2 - l.Vt) > 0.1));

    CHECK(l.Vt >= 0.0);
    CHECK(l.Vt <= 1.0e6);

   //cout<<"(vrel_func):   l.deltao:  "<<l.deltao<<"FI: "<<s.FI<<endl;
}
