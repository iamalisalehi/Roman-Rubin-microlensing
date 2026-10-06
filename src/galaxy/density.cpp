// Density along the line of sight (Disk_model) and optical depth.
#include "galaxy/density.h"
#include "util/random.h"

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
    double nnf = THICK_NNF;
    double mBarre; //stars/pc^3
    double Rx0, Ry0, Rz0, Rc, cp, cn, rhoE, rhoS;
    double alfa = BAR_ANGLE_DEG / RAa;
    double xf, yf, zf, rho;

    s.Romaxs = s.Nstart = s.Rostart = 0.0;
    s.Romins = 10000000000.0;

    double fd    = DENS_FD;  //see config/parameters.h, section 2a
    double fb    = DENS_FB;
    double fh    = DENS_FH;
    double Rdd   = THIN_RDD;
    double Rhh   = THIN_RHH;

//    double frac = 0.05; //fraction of halo in the form of compact objects
//    frac was replaced with 1.0 and some rescaling was applyed manually

    char filename[40];
    FILE *filj;
    FILE *fill;

    int flagf = 0;
    if (numt < 10) {
       flagf = 1;
       filj  = fopen((std::string(PATH_DENSITY_DIR) + "lb.txt").c_str(), "a+");
       sprintf(filename, "%s%c%d.dat", PATH_DENSITY_DIR, 'D', numt);
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
                rho = std::exp(-rdi / THIN_YOUNG_L1) - std::exp(-rdi / THIN_YOUNG_L2);
            }
            else if (ii > 0) {
                rho = std::exp(-std::sqrt(THIN_CORE + rdi / (Rdd * Rdd))) - std::exp(-std::sqrt(THIN_CORE + rdi / (Rhh * Rhh)));
            }
            rho *= THIN_NORM;
            s.rho_thin[i] = std::fabs(s.rho_thin[i] + rho0[ii] * corr[ii] * 0.001 * rho/d0[ii]);
        } //Msun/pc^3

///========== Galactic Thick Disk =====================
        double rho00 = THICK_RHO00;
        if (std::fabs(zb) < THICK_H1) {
            s.rho_thick[i] = std::fabs((rho00 / THICK_RHO_DIV) * std::exp(-(Rb - Dsun) / THICK_SCALE_LEN) * (1.0 - zb * zb / (THICK_H1 * THICK_H2 * (2.0 + nnf))));
        }
        else {
            s.rho_thick[i] = std::fabs((rho00 / THICK_RHO_DIV) * std::exp(-(Rb - Dsun) / THICK_SCALE_LEN) * std::exp(nnf) * std::exp(-std::fabs(zb) / THICK_H2) / (1.0 + 0.5 * nnf)); //Msun/pc^3
        }
        s.rho_thick[i] *= THICK_NORM;

///========== Galactic Stellar Halo=================
        rdi = std::sqrt(Rb * Rb + zb * zb / (HALO_FLATTEN * HALO_FLATTEN));
        if (rdi <= HALO_CORE) {
            s.rho_halo[i] = std::fabs(1.0 * HALO_RHO0 * std::pow(HALO_CORE / Dsun, HALO_SLOPE));
        }
        else {
            s.rho_halo[i] = std::fabs(1.0 * HALO_RHO0 * std::pow(rdi / Dsun, HALO_SLOPE)); //Msun/pc^3
        }
        s.rho_halo[i] *= HALO_NORM;

///========== Galactic bulge =====================
        constexpr double barMassRescale = BAR_MASS_RESCALE;   // see config/parameters.h, section 2a
        xf =  xb * std::cos(alfa) + yb * std::sin(alfa);
        yf = -xb * std::sin(alfa) + yb * std::cos(alfa);
        zf =  zb;

        Rx0    = BAR_S_RX0;
        Ry0    = BAR_S_RY0;
        Rz0    = BAR_S_RZ0;
        Rc     = BAR_S_RC;
        cp     = BAR_S_CP;
        cn     = BAR_S_CN;
        mBarre = BAR_S_MASS_NUM / BAR_S_MASS_DEN * barMassRescale;

        r4  = std::pow(std::pow(std::fabs(xf / Rx0), cn)
                     + std::pow(std::fabs(yf / Ry0), cn), cp / cn)
            + std::pow(std::fabs(zf / Rz0), cp);
        r4  = std::pow(std::fabs(r4), 1.0 / cp);
        r2  = std::sqrt(std::fabs(xf * xf + yf * yf));

        if (r2 <= Rc) {
            rhoS = mBarre * 1.0 / (std::cosh(-r4) * std::cosh(-r4));
        }
        else {
            rhoS = mBarre * 1.0 / (std::cosh(-r4) * std::cosh(-r4)) * std::exp(-BAR_CUTOFF_K * (r2 - Rc) * (r2 - Rc));
        }

        Rx0    = BAR_E_RX0;
        Ry0    = BAR_E_RY0;
        Rz0    = BAR_E_RZ0;
        Rc     = BAR_E_RC;
        cp     = BAR_E_CP;
        cn     = BAR_E_CN;
        mBarre = BAR_E_MASS_NUM / BAR_E_MASS_DEN * barMassRescale; // normalized

        r4  = std::pow(std::fabs(std::pow(std::fabs(xf / Rx0), cn)
                               + std::pow(std::fabs(yf / Ry0), cn)), cp / cn)
            + std::pow(std::fabs(zf / Rz0), cp);
        r4  = std::pow(r4, 1.0 / cp);
        r2  = std::sqrt(std::fabs(xf * xf + yf * yf));

        if (r2 <= Rc) {
            rhoE = mBarre * std::exp(-r4);
        }
        else {
            rhoE = mBarre * std::exp(-r4) * std::exp(-BAR_CUTOFF_K * (r2 - Rc) * (r2 - Rc));
        }

        s.rho_bulge[i]  = std::fabs(rhoS) + std::fabs(rhoE); ///Msun/pc^3
        s.rho_bulge[i] *= BAR_NORM;
///==================================================================

        s.Rostar0[i] = std::fabs(s.rho_thin[i] + s.rho_thick[i] + s.rho_bulge[i] + s.rho_halo[i]); //[Msun/pc^3]
        s.Rostari[i] = s.Rostar0[i] * x * x * step * 1.0e9 * (M_PI / 180.0) * (M_PI / 180.0); //[Msun/deg^2]
        s.Nstari[i]  = binary_fraction
                     * (s.rho_thin[i] * fd / MEANMASS_THIN
                      + s.rho_thick[i] * fh / MEANMASS_THICK
                      + s.rho_halo[i] * fh / MEANMASS_HALO
                      + s.rho_bulge[i] * fb / MEANMASS_BULGE); //[Nt/pc^3]
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
