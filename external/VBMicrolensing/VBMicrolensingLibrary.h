// VBMicrolensing v5.6 (2026)
//
// This code has been developed by Valerio Bozza (University of Salerno) and collaborators.
// Check the repository at https://github.com/valboz/VBMicrolensing
// for the newest version.
// Any use of the code for scientific publications should be acknowledged by a citation
// to the appropriate publication, as detailed in the repository page.
//
// The code relies on the root solving algorithm by Jan Skworon and Andy Gould
// described in Skowron & Gould arXiv:1203.1034.
// Please also cite this paper if specifically relevant in your scientific publication.
// The original Fortran code available on http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/
// has been translated to C++ by Tyler M. Heintz and Ava R. Hoag (2017)
//
// The Multipoly method for the calculation of multiple-lens microlensing was developed
// with Vito Saggese (2024).
// 
// Optimizations for complex functions and for high-mag regime have been developed 
// by Jiyuan Zhang (2025).
//
// GNU Lesser General Public License applies to all parts of this code.
// Please read the separate LICENSE.txt file for more details.


#ifndef __multilens
#define __multilens
#define __unmanaged

#define _L_1 x1-((x1+a/2.0)/((x1+a/2.0)*(x1+a/2.0)+x2*x2)+q*(x1-a/2.0)/((x1-a/2.0)*(x1-a/2.0)+x2*x2))/(1.0+q) // Used in PlotCrits
#define _L_2 x2-(x2/((x1+a/2.0)*(x1+a/2.0)+x2*x2)+q*x2/((x1-a/2.0)*(x1-a/2.0)+x2*x2))/(1.0+q)
#define _LL (y-z)+coefs[21]/(zc-coefs[20])+coefs[22]/zc //Lens equation test
#define _J1c coefs[21]/((zc-coefs[20])*(zc-coefs[20]))+coefs[22]/(zc*zc) //#define _J1 m1/((zc-0.5*a)*(zc-0.5*a))+m2/((zc+0.5*a)*(zc+0.5*a))
#define _J2 -2.0*(coefs[21]/((z-coefs[20])*(z-coefs[20])*(z-coefs[20]))+coefs[22]/(z*z*z))
#define _J3 6.0*(coefs[21]/((z-coefs[20])*(z-coefs[20])*(z-coefs[20])*(z-coefs[20]))+coefs[22]/(z*z*z*z))
#define _skew(p1,p2,q1,q2) p1*q2-p2*q1
#define _NP 200.0
#define __rsize_ESPL 151
#define __zsize_ESPL 101

#define _sign(x) ((x>0)? +1 : -1)

#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <random>
class _sols_for_skiplist_curve;
class _skiplist_curve;
/*******************************************   end   *******************************************/

class _curve;
class _sols;
class _theta;
class complex;
struct annulus;


class complex {
public:
	double re;
	double im;
	complex(double a, double b) { re = a; im = b; }
	complex(double a) { re = a; im = 0; }
	complex(void) { re = 0; im = 0; }
};



class VBMicrolensing
{
	const int maxiter = 10000;
	const double pseudorandom[12] = { 0.5, 0.4, 0.3, 0.7, 0.6, 0.5, 0.4, 0.7, 0.5, 0.6, 0.3,
		0.3 };

	int* ndatasat;
	double** tsat, *** possat;
	double** posEar, startEar, stepEar;
	int ndataEar;
	double Mag0;
	double* dist_mp, * q;
	int nim0, n, n2, nnm1, nroots, nrootsmp, * nrootsmp_mp;
	complex* zr, * zcr, ** pmza, ** pyaza, ** ppmy, * pza, * pza2, ** pmza2, * pdum, * ppy, * a, * s_offset, * pert, y, yc, * s;
	complex* y_mp, *** pmza_mp, ** pza_mp, *** pyaza_mp, *** ppmy_mp, ** ppy_mp, ** zr_mp;
	complex* zaltc, * J1, * J1c, ** za, ** za2;
	complex* coefs, ** coefs_mp;
	complex** a_mp, * s_sort;

	double* prodevs, * errs, err, L0f, Jacf;
	complex* devs, * init, * centralimages, * newseeds, * grads, zf, S2f, * S2s, * S3s, * S4s;
	int lencentralimages, lennewseeds, ngoodold, ngood, iter, iter2;
	////
	double* good, * Jacs, rho, rho2, * m;
	double** m_mp, * q_sort;
	int* worst;
	double e, phi, phip, phi0, Om, inc, u0, tE_inv, t0, alpha, pai1, pai2, PosAng, dPosAng, thetaE, d3, v3, GM, flagits;
	int iastro;
	double Obj[3], rad[3], tang[3], t0old, t0parold;
	double Eq2000[3], Quad2000[3], North2000[3];
	double Et0[2], vt0[2], Et[2], Ehel[2], lighttravel, lighttravel0, lighttravel0par;
	double ESPLout[__rsize_ESPL][__zsize_ESPL], ESPLin[__rsize_ESPL][__zsize_ESPL], ESPLoutastro[__rsize_ESPL][__zsize_ESPL], ESPLinastro[__rsize_ESPL][__zsize_ESPL];
	bool coordinates_set;
	bool multidark;
	double* LDtab, * rCLDtab, * CLDtab;
	double scr2, sscr2;
	int npLD;
	annulus* annlist;
	_skiplist_curve** cprec, ** cpres, ** cfoll;
	double** A;

	void ComputeCentroids(double* pr, double t, double* c1s, double* c2s, double* c1l, double* c2l);
	void ComputeParallax(double, double);
	double du_par(double t);
	double LDprofile(double r);
	double rCLDprofile(double tc, annulus*, annulus*);
	void initroot();
	int froot(complex);
	bool checkroot(_theta*);

	void SetLensGeometry_spnp(int n, double* q, complex* s);
	void SetLensGeometry_multipoly(int n, double* q, complex* s);
	void initrootpoly();
	_curve* NewImages(complex, complex*, _theta*);
	_curve* NewImages(_theta*);
	_curve* NewImagespoly(_theta*);
	_curve* NewImagesmultipoly(_theta*);
	_curve* NewImages_shear(complex, complex*, _theta*);
	double BinaryMagSafe(double s, double q, double y1, double y2, double rho, _sols_for_skiplist_curve** images);
	double MultiMagSafe(double y1, double y2, double rho, _sols_for_skiplist_curve** images);
	void OrderImages(_sols_for_skiplist_curve*, _curve*);

	void OrderMultipleImages(_sols_for_skiplist_curve*, _curve*);
	void cmplx_laguerre(complex*, int, complex*, int&, bool&);
	void cmplx_newton_spec(complex*, int, complex*, int&, bool&);
	void cmplx_laguerre2newton(complex*, int, complex*, int&, bool&, int);
	void solve_quadratic_eq(complex&, complex&, complex*);
	void solve_cubic_eq(complex&, complex&, complex&, complex*);
	void polyproduct(complex* p1, int n1, complex* p2, int n2, complex* pdest);
	void copypol(complex* p1, int n1, complex* pdest);
	void change_n(int nn);
	void change_n_mp(int nn);
	void polycoefficients();
	void polycoefficients_multipoly();
	void polycritcoefficients(complex eiphi);

public:

	double rootaccuracy;
	double samplingfactor;
	bool squarecheck;
	bool astrometry;
	bool turn_off_secondary_source;
	bool turn_off_secondary_lens;
	bool block_tertiary_lens;
	bool ESPLoff;
	bool t_in_HJD;

	static char ESPLtablefile[1024];
	static void SetESPLtablefile(char* instring) { strcpy(ESPLtablefile, instring); }
	static char Suntablefile[1024];
	static void SetSuntablefile(char* instring) { strcpy(Suntablefile, instring); }
	double Tol, RelTol, a1, a2, corrquad, corrquad2, safedist;
	double mass_radius_exponent, mass_luminosity_exponent, lens_mass_luminosity_exponent;
	int satellite, parallaxsystem, t0_par_fixed, nsat;
	double t0_par;
	void t0_from_t0_par(double t0, double tE, double u0, double pai1, double pai2);
	double t0_out, u0_out, alpha_out, pai1_out, pai2_out, tE_out;
	bool suntable, parallaxephemeris;
	int parallaxextrapolation;
	int minannuli, maxannuli, nannuli, NPS, NPcrit;
	int newtonstep;
	double y_1, y_2, av, therr, astrox1, astrox2;
	double (*CumulativeFunction)(double r, double* LDpars);

	// Critical curves and caustics calculation
	_sols* PlotCrit();
	_sols* PlotCrit(double a, double q);
	// Initialization for parallax calculation
	void SetObjectCoordinates(char* Coordinates_file, char* Directory_for_satellite_tables);
	void SetObjectCoordinates(char* CoordinateString);
	bool AreCoordinatesSet();

	// ====================================================================================================
	// LOCAL ADDITION to VBMicrolensing v5.6.1 (upstream commit db4f95267d74): NOT part of the library as
	// distributed by its authors. Added 2026-10-11 for the Roman-Rubin simulator, which vendors this copy.
	// This is one of two local changes to the upstream files (the other is LOCAL ADDITION 2, below); remove the
	// method and that part of the library is unmodified.
	//
	// EarthVelocityOnSky(t0, v) returns in v[0], v[1] the Earth's velocity at time t0 projected on the sky
	// at the target set by SetObjectCoordinates, along celestial South and West, in AU/day (t0 in the
	// units of the light curves: JD - 2450000). It is exactly the member vt0 that ComputeParallax builds
	// for t0_par = t0, the epoch of the geocentric frame (Gould 2004) in which every *Parallax light
	// curve is computed; vt0 is private and not otherwise reachable.
	//
	// Why: the simulator draws events with heliocentric proper motions, while the parallax light curves
	// are in the geocentric frame, where the lens-source proper motion differs by the Earth's projected
	// velocity times the relative parallax (mu_geo = mu_hel - v_E,perp * pi_rel; the same relation
	// ComputeCentroids uses). Converting needs v_E,perp at t0, from the very table the light curve uses.
	//
	// Needs SetObjectCoordinates called, and satellite <= nsat (otherwise ComputeParallax returns without
	// computing and v is not set). The satellite selection does not enter vt0. The method restores
	// t0_par_fixed and t0_par and leaves the parallax caches (t0old, t0parold) invalidated, so the next
	// light-curve call recomputes everything; it changes no result of any existing function.
	// ====================================================================================================
	void EarthVelocityOnSky(double t0, double* v) {
		const int keepFixed = t0_par_fixed;
		const double keepPar = t0_par;
		t0old = t0parold = 1.e200;	// the parallax caches depend on the target: drop them
		t0_par_fixed = 0;			// t0_par = t0 inside ComputeParallax
		ComputeParallax(t0, t0);	// builds vt0, the geocentric frame's Earth velocity at t0
		t0_par_fixed = keepFixed;
		t0_par = keepPar;
		v[0] = vt0[0]; v[1] = vt0[1];	// South, West [AU/day]
		t0old = t0parold = 1.e200;
	}
	// ====================================================================================================
	// ====================================================================================================
	// LOCAL ADDITION 2 (2026-11-13) to VBMicrolensing v5.6.1: NOT part of the library as distributed by its
	// authors. Added for the Roman-Rubin simulator, which vendors this copy. Only lines were added; none of
	// upstream's was changed or removed (external/VBMicrolensing/README.md, "Local modifications").
	//
	// ComputeParallax(t, t0) works out, at every call, where the observers are at time t: the Earth's
	// ecliptic position Ear(t), interpolated linearly in the daily Sun table, and, for a satellite, its
	// position Spit(t), found by a binary search of its table and interpolated. Neither depends on the
	// target or on the event (t0, tE, ...), only on t and on the loaded tables, and a simulation evaluates
	// the same epochs thousands of times. SetObserverGrid(t, n) computes them once for n strictly increasing
	// times t[] (JD - 2450000, the units of the light curves), by the very expressions ComputeParallax uses,
	// and stores them with the value those expressions would give parallaxextrapolation (1: outside the Sun
	// table, 2: outside the satellite's). ComputeParallax then, when t is a grid time bit for bit, writes
	// the stored doubles where it would have computed them and runs every other line unchanged and in its
	// order (the light-travel term, the projections on rad and tang, the geocentric-frame offsets Et): the
	// results are identical, bit for bit. A t that is not on the grid takes upstream's lines.
	//
	// Used only when SetParallaxCache(true) has been called (off by default: the library is then upstream's),
	// parallaxephemeris is set and t_in_HJD is false. LoadSunTable and the two-argument SetObjectCoordinates
	// (which reload the tables) empty the grid; a change of target does not, nothing stored depends on it.
	// An unsorted or non-finite t[] leaves the grid empty. Not thread-safe, like the rest of the library.
	// pc_hits and pc_misses count the lookups, for tests.
	// ====================================================================================================
	long long pc_hits = 0, pc_misses = 0;
	void SetParallaxCache(bool on) { pc_on = on; }
	void SetObserverGrid(const double* t, int n) {
		PCClear();
		if (n <= 0) return;
		for (int i = 0; i < n; i++) {
			if (!(t[i] - t[i] == 0.0)) return;			// NaN or infinite
			if (i > 0 && !(t[i] > t[i - 1])) return;	// not strictly increasing
		}
		if (!suntable) LoadSunTable(Suntablefile);	// as ComputeParallax does
		if (!suntable) return;
		pc_t.assign(t, t + n);
		pc_ear.assign(3 * (size_t)n, 0.0);
		pc_flagE.assign((size_t)n, 0);
		for (int k = 0; k < n; k++) {
			// Exactly the Earth's interpolation of ComputeParallax, clamps included.
			double ty = (t[k] - startEar) / stepEar;
			int ic = (int)floor(ty);
			ty -= ic;
			if (ic < 0) { pc_flagE[k] = 1; ic = 0; ty = 0; }
			if (ic >= ndataEar) { pc_flagE[k] = 1; ic = ndataEar - 1; ty = 1; }
			for (int i = 0; i < 3; i++) pc_ear[3 * (size_t)k + i] = posEar[ic][i] * (1 - ty) + posEar[ic + 1][i] * ty;
		}
		pc_sp.assign((size_t)nsat, std::vector<double>());
		pc_flagS.assign((size_t)nsat, std::vector<unsigned char>());
		for (int s = 0; s < nsat; s++) {
			if (ndatasat[s] <= 2) continue;	// ComputeParallax skips such a satellite
			pc_sp[s].assign(3 * (size_t)n, 0.0);
			pc_flagS[s].assign((size_t)n, 0);
			for (int k = 0; k < n; k++) {
				// Exactly the satellite's search and interpolation of ComputeParallax.
				int ic, left, right;
				if (t[k] < tsat[s][0]) {
					ic = 0;
					pc_flagS[s][k] = 1;
				}
				else {
					if (t[k] > tsat[s][ndatasat[s] - 1]) {
						ic = ndatasat[s] - 2;
						pc_flagS[s][k] = 1;
					}
					else {
						left = 0;
						right = ndatasat[s] - 1;
						while (right - left > 1) {
							ic = (right + left) / 2;
							if (tsat[s][ic] > t[k]) {
								right = ic;
							}
							else {
								left = ic;
							}
						}
						ic = left;
					}
				}
				double ty = (t[k] - tsat[s][ic]) / (tsat[s][ic + 1] - tsat[s][ic]);
				for (int i = 0; i < 3; i++) pc_sp[s][3 * (size_t)k + i] = possat[s][ic][i] * (1 - ty) + possat[s][ic + 1][i] * ty;
			}
		}
		pc_cur.assign((size_t)nsat + 1, 0);
		pc_n = n;
	}
private:
	bool pc_on = false;
	int pc_n = 0;						// number of grid times; 0 = no grid
	std::vector<double> pc_t, pc_ear;	// grid times; Ear(t) of each, 3 per time
	std::vector<unsigned char> pc_flagE;	// 1 where the Earth's interpolation was clamped (parallaxextrapolation = 1)
	std::vector<std::vector<double> > pc_sp;	// per satellite: Spit(t) of each grid time, 3 per time
	std::vector<std::vector<unsigned char> > pc_flagS;	// 1 where the satellite's table was extrapolated (= 2)
	std::vector<int> pc_cur;			// per observer (0 = Earth, k = satellite k): the last grid index found
	void PCClear() {
		pc_n = 0;
		pc_t.clear(); pc_ear.clear(); pc_flagE.clear(); pc_sp.clear(); pc_flagS.clear(); pc_cur.clear();
	}
	static bool PCSame(double a, double b) {
		unsigned long long x, y;
		memcpy(&x, &a, sizeof x);
		memcpy(&y, &b, sizeof y);
		return x == y;
	}
	// The grid index of the time t, found by the equality of its bits, or -1. The next grid time after the
	// last one found (the order of a time loop) costs one comparison; a later time is found by galloping
	// forward from there; an earlier one by bisection of everything before.
	int PCFind(double t, int observer) {
		if (pc_n == 0 || observer < 0 || observer >= (int)pc_cur.size()) { pc_misses++; return -1; }
		const double* g = pc_t.data();
		const int n = pc_n;
		const int cur = pc_cur[observer];
		int lo, hi;
		if (t >= g[cur]) {
			if (PCSame(g[cur], t)) { pc_hits++; return cur; }
			if (cur + 1 < n && PCSame(g[cur + 1], t)) { pc_hits++; return pc_cur[observer] = cur + 1; }
			lo = cur;
			int step = 1;
			hi = cur + 1;
			while (hi < n && g[hi] < t) { lo = hi; step *= 2; hi = lo + step; }
			if (hi > n) hi = n;
		}
		else {
			lo = 0;
			hi = cur;
		}
		// the first index in [lo, hi) with g >= t, or hi (for a NaN t every comparison fails: lo)
		while (lo < hi) {
			const int mid = lo + (hi - lo) / 2;
			if (g[mid] < t) lo = mid + 1; else hi = mid;
		}
		if (lo < n && PCSame(g[lo], t)) { pc_hits++; return pc_cur[observer] = lo; }
		pc_misses++;
		return -1;
	}
public:
	// END LOCAL ADDITION 2
	// ====================================================================================================
	// Skowron & Gould root calculation
	void cmplx_roots_gen(complex*, complex*, int, bool, bool);
	void cmplx_roots_multigen(complex*, complex**, int, bool, bool);
	// Bozza optimization
	int findimagepoly(int iroot);
	int findimagemultipoly(int iroot);

	// Set Lens Geometry
	void SetLensGeometry(int n, double* q, complex* s);
	void SetLensGeometry(int n, double* pr);

	// Magnification calculation functions.

	double BinaryMag0(double s, double q, double y1, double y2, _sols_for_skiplist_curve** Images);
	double BinaryMag0(double s, double q, double y1, double y2);

	double BinaryMag(double s, double q, double y1, double y2, double rho, double accuracy, _sols_for_skiplist_curve** Images);
	double BinaryMag(double s, double q, double y1, double y2, double rho, double accuracy);
	double BinaryMag2(double s, double q, double y1, double y2, double rho);
	double BinaryMagDark(double s, double q, double y1, double y2, double rho, double accuracy);
	void BinaryMagMultiDark(double s, double q, double y1, double y2, double rho, double* a1_list, int n_filters, double* mag_list, double accuracy);
	double MultiMag0(double y1, double y2, _sols_for_skiplist_curve** Images);
	double MultiMag0(double y1, double y2);
	double MultiMag(double y1, double y2, double rho, double accuracy, _sols_for_skiplist_curve** Images);
	double MultiMag(double y1, double y2, double rho, double accuracy);
	double MultiMag(double y1, double y2, double rho);
	double MultiMag2(double y1, double y2, double rho);
	double MultiMagDark(double y1, double y2, double rho, double accuracy);

	double BinaryMag0_shear(double s, double q, double y1, double y2, double K1, double G1, double Gi, _sols** Images);
	double BinaryMag0_shear(double s, double q, double y1, double y2, double K1, double G1, double Gi);

	// Limb Darkening control
	enum LDprofiles { LDlinear, LDquadratic, LDsquareroot, LDlog, LDuser };
	void SetLDprofile(double(*UserLDprofile)(double), int tablesampling);
	void SetLDprofile(LDprofiles);

	// Method control
	enum class Method { Singlepoly, Multipoly, Nopoly };
	void SetMethod(Method);

	//ESPL functions
	void LoadESPLTable(const char* tablefilename);
	double ESPLMag(double u, double rho);
	double ESPLMag2(double u, double rho);
	double ESPLMagDark(double u, double rho);
	double PSPLMag(double u);


	// New (v2) light curve functions, operating on arrays

	void LoadSunTable(char* tablefilename);
	void PSPLLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void PSPLLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void ESPLLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void ESPLLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);

	void BinaryLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinaryLightCurveW(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinaryLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinaryLightCurveOrbital(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* sep_array, int np);
	void BinaryLightCurveKepler(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* sep_array, int np);

	void BinSourceLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinSourceLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinSourceLightCurveXallarap(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* sep_array, int np);
	void BinSourceExtLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinSourceExtLightCurveXallarap(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* y1_array2, double* y2_array2, int np);
	void BinSourceSingleLensXallarap(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* y1_array2, double* y2_array2, int np);
	void BinSourceBinLensXallarap(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void BinSourceBinLensLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* y1_array2, double* y2_array2, double* seps_array, int np);

	void TripleLightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void TripleLightCurveParallax(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np);
	void TripleLightCurveOrbital(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, double* seps_array, double* seps2_array, double* psi_array, int np);
	void LightCurve(double* parameters, double* t_array, double* mag_array, double* y1_array, double* y2_array, int np, int nl);

	// Old (v1) light curve functions, for a single calculation
	double PSPLLightCurve(double* parameters, double t);
	double PSPLLightCurveParallax(double* parameters, double t);
	double ESPLLightCurve(double* parameters, double t);
	double ESPLLightCurveParallax(double* parameters, double t);

	double BinaryLightCurve(double* parameters, double t);
	double BinaryLightCurveW(double* parameters, double t);
	double BinaryLightCurveParallax(double* parameters, double t);
	double BinaryLightCurveOrbital(double* parameters, double t);
	double BinaryLightCurveKepler(double* parameters, double t);

	double BinSourceLightCurve(double* parameters, double t);
	double BinSourceLightCurveParallax(double* parameters, double t);
	double BinSourceLightCurveXallarap(double* parameters, double t);
	double BinSourceExtLightCurve(double* parameters, double t);
	double BinSourceExtLightCurveXallarap(double* parameters, double t);
	double BinSourceBinLensXallarap(double* parameters, double t);
	double BinSourceSingleLensXallarap(double* parameters, double t);
	double BinSourceBinLensLightCurve(double* parameters, double t);


	double TripleLightCurve(double* parameters, double t);
	double TripleLightCurveParallax(double* parameters, double t);
	double TripleLightCurveOrbital(double* parameters, double t);

	// Astrometric functions
	void CombineCentroids(double* mags, double* c1s, double* c2s, double* c1l, double* c2l, double* c1ltot, double* c2tot, double g, int np);
	void PSPLAstroLightCurve(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, int np);
	void ESPLAstroLightCurve(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, int np);
	void BinaryAstroLightCurve(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, int np);
	void BinaryAstroLightCurveOrbital(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, double* seps_array, int np);
	void BinaryAstroLightCurveKepler(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, double* seps_array, int np);
	void BinSourceAstroLightCurveXallarap(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, double* y1_array2, double* y2_array2, int np);
	void BinSourceBinLensAstroLightCurve(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, double* y1_array2, double* y2_array2, double* seps_array, int np);
	void TripleAstroLightCurve(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, int np);
	void TripleAstroLightCurveOrbital(double* parameters, double* t_array, double* mag_array, double* centroid_s1_array, double* centroid_s2_array, double* centroid_l1_array, double* centroid_l2_array, double* y1_array, double* y2_array, double* seps_array, double* seps2_array, double* psi_array, int np);

	// Constructor and destructor

	VBMicrolensing();
	~VBMicrolensing();

private: // Must be declared here at the end
	LDprofiles curLDprofile;
	Method SelectedMethod;

};

//std::string VBMicrolensing::ESPLtablefile = "";

double VBDefaultCumulativeFunction(double r, double* a1);

struct annulus {
	double bin;
	double cum;
	double Mag;
	double err;
	double f;
	int nim;
	double LDastrox1, LDastrox2;
	annulus* prev, * next;
};



class _theta {
public:
	double th, maxerr, Mag, errworst, astrox1, astrox2;
	int imlength;
	_theta* prev, * next;

	_theta(double);

};

class _thetas {
public:
	_theta* first, * last;
	int length;

	_thetas(void);
	~_thetas(void);
	_theta* insert(double);
	_theta* insert_at_certain_position(_theta*, double);
	// method: this method can only be used when inserting an element in the middle of linked list
	// 		   i.e. *first's 'th' < current 'th' < *last's 'th'
	// 		   and the new element is forced to be inserted between itheta and itheta->next, 
	// 		   which means it's the programmer's responsibility to guarantee itheta->th < th < itheta->next->th holds
	//         (O(1) complexity)
	void remove(_theta*);

};

#define max_skiplist_level 2

class _point {
public:
	double x1;
	double x2;
	double parab, ds, dJ, Mag, err, parabastrox1;
	complex d;								  	  // d is z'(theta) at this point
	_theta* theta;							  // pointer
	_point* next, * prev;                       // pointers that point to _point variable
	_point* next_array[max_skiplist_level + 1];

	double parabastrox2;
	_point(double, double, _theta*);
	double operator-(_point);
};

class _curve {
public:
	int length;
	_point* first, * last;
	_curve* next, * prev;
	_curve* partneratstart, * partneratend;
	double parabstart, Magstart, errstart, parabastrox1, parabastrox2;

	_curve(_point*);
	_curve(void);
	~_curve(void);

	_curve* divide(_point*);
	void drop(_point*);
	void append(double, double);
	void append(_point*);
	void prepend(double, double);
	//	void prepend(_point *);
	_curve* join(_curve*);
	_curve* joinbefore(_curve*);
	_curve* reverse(void);
	double closest(_point*, _point**);
	double closest2(_point*, _point**);
	void complement(_point**, int, _point**, int);
};

class _sols {
public:
	int length;
	_curve* first, * last;

	_sols(void);
	~_sols(void);
	void drop(_curve*);
	void append(_curve*);
	void prepend(_curve*);
	void join(_sols*);
};


class _skiplist_curve {							 // a _skiplist_curve class variable is a skip list of _point variables
public:
	_point* first, * last;

	_point* head;
	_point* last_array[max_skiplist_level + 1];
	int Level;

	int length_notation;

	_skiplist_curve* next, * prev;
	_skiplist_curve* partneratstart, * partneratend;
	double parabstart, Magstart, errstart, parabastrox1, parabastrox2;

	_skiplist_curve(_point* p1, int new_Level);
	_skiplist_curve(void);
	~_skiplist_curve(void);
	_skiplist_curve* join(_skiplist_curve* new_curve);
	void append(_point* pp, int append_Level);
	void append(double x1, double x2, int append_Level);
	_skiplist_curve* find_prev_then_divide(double th);

};




class _sols_for_skiplist_curve {                 //       a  _sols_for_skiplist_curve class variable is a linked list of _skiplist_curve variables, 
	// while a _skiplist_curve class variable           is a skip list   of _point variables
public:
	int length;
	_skiplist_curve* first, * last;

	_sols_for_skiplist_curve(void);
	~_sols_for_skiplist_curve(void);
	void drop(_skiplist_curve* ref);
	void append(_skiplist_curve* cc);

};

#define MR 8
#define MT 10
#define MAXIT (MT*MR)
#define MAXM 101

//////////////////////////////
//////////////////////////////
////////complex methods and operators
//////////////////////////////
//////////////////////////////

inline double abs2(complex z) {
	return (z.re * z.re + z.im * z.im);
}

inline double abs(complex z) {
	return sqrt(z.re * z.re + z.im * z.im);
}

inline complex conj(complex z) {
	return complex(z.re, -z.im);
}

inline complex sqrt(complex z) {
	double md = sqrt(z.re * z.re + z.im * z.im);
	return (md > 0) ? complex(sqrt((md + z.re) / 2), (sqrt((md - z.re) / 2) * ((z.im > 0) ? 1 : -1))) : 0.0;
}



inline double real(complex z) {
	return z.re;
}

inline double imag(complex z) {
	return z.im;
}

inline complex operator+(complex p1, complex p2) {
	return complex(p1.re + p2.re, p1.im + p2.im);
}

inline complex operator-(complex p1, complex p2) {
	return complex(p1.re - p2.re, p1.im - p2.im);
}

inline complex operator*(complex p1, complex p2) {
	return complex(p1.re * p2.re - p1.im * p2.im, p1.re * p2.im + p1.im * p2.re);
}

inline complex operator/(complex p1, complex p2) {
	double md = p2.re * p2.re + p2.im * p2.im;
	return complex((p1.re * p2.re + p1.im * p2.im) / md, (p1.im * p2.re - p1.re * p2.im) / md);
}

inline complex operator+(complex z, double a) {
	return complex(z.re + a, z.im);
}

inline complex operator-(complex z, double a) {
	return complex(z.re - a, z.im);
}

inline complex operator*(complex z, double a) {
	return complex(z.re * a, z.im * a);
}

inline complex operator/(complex z, double a) {
	return complex(z.re / a, z.im / a);
}

inline complex operator+(double a, complex z) {
	return complex(z.re + a, z.im);
}

inline complex operator-(double a, complex z) {
	return complex(a - z.re, -z.im);
}

inline complex operator*(double a, complex z) {
	return complex(a * z.re, a * z.im);
}

inline complex operator/(double a, complex z) {
	double md = z.re * z.re + z.im * z.im;
	return complex(a * z.re / md, -a * z.im / md);
}

inline complex operator+(complex z, int a) {
	return complex(z.re + a, z.im);
}

inline complex operator-(complex z, int a) {
	return complex(z.re - a, z.im);
}

inline complex operator*(complex z, int a) {
	return complex(z.re * a, z.im * a);
}

inline complex operator/(complex z, int a) {
	return complex(z.re / a, z.im / a);
}

inline complex operator+(int a, complex z) {
	return complex(z.re + a, z.im);
}

inline complex operator-(int a, complex z) {
	return complex(a - z.re, -z.im);
}

inline complex operator*(int a, complex z) {
	return complex(a * z.re, a * z.im);
}

inline complex operator/(int a, complex z) {
	double md = z.re * z.re + z.im * z.im;
	return complex(a * z.re / md, -a * z.im / md);
}

inline complex operator-(complex z) {
	return complex(-z.re, -z.im);
}

inline bool operator==(complex p1, complex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return true;
	return false;
}

inline bool operator!=(complex p1, complex p2) {
	if (p1.re == p2.re && p1.im == p2.im) return false;
	return true;
}

inline complex expcmplx(complex p1) {
	double r = exp(p1.re);
	double theta = atan2(p1.im, p1.re);
	return complex(r * cos(theta), r * sin(theta));
}

inline complex cbrt(complex z) {
	complex zout;
	double r, r_cube, theta, theta_cube;
	r = abs(z);
	r_cube = pow(r, 0.333333333333);
	theta = atan2(z.im, z.re);
	theta_cube = theta / 3.;
	return 	complex(r_cube * cos(theta_cube), r_cube * sin(theta_cube));
}


inline _point::_point(double x, double y, _theta* theta1) {
	x1 = x;
	x2 = y;
	theta = theta1;
	next = 0;
	prev = 0;
	for (int i = 0; i < (max_skiplist_level + 1); i++)
	{
		next_array[i] = 0;
	}
}

inline double _point::operator-(_point p2) {
	static double dx1, dx2;
	dx1 = x1 - p2.x1;
	dx2 = x2 - p2.x2;
	return dx1 * dx1 + dx2 * dx2;
}


#endif


