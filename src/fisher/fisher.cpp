// FisherM (information matrices by finite differences) and ErrorCal (sigmas from the inverses).
#include "fisher/fisher.h"
#include "events/lightcurve.h"
#include "surveys/noise.h"


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
    co.Delta1[0] = FD_STEP_U0;//u0
    co.Delta1[1] = l.tE * FD_STEP_TE_FRAC;//tE [days]
    co.Delta1[3] = FD_STEP_PIE_FRAC * l.piE;//piE
    co.Delta1[4] = FD_STEP_XI_DEG * M_PI / 180.0; //for xi [radian]
    co.Delta1[5] = l.tE * FD_STEP_TE_FRAC; // t0 [days] -- same fractional-of-tE step as tE itself,
    // since t0's local curvature scale is set by the same characteristic timescale.
    // Placeholder pending Step C3's step-size convergence sweep, like the other four.
    co.Delta1[6] = FD_STEP_MBS; // mbs0 [mag] -- Rubin baseline magnitude
    co.Delta1[8] = FD_STEP_MBS; // mbs1 [mag] -- Roman baseline magnitude

    // Step C3. The five constants above came from the legacy LMC codebase and are ~25% of the
    // parameter value -- u0 is perturbed by 0.15 on a u0 of ~0.3, tE and t0 by 0.25*tE. That is
    // enormous for a derivative step, and the sweep (./fishertest --sweep, plotted by
    // tests/c3_step_sweep.py) showed every nonlinear parameter sitting far up the
    // truncation-error branch: sigma(u0) off by ~57% median and sigma(t0) by ~71%.
    //
    // Scaling them into the convergence plateau, where sigma is flat to <0.2% over the decades
    // 1e-6..1e-3 of the legacy values. Below ~1e-6 round-off takes over (subtracting two model
    // magnitudes that agree to ~1e-9 leaves too few significant digits); above ~1e-3 truncation
    // does. 1e-4 sits two decades clear of the round-off wall, and 1e-5 reproduces every sigma
    // to 0.3%, so the choice is not delicate.
    //
    // Expressed as a scale factor rather than five rewritten constants so the change stays
    // auditable against the legacy values and the sweep -- which is defined in these units --
    // remains directly comparable.
    //
    // Consequence worth knowing: this makes the sigmas MUCH larger for short events, and that is
    // the correct answer, not a regression. The oversized steps were manufacturing information
    // in the near-degenerate directions of the Fisher matrix -- for a 5-day event the smallest
    // normalized eigenvalue was inflated by ~4e7, reporting a 0.2% parallax measurement from a
    // light curve far too short to measure annual parallax at all. Per-event condition numbers
    // rise correspondingly and are now genuine. The joint-vs-single-survey RATIOS, which the
    // thesis rests on, are unchanged in ordering and magnitude (Step C5's paired comparison
    // cancels the common inflation); see DEVIATIONS.md.
    // kFDStepScale (= 1.0e-4): see config/parameters.h, section 7.
    for (int q = 0; q < Nx; ++q) co.Delta1[q] *= kFDStepScale * co.deltaScale[q];
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



    // Deviation 76: a telescope with fewer than kMinTeleEpochs epochs is left out (see config/parameters.h).
    std::array<int, 2> nTele{0, 0};
    for (int i = 0; i < ndw; ++i) nTele[int(l.tele[i]) == 1 ? 1 : 0] += 1;

    for (int i = 0; i < ndw; ++i) {//data
        tt = int(l.tele[i]);//telescope[0,1] LSST, ELT
        if (nTele[tt == 1 ? 1 : 0] < kMinTeleEpochs) continue;
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
        if (s.fb[tt] < FB_BIN_LO)      {co.bb[0] =+ FB_STEP_SMALL; co.bb[1] =+ FB_STEP_LARGE;}
        else if (s.fb[tt] < FB_BIN_HI) {co.bb[0] =- FB_STEP_SMALL; co.bb[1] =+ FB_STEP_SMALL;}
        else                      {co.bb[0] =- FB_STEP_SMALL; co.bb[1] =- FB_STEP_LARGE;}

        // Step C3: apply the same plateau scaling as Delta1[] above, then clamp so the step
        // still cannot leave the physical range [0,1]. The bin edges above guarantee
        // bound-safety only at the unscaled sizes, so a scaled-UP sweep point would otherwise
        // trip CHECK(s.fb[tt] <= 1.0) and abort. Clamping rather than skipping keeps every
        // sweep point usable; the plot shows where the clamp bites as a flattening at large
        // scale. At the production scale the clamp is now far from binding.
        //
        // Note on the stencil: in the middle bin bb = {-0.07, +0.07} is a central difference,
        // but the outer bins are {+0.07, +0.15} and {-0.07, -0.15} -- two forward (or two
        // backward) differences, i.e. the same first-order bias as sig2 (Bulge.h). fb's
        // accuracy therefore depends on which bin it lands in. kFDStepScale shrinks the steps
        // enough that the residual bias is negligible (the sweep already showed fb0/fb1 flat to
        // <0.3% even unscaled), so this is left as-is rather than restructured here; recorded
        // in OPEN_ITEMS.md.
        {
            const double fscale = kFDStepScale * co.deltaScale[(tt == 0) ? 2 : 7];
            for (int b = 0; b < 2; ++b) {
                double step = co.bb[b] * fscale;
                const double lo = 1.0e-6 - s.fb[tt];        //keeps fb strictly above 0
                const double hi = 1.0 - 1.0e-6 - s.fb[tt];  //keeps fb strictly below 1
                if (step < lo) step = lo;
                if (step > hi) step = hi;
                co.bb[b] = step;
            }
        }

        // Step H3b. The unperturbed model magnitude for THIS epoch, recomputed under the
        // observer configuration currently in force -- NOT the cached l.magn[i].
        //
        // Why this is physics and not tidying. Every derivative below is formed as
        //     (model(theta + Delta) - reference) / Delta
        // and for that to BE a derivative the reference must be model(theta) under the SAME
        // observing geometry. l.magn[i] was written when the light curve was generated, at
        // whatever as.satScale the run itself used. While satScale never changes the two agree
        // to the last bit and this line is a no-op -- that is the production case, and it is
        // why the Fisher fixture is unchanged by this.
        //
        // Step H3 breaks that assumption on purpose: it re-characterises one event with Roman
        // moved to Earth (satScale = 0) to ask what that geometry alone could measure. Then
        // model(theta) != l.magn[i], and their difference dm_sat -- the satellite-parallax
        // perturbation itself, the very signal being looked for -- leaks into every derivative
        // as a CONSTANT dm_sat/Delta. With Delta ~ 1e-4 of the parameter that constant is
        // enormous, and it is information the data do not contain.
        //
        // It cancels wherever the stencil is symmetric (sig = {+1,-1}: the +Delta and -Delta
        // terms carry +dm_sat/Delta and -dm_sat/Delta), which is what hid it for u0, tE, piE,
        // xi and t0. It does NOT cancel in the places that decide the answer:
        //   - the cross-telescope rows, where co.diff = 1.0 is a dummy and the derivative must
        //     be exactly ZERO -- Rubin's data says nothing about Roman's blend fraction -- but
        //     becomes dm_sat instead, coupling fb0/mbs0 to Roman epochs and fb1/mbs1 to Rubin
        //     ones and destroying the block structure the partitioned matrices rely on;
        //   - fb's outer bins, where co.bb is two same-signed forward steps rather than a
        //     central pair, so the constant survives multiplied by ~1e5.
        // The astrometric block below is worse still: tetE and piE use sig2 = {+0.5,+1.0}, two
        // forward differences, where nothing cancels for any parameter at all.
        //
        // This explains the measured failure exactly, which is why it is the cause and not a
        // candidate: the satScale = 0 matrix gained spurious information that GREW with the
        // observer separation, so its sigma FELL as du_sat rose (the wrong-sign correlation),
        // sigma_tE looked five-fold better for it, and events with no Roman epochs near the
        // peak -- where dm_sat is identically zero -- returned ratio 1.000000 and made the
        // plumbing look verified. DEVIATIONS.md 36.
        lightcurve(s, l, as, l.timn[i], tt);
        s.Astar = (s.ut * s.ut + 2.0) / std::sqrt(s.ut * s.ut * (s.ut * s.ut + 4.0));
        const double magn0 = s.mbs[tt] - 2.5 * std::log10(s.Astar * s.fb[tt] + 1.0 - s.fb[tt]);

        for (int j = 0; j < Nx; ++j) {
            for (int h = 0; h < 2; ++h) {

                if (j == 0) {co.diff = double(+co.Delta1[j] * sig[h]) ;      l.u0 += co.diff;}
                if (j == 1) {co.diff = double(+co.Delta1[j] * sig[h]) ;      l.tE += co.diff;}
                // fb0 is RUBIN's source-flux fraction. It must perturb s.fb[0] on Rubin epochs
                // only -- never s.fb[tt]. Perturbing s.fb[tt] here (correct before Step C2b, when
                // index 2 was a single telescope-selected fb) makes parameters 2 and 7 both move
                // s.fb[1] on Roman epochs, i.e. exactly degenerate, which corrupts the joint
                // matrix with a duplicated direction.
                if (j == 2) {co.diff = (tt == 0) ? double(co.bb[h]) : 1.0;
                             if (tt == 0) s.fb[0] += co.diff;}
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
                if (j == 3) {co.diff = double(+co.Delta1[j] * sig[h]) ;     l.piE += co.diff;}
                if (j == 4) {co.diff = double(+co.Delta1[j] * sig[h]) ;      s.xi += co.diff;}
                if (j == 5) {co.diff = double(+co.Delta1[j] * sig[h]) ;      l.t0 += co.diff;}

                // Step H1: the observer that produced THIS datum. A derivative evaluated
                // with a different observer than its datum makes the Fisher matrix
                // inconsistent, and the error it produces is not a forecast of anything.
                lightcurve(s, l, as, l.timn[i], int(l.tele[i]));
                s.Astar = (s.ut * s.ut + 2.0) / std::sqrt(s.ut * s.ut * (s.ut * s.ut + 4.0));
                co.magw = s.mbs[tt] - 2.5 * std::log10(s.Astar * s.fb[tt] + 1.0 - s.fb[tt]);
                co.derm1[h] = double(co.magw - magn0) / co.diff;

                CHECK(l.tE > 0.0);
                CHECK(s.fb[tt] > 0.0);
                CHECK(s.fb[tt] <= 1.0);
                CHECK(l.piE > 0.0);
                CHECK(l.u0 != 0.0);
                CHECK(co.diff != 0.0);
                CHECK(l.errm[i] > 0.0);

                if (j == 0)      l.u0 -= co.diff;
                if (j == 1)      l.tE -= co.diff;
                if (j == 2 and tt == 0)  s.fb[0] -= co.diff;
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
                    if (k == 1) {co.diff = double(+co.Delta1[k] * sig[h]) ;     l.tE += co.diff;}
                    if (k == 2) {co.diff = (tt == 0) ? double(co.bb[h]) : 1.0;
                                 if (tt == 0) s.fb[0] += co.diff;}
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
                    if (k == 3) {co.diff = double(+co.Delta1[k] * sig[h]) ;    l.piE += co.diff;}
                    if (k == 4) {co.diff = double(+co.Delta1[k] * sig[h]) ;     s.xi += co.diff;}
                    if (k == 5) {co.diff = double(+co.Delta1[k] * sig[h]) ;     l.t0 += co.diff;}

                    lightcurve(s, l, as, l.timn[i], int(l.tele[i])); //Step H1: same observer as the datum
                    s.Astar = (s.ut * s.ut + 2.0) / std::sqrt(s.ut * s.ut * (s.ut * s.ut + 4.0));
                    co.magw = s.mbs[tt] - 2.5 * std::log10(s.Astar * s.fb[tt] + 1.0 - s.fb[tt]);
                    co.derm2[h] = double(co.magw - magn0) / co.diff;

                    CHECK(l.tE > 0.0);
                    CHECK(s.fb[tt] >= 0.0);
                    CHECK(s.fb[tt] <= 1.0);
                    CHECK(l.piE > 0.0);
                    CHECK(l.u0 != 0.0);
                    CHECK(co.diff != 0.0);
                    CHECK(l.errm[i] > 0.0);

                    if (k == 0)      l.u0 -= co.diff;
                    if (k == 1)      l.tE -= co.diff;
                    if (k == 2 and tt == 0)  s.fb[0] -= co.diff;
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

    co.Delta2[0] = l.tetE * 0.2500467345692345   * kFDStepScaleB[0] * co.deltaScaleB[0];
    co.Delta2[1] = s.mus1 * 0.2500465923465493   * kFDStepScaleB[1] * co.deltaScaleB[1];
    co.Delta2[2] = s.mus2 * 0.25000465936443     * kFDStepScaleB[2] * co.deltaScaleB[2];
    co.Delta2[3] = 0.25000375423543264 * l.piE   * kFDStepScaleB[3] * co.deltaScaleB[3];
    for (int q = 0; q < NSURV; ++q)
    for(int j = 0; j < Ny; ++j){
        for(int k = 0; k < Ny; ++k){
            gsl_matrix_set(co.inputB[q].get(), j, k, 0.0);
            gsl_matrix_set(co.inverB[q].get(), j, k, 0.0);
        }
    }

    // ---- Step 3c (Deviation 71): per-coordinate weights, free reference positions, and the
    // three noise variants (AST_SIGC in config/parameters.h), from ONE pass over the epochs. ----
    //
    // Per epoch the 4-vector of derivatives of each sky coordinate is computed once (Ny x 2 model
    // evaluations). Until Deviation 71 the matrix was built element by element and re-evaluated
    // the k-derivative inside the j loop: 28 evaluations per epoch for the same numbers.
    //
    // Sums kept per coordinate c (0 = x, 1 = y): S_w = sum w, S_wd[c][j] = sum w d_cj, S_wdd[c][jk].
    // Rubin's epochs go into one white group. Roman's go into day blocks (the sigma_c correlation
    // unit), each tagged with its season and roll so the variants can group them differently.
    struct AstSums {
        double Sw = 0.0;
        std::array<std::array<double, Ny>, 2>      Swd{};
        std::array<std::array<double, Ny * Ny>, 2> Swdd{};
        int season = -1, roll = -1;
        void add(double w, const std::array<double, Ny>& dx, const std::array<double, Ny>& dy) {
            Sw += w;
            for (int c = 0; c < 2; ++c) {
                const auto& d = c ? dy : dx;
                for (int j = 0; j < Ny; ++j) {
                    Swd[c][j] += w * d[j];
                    for (int k = 0; k < Ny; ++k) Swdd[c][j * Ny + k] += w * d[j] * d[k];
                }
            }
        }
    };
    AstSums rubin;
    std::map<long, AstSums> romanDays;
    std::array<double, Ny> dx{}, dy{};

    for (int i = 0; i < ndw; ++i) {
        // Step H3b / H1: reference and every perturbation evaluated from the observer that
        // produced THIS datum (see the history of this loop in DEVIATIONS 27 and 37).
        lightcurve(s, l, as, l.timn[i], int(l.tele[i]));
        const double soux0 = s.pos1c, souy0 = s.pos2c;

        for (int j = 0; j < Ny; ++j) {
            for (int h = 0; h < 2; ++h) {
                // Deviation 75: all four astrometric parameters use the CENTRAL stencil (tetE and
                // piE used sig2, two forward differences with an O(h) bias; Step C3 fixed the same
                // thing on the photometric side).
                co.diff = double(co.Delta2[j] * sig[h]);
                if (j == 0) l.tetE += co.diff;
                if (j == 1) s.mus1 += co.diff;
                if (j == 2) s.mus2 += co.diff;
                if (j == 3) l.piE  += co.diff;

                lightcurve(s, l, as, l.timn[i], int(l.tele[i]));
                co.dera1[h] = double(s.pos1c - soux0) / co.diff;
                co.derb1[h] = double(s.pos2c - souy0) / co.diff;

                CHECK(l.tetE > 0.0);
                CHECK(l.piE > 0.0);
                CHECK(co.diff != 0.0);
                // A null derivative is legitimate (an epoch at t0 or t = 0, or far from the peak);
                // it contributes nothing to that parameter, which is right. See OPEN_ITEMS.md.

                if (j==0) l.tetE -= co.diff;
                if (j==1) s.mus1 -= co.diff;
                if (j==2) s.mus2 -= co.diff;
                if (j==3)  l.piE -= co.diff;
            }
            dx[j] = (co.dera1[0] + co.dera1[1]) * 0.5;
            dy[j] = (co.derb1[0] + co.derb1[1]) * 0.5;
        }
        CHECK(l.erra[i] > 0.0);
        const double w = 1.0 / (l.erra[i] * l.erra[i]);       // per coordinate (Deviation 71)
        if (int(l.tele[i]) == 0) {
            rubin.add(w, dx, dy);
        } else {
            CHECK(l.rseas[i] >= 0 and l.rseas[i] < AST_MAX_SEASONS);
            CHECK(l.rroll[i] == 0 or l.rroll[i] == 1);
            AstSums& d = romanDays[long(std::floor(l.timn[i]))];
            if (d.season < 0) { d.season = l.rseas[i]; d.roll = l.rroll[i]; }
            CHECK(d.season == l.rseas[i] and d.roll == l.rroll[i]);
            d.add(w, dx, dy);
        }
    }

    // A frame group: sum of block matrices, the block offset-vectors b per coordinate, and c.
    struct Group {
        std::array<double, Ny * Ny> F{};
        std::array<std::array<double, Ny>, 2> b{};
        double c = 0.0;
        void fold(const AstSums& k, double sc) {          // one block, correlated at sigma_c = sc
            const double s2 = sc * sc, den = 1.0 + s2 * k.Sw;
            for (int cc = 0; cc < 2; ++cc)
                for (int j = 0; j < Ny; ++j) {
                    b[cc][j] += k.Swd[cc][j] / den;
                    for (int m = 0; m < Ny; ++m)
                        F[j * Ny + m] += k.Swdd[cc][j * Ny + m] - s2 * k.Swd[cc][j] * k.Swd[cc][m] / den;
                }
            c += k.Sw / den;
        }
        void marginaliseInto(std::array<double, Ny * Ny>& out) const {   // free offset
            for (int j = 0; j < Ny; ++j)
                for (int m = 0; m < Ny; ++m)
                    out[j * Ny + m] += F[j * Ny + m]
                        - (c > 0.0 ? (b[0][j] * b[0][m] + b[1][j] * b[1][m]) / c : 0.0);
        }
    };

    std::array<double, Ny * Ny> FL{};
    { Group g; g.fold(rubin, 0.0); g.marginaliseInto(FL); }
    for (int v = 0; v < NAVAR; ++v) {
        std::vector<Group> groups(v == AV_W ? 1 : (v == AV_N ? 2 : AST_MAX_SEASONS));
        for (const auto& [day, k] : romanDays)
            groups[v == AV_W ? 0 : (v == AV_N ? k.roll : k.season)].fold(k, AST_SIGC[v]);
        std::array<double, Ny * Ny> FR{};
        for (const auto& g : groups) g.marginaliseInto(FR);
        for (int jm = 0; jm < Ny * Ny; ++jm) {
            co.FBV[v][SRUBIN][jm] = FL[jm];
            co.FBV[v][SROMAN][jm] = FR[jm];
            co.FBV[v][SJOINT][jm] = FL[jm] + FR[jm];
        }
    }
    // The main matrices are variant W.
    for (int q = 0; q < NSURV; ++q)
        for (int j = 0; j < Ny; ++j)
            for (int m = 0; m < Ny; ++m)
                gsl_matrix_set(co.inputB[q].get(), j, m, co.FBV[AV_W][q][j * Ny + m]);

    for (int q = 0; q < NSURV; ++q) {
        co.okB[q] = (co.nepochA[q] >= Ny) ? invert_matrix(co, 1, q) : 0;
    }

    // Variants N and P, through the same normalised inversion and condition cut.
    {
        static const std::vector<int> kAllAst = {0, 1, 2, 3};
        gsl_matrix* Fin  = gsl_matrix_alloc(Ny, Ny);
        gsl_matrix* Fout = gsl_matrix_alloc(Ny, Ny);
        for (int v = 0; v < NAVAR; ++v)
            for (int q = 0; q < NSURV; ++q) {
                co.okBV[v][q] = 0; co.condBV[v][q] = -1.0;
                for (int k = 0; k < Ny; ++k) co.ErbV[v][q][k] = -1.0;
                if (v == AV_W or co.nepochA[q] < Ny) continue;      // W is filled in ErrorCal
                for (int j = 0; j < Ny; ++j)
                    for (int m = 0; m < Ny; ++m) gsl_matrix_set(Fin, j, m, co.FBV[v][q][j * Ny + m]);
                co.okBV[v][q] = invertNormalized(Fin, Fout, kAllAst, co.condBV[v][q], nullptr);
                if (co.okBV[v][q])
                    for (int k = 0; k < Ny; ++k)
                        co.ErbV[v][q][k] = std::sqrt(std::fabs(gsl_matrix_get(Fout, k, k)));
            }
        gsl_matrix_free(Fin);
        gsl_matrix_free(Fout);
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

  // Fractional 1-sigma on the lens mass, per survey (Step D1).
  //
  // Ml is never fitted. It follows from the two Einstein-radius observables,
  //     Ml = tetE / (kappa * piE),   kappa = 8.144 mas/Msun,
  // and because that is a pure ratio the fractional errors add in quadrature. Computed
  // per survey because its two ingredients come from different instruments' strengths:
  // tetE from the astrometric matrix (sub-mas centroid motion -- Roman) and piE from the
  // photometric one over a long time baseline (annual parallax distortion -- Rubin). A
  // mass the joint fit measures and neither survey measures alone is the whole point.
  for (int q = 0; q < NSURV; ++q) {
      // Sentinel-aware throughout: -1.0 means "not measured", never "measured to be -1".
      // piE has two independent routes -- the photometric matrix (Era[3]) and the
      // astrometric one (Erb[3]) -- so take the better of whichever is actually
      // available, and refuse to report a mass at all if either ingredient is missing.
      const double rp1 = (co.Era[q][3] >= 0.0) ? co.Era[q][3] / (std::fabs(l.piE)  + eps) : -1.0;
      const double rp2 = (co.Erb[q][3] >= 0.0) ? co.Erb[q][3] / (std::fabs(l.piE)  + eps) : -1.0;
      const double rp  = (rp1 >= 0.0 and rp2 >= 0.0) ? MIN(rp1, rp2) : std::max(rp1, rp2);
      const double rt  = (co.Erb[q][0] >= 0.0) ? co.Erb[q][0] / (std::fabs(l.tetE) + eps) : -1.0;
      co.relMl[q] = (rp >= 0.0 and rt >= 0.0) ? std::sqrt(rp * rp + rt * rt) : -1.0;
  }

  // The astrometric noise variants (Deviation 71): W is the main result; N and P were inverted
  // in FisherM. Their lens-mass errors follow the same rule, with each variant's own tetE and
  // astrometric piE and the (variant-independent) photometric piE.
  for (int q = 0; q < NSURV; ++q) {
      co.okBV[AV_W][q] = co.okB[q];
      co.condBV[AV_W][q] = co.condB[q];
      for (int k = 0; k < Ny; ++k) co.ErbV[AV_W][q][k] = co.Erb[q][k];
      for (int v = 0; v < NAVAR; ++v) {
          const auto& E = co.ErbV[v][q];
          const double rp1 = (co.Era[q][3] >= 0.0) ? co.Era[q][3] / (std::fabs(l.piE) + eps) : -1.0;
          const double rp2 = (E[3] >= 0.0) ? E[3] / (std::fabs(l.piE) + eps) : -1.0;
          const double rp  = (rp1 >= 0.0 and rp2 >= 0.0) ? MIN(rp1, rp2) : std::max(rp1, rp2);
          const double rt  = (E[0] >= 0.0) ? E[0] / (std::fabs(l.tetE) + eps) : -1.0;
          co.relMlV[v][q] = (rp >= 0.0 and rt >= 0.0) ? std::sqrt(rp * rp + rt * rt) : -1.0;
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

  // piE is measured twice over: photometrically (resu[3], from Era[3]) and astrometrically
  // (resu[8], from Erb[3]). Take the better -- but ONLY among the routes that actually
  // produced a measurement. ErrorCal writes -1.0 as an explicit "not characterizable"
  // sentinel, and an unguarded MIN picks that -1 in preference to a real sigma whenever
  // the astrometric matrix is singular. That is not hypothetical: it is the source of the
  // resu[3] = -697 event recorded in DEVIATIONS entry 8, which was worked around by
  // dropping such events from the field average rather than fixed. It propagates into
  // resu[9] (mass) and resu[10] (distance), so both were wrong on those events too.
  // Same rule as relMl[] above; the two must agree, since they are the same quantity.
  if      (co.resu[3] < 0.0 and co.resu[8] < 0.0) co.resu[3] = -1.0;
  else if (co.resu[3] < 0.0)                      co.resu[3] = co.resu[8];
  else if (co.resu[8] >= 0.0)                     co.resu[3] = MIN(co.resu[3], co.resu[8]);

  //Lens Mass: Ml = tetE/(kappa*piE), a pure ratio, so the fractional errors add in
  //quadrature. -1 if either ingredient is missing -- see relMl[] above.
  co.resu[9]  = (co.resu[3] >= 0.0 and co.resu[5] >= 0.0)
              ? std::sqrt(co.resu[3] * co.resu[3] + co.resu[5] * co.resu[5]) : -1.0;
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
