// The per-event output table header.
#include "run/outputs.h"

// ---------------------------------------------------------------------------
// Column names of the per-event table, in exactly the order the `filg_in <<` block in
// main() writes them (Step D1).
//
// Written once at the top of the file so a table is self-describing: an unlabelled
// 88-column matrix is unusable six months later, and mis-numbering a column by one is
// the kind of error that produces a plausible plot of the wrong quantity. Anything
// added to the row MUST be appended both here and there, in the same place -- the
// verification step is `head -1 | wc -w` against a data line's `wc -w`.
//
// Note the two indexing systems, which look alike and are not (see CLAUDE.md):
// magb_*/blend_* are per FILTER (ugrizy, F146); mbs0/fb0 and mbs1/fb1 are per TELESCOPE
// (0 = Rubin, 1 = Roman).
// ---------------------------------------------------------------------------
const char* eventTableHeader()
{
    return
        "# icon FFG0 tE RE_AU piE tetE Vt u0 Ml opt_1e6 Dl Ds vl vs mbs0 fb0 gg struc "
        "FWHM_yr vsave_mus DeltaT_errA murel_yr "
        "rel_u0 rel_tE rel_fb0 rel_piE rel_tetE rel_Ml rel_Dl rel_mul rel_mus "
        "Map_r nsbl_r flagi Ai_r "
        "ndw_L ndw_R detL detR detJ okA_J okA_L okA_R "
        "sigtE_J sigtE_L sigtE_R sigpiE_J sigpiE_L sigpiE_R sigtetE_J sigtetE_L sigtetE_R "
        "detCls synClass condA_J condA_L condA_R "
        "t0 xi lon lat mbs1 fb1 "
        "magb_u magb_g magb_r magb_i magb_z magb_y magb_F146 "
        "blend_u blend_g blend_r blend_i blend_z blend_y blend_F146 "
        "relMl_J relMl_L relMl_R okB_J okB_L okB_R condB_J condB_L condB_R "
        "dt_edge t0zone w_area du_sat nepL_pk nepR_pk "
        // Step R1. Epoch counts at which the two lensing-induced images were separately
        // detectable AND separated by more than the bar named in the suffix; dsep_max is the
        // largest separation reached at such an epoch, -1 if there was none.
        "nres5_L nres20_L nresPSF_L dsep_max_L nres5_R nres20_R nresPSF_R dsep_max_R "
        // Deviation 71. The astrometric forecast under the N (nominal) and P (pessimistic)
        // noise variants, joint and Roman partitions; the main sigtetE_*/relMl_*/okB_* columns
        // are variant W (white). Rubin's partition is the same in all three. See AST_SIGC.
        "sigtetE_NJ sigtetE_NR sigtetE_PJ sigtetE_PR relMl_NJ relMl_NR relMl_PJ relMl_PR "
        "okB_NJ okB_NR okB_PJ okB_PR "
        // Deviation 74: luminous lens (1 = a main-sequence star whose light is blended) and the
        // lens's share of the baseline flux in Rubin's reference band and in F146.
        "lensLum fLens_L fLens_R "
        // Deviation 76: the observed (Earth-frame, parallax-bent) peak time and impact parameter.
        // t0zone, dt_edge and nep_pk_* are now measured from t0obs, not from t0.
        "t0obs umin_obs";
}
