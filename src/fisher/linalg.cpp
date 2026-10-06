// Inversion of the Fisher information matrices (normalised, with a condition-number cut).
#include "fisher/linalg.h"
#include "fisher/fisher.h"

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
    return invertNormalized(in, out, act, cond, &co.deter);
}

// The core of invert_matrix, on any information matrix (Deviation 71: the astrometric noise
// variants are inverted through exactly the same normalisation and condition cut as the main
// matrices). `in` is read, `out` receives the inverse scattered to full size (zero elsewhere).
int invertNormalized(const gsl_matrix* in, gsl_matrix* out, const std::vector<int>& act,
                     double& cond, double* deter)
{
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
    if (deter) *deter = det;

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
