// Inversion of the Fisher information matrices (normalised, with a condition-number cut).
#include "fisher/linalg.h"
#include "fisher/fisher.h"

// Largest condition number of the NORMALIZED matrix from which an inverse is still trusted.
// Inverting at condition 10^k costs about k of double's ~16 digits, and the few-percent precision
// gains being measured need 3-4 good digits, so 1e12 is the limit. The per-event condition number
// is stored regardless, so a different cut can be applied in analysis.
constexpr double kMaxCondition = 1.0e12;

int invert_matrix(covarian & co, int flag, int surv)
{
    // Returns 1 on a usable inverse, 0 if the matrix is singular or too ill-conditioned (the
    // partition cannot characterize the event). The photometric matrix is inverted only over
    // activePhotParams; the inverse is scattered back to full size with inactive rows and columns
    // at zero, and ErrorCal reports sigma = -1 for those.
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

// The core of invert_matrix, usable on any information matrix (the astrometric noise variants use
// it too). `in` is read; `out` receives the inverse scattered to full size (zero elsewhere).
int invertNormalized(const gsl_matrix* in, gsl_matrix* out, const std::vector<int>& act,
                     double& cond, double* deter)
{
    const int dim = static_cast<int>(act.size());

    cond = -1.0;
    gsl_matrix_set_zero(out);

    // ---- 1. Normalizing scale D = diag(1/sqrt(F_ii)) over the active parameters ----
    // F_jk has units 1/(theta_j theta_k), so entries span many orders of magnitude from the choice
    // of units alone (tE ~30 d, u0 ~0.3, xi in rad, mbs in mag). A non-positive diagonal means
    // the parameter has no information.
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

    // ---- 3. Condition number of the normalized matrix (lambda_max/lambda_min, symmetric PSD) ----
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

    // ---- 5. Undo the scaling and scatter to full size: F^-1 = D Ftilde^-1 D (exact) ----
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
