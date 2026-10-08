// GSL matrix helpers: an owning pointer, the normalised/condition-checked inversion, and a debug print.
#ifndef ROMAN_FISHER_LINALG_H
#define ROMAN_FISHER_LINALG_H

#include "common.h"

///============================================================================
struct GSLMatrixDeleter {
        void operator()(gsl_matrix* m) const {
        gsl_matrix_free(m);
    }
};

using gsl_matrix_uptr = std::unique_ptr<gsl_matrix,GSLMatrixDeleter>;

struct covarian;

// Inverts one of the Fisher matrices in place. `flag` selects photometric (0, Nx) or
// astrometric (1, Ny); `surv` selects which SurveyIdx partition. Returns 1 if the matrix was
// non-singular and the inverse is usable, 0 if it was singular -- in which case the caller must
// treat that partition as not-characterizable rather than reading numbers out of it.
int    invert_matrix(covarian & co, int flag, int surv);
int    invertNormalized(const gsl_matrix* in, gsl_matrix* out, const std::vector<int>& act,
                        double& cond, double* deter);   // invert_matrix's core
void   print_mat_contents(gsl_matrix *matrix, int);

#endif // ROMAN_FISHER_LINALG_H
