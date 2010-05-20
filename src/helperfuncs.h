#ifndef CGHSEG_HELPERFUNCS_H
#define CGHSEG_HELPERFUNCS_H

#include "numlib.h"

namespace cghseg {

typedef double (*basicfunc)(const double);

numlib_matrix * apply_basicfunc_numlib_matrix(numlib_matrix *, basicfunc);

numlib_vector * apply_basicfunc_numlib_vector(numlib_vector *, basicfunc);

numlib_vector *copy_range_numlib_vector(numlib_vector *, int, int);

numlib_vector *build_range_numlib_vector(int, int);


numlib_vector * rowsum_numlib_matrix(numlib_matrix *);

numlib_vector * colsum_numlib_matrix(numlib_matrix *);

double sum_numlib_matrix(numlib_matrix *);

double sum_numlib_vector(numlib_vector *);

numlib_vector * cumsum_numlib_matrix(numlib_matrix *);

numlib_vector * cumsum_numlib_vector(numlib_vector *);

numlib_vector *colmax_numlib_matrix(numlib_matrix *);

numlib_vector *rowmax_numlib_matrix(numlib_matrix *);

numlib_vector *colmin_numlib_matrix(numlib_matrix *);

numlib_vector *rowmin_numlib_matrix(numlib_matrix *);


double mean_numlib_vector(numlib_vector *);

double var_numlib_vector(numlib_vector *);

double stddev_numlib_vector(numlib_vector *);

numlib_vector *diff_numlib_vector(numlib_vector *);

int sum_diffs_numlib_vectors(numlib_vector *, numlib_vector *);

int sum_diffs_numlib_matrices(numlib_matrix *, numlib_matrix *);

double checked_max_numlib_vector(numlib_vector *);

int checked_max_index_numlib_vector(numlib_vector *);


double checked_min_numlib_vector(numlib_vector *);

int checked_min_index_numlib_vector(numlib_vector *);

}
#endif
