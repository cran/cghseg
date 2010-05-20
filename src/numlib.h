#ifndef CGHSEG_CGHSEG_NUMLIB_H
#define CGHSEG_CGHSEG_NUMLIB_H

#ifdef USE_GSL

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_int.h>

typedef gsl_matrix numlib_matrix;
typedef gsl_matrix_view numlib_matrix_view;
typedef gsl_permutation numlib_permutation;
typedef gsl_vector numlib_vector;
typedef gsl_vector_int numlib_vector_int;

#define NUMLIB_MAX GSL_MAX
#define NUMLIB_MIN GSL_MIN
#define NUMLIB_NAN GSL_NAN
#define NUMLIB_NEGINF GSL_NEGINF
#define NUMLIB_POSINF GSL_POSINF

#define numlib_blas_ddot(v1,v2,r) (gsl_blas_ddot((v1),(v2),(r)))

#define numlib_finite(v) (gsl_finite((v)))
#define numlib_isinf(v) (gsl_isinf((v)))
#define numlib_isnan(v) (gsl_isnan((v)))

#define numlib_matrix_add(m1,m2) (gsl_matrix_add((m1),(m2)))
#define numlib_matrix_alloc(size1,size2) (gsl_matrix_alloc((size1),(size2)))
#define numlib_matrix_calloc(size1,size2) (gsl_matrix_alloc((size1),(size2)))
#define numlib_matrix_div_elements(m1,m2) (gsl_matrix_div_elements((m1),(m2)))
#define numlib_matrix_free(m) (gsl_matrix_free((m)))
#define numlib_matrix_get(m,i,j) (gsl_matrix_get((m),(i),(j)))
#define numlib_matrix_get_col(c,m,i) (gsl_matrix_get_col((c),(m),(i)))
#define numlib_matrix_get_row(r,m,i) (gsl_matrix_get_row((r),(m),(i)))
#define numlib_matrix_memcpy(m1,m2) (gsl_matrix_memcpy((m1),(m2)))
#define numlib_matrix_set(m,i1,i2,v) (gsl_matrix_set((m),(i1),(i2),(v)))
#define numlib_matrix_set_all(m,v) (gsl_matrix_set_all((m),(v)))
#define numlib_matrix_set_col(m,i,v) (gsl_matrix_set_col((m),(i),(v)))
#define numlib_matrix_set_row(m,i,v) (gsl_matrix_set_row((m),(i),(v)))
#define numlib_matrix_sub(m1,m2) (gsl_matrix_sub((m1),(m2)))
#define numlib_matrix_view_array(a,size1,size2) (gsl_matrix_view_array((a),(size1),(size2)))

#define numlib_permutation_calloc(n) (gsl_permutation_calloc((n)))
#define numlib_permutation_free(p) (gsl_permutation_free((p)))
#define numlib_permutation_get(p,i) (gsl_permutation_get((p),(i)))

#define numlib_sort_vector_index(p,v) (gsl_sort_vector_index((p),(v)))
#define numlib_sort_vector(v) (gsl_sort_vector((v)))

#define numlib_vector_add(v1,v2) (gsl_vector_add((v1),(v2)))
#define numlib_vector_add_constant(v,c) (gsl_vector_add_constant((v),(c)))
#define numlib_vector_alloc(n) (gsl_vector_alloc((n)))
#define numlib_vector_calloc(n) (gsl_vector_calloc((n)))
#define numlib_vector_div(v1,v2) (gsl_vector_div((v1),(v2)))
#define numlib_vector_free(v) (gsl_vector_free((v)))
#define numlib_vector_get(v,i) (gsl_vector_get((v),(i)))
#define numlib_vector_memcpy(v1,v2) (gsl_vector_memcpy((v1),(v2)))
#define numlib_vector_mul(v1,v2) (gsl_vector_mul((v1),(v2)))
#define numlib_vector_scale(v,c) (gsl_vector_scale((v),(c)))
#define numlib_vector_set(v,i,x) (gsl_vector_set((v),(i),(x)))
#define numlib_vector_set_all(v,x) (gsl_vector_set_all((v),(x)))
#define numlib_vector_sub(v1,v2) (gsl_vector_sub((v1),(v2)))

#define numlib_vector_int_get(v,i) (gsl_vector_int_get((v),(i)))

#define numlib_sf_exp(x) (gsl_sf_exp((x)))
#define numlib_sf_log(x) (gsl_sf_log((x)))
#define numlib_pow_2(x) (gsl_sf_pow_2((x)))

#else

#include <cfloat>
#include <cmath>
#include <limits>


#define M_LNPI     1.14472988584940017414342735135 // Cut & Paste from GSL.

#define NUMLIB_MAX(a,b) ((a)>=(b)?(a):(b))
#define NUMLIB_MIN(a,b) ((a)<=(b)?(a):(b))


// The following constants are defined as their GSL counterparts.
static const double NUMLIB_NAN=(0.0/0.0) ;
static const double NUMLIB_NEGINF=(-1.0/0.0);
static const double NUMLIB_POSINF=(+1.0/0.0);


#define numlib_finite(v) (std::isfinite((v)))
#define numlib_isinf(v) (std::isinf((v)))
#define numlib_isnan(v)(std::isnan((v)))

typedef struct {
  int size1, size2;
  double *data;
} numlib_matrix;

typedef struct {
  numlib_matrix matrix;
} numlib_matrix_view;

typedef struct {
  int size;
  double *data;
} numlib_vector;

typedef struct {
  int size;
  int *data;
} numlib_vector_int;

typedef numlib_vector_int numlib_permutation;



void numlib_blas_ddot(numlib_vector *v1,numlib_vector *v2, double *r);

int numlib_matrix_add(numlib_matrix *m1, numlib_matrix *m2);
numlib_matrix *numlib_matrix_alloc(int size1,int size2);
numlib_matrix *numlib_matrix_calloc(int size1, int size2);
int numlib_matrix_div_elements(numlib_matrix *m1, numlib_matrix *m2);
void numlib_matrix_free(numlib_matrix *m);
double numlib_matrix_get(numlib_matrix *m,int i,int j);
void numlib_matrix_get_col(numlib_vector *c,numlib_matrix *m,int i);
void numlib_matrix_get_row(numlib_vector *r,numlib_matrix *m,int i);
void numlib_matrix_memcpy(numlib_matrix *m1,numlib_matrix *m2);
void numlib_matrix_set(numlib_matrix *m,int i1,int i2,double v);
void numlib_matrix_set_all(numlib_matrix *m,double v);
void numlib_matrix_set_col(numlib_matrix *m,int i,numlib_vector *v);
void numlib_matrix_set_row(numlib_matrix *m,int i,numlib_vector * v);
int numlib_matrix_sub(numlib_matrix *m1,numlib_matrix *m2);
numlib_matrix_view numlib_matrix_view_array(double *a,int size1,int size2);

numlib_permutation *numlib_permutation_calloc(int n);
void  numlib_permutation_free(numlib_permutation *p);
int numlib_permutation_get(numlib_permutation *p,int i);

void  numlib_sort_vector_index(numlib_permutation *p,numlib_vector *v);
void  numlib_sort_vector(numlib_vector *v);

int numlib_vector_add(numlib_vector *v1,numlib_vector *v2);
int numlib_vector_add_constant(numlib_vector *v,double c);
numlib_vector *numlib_vector_alloc(int n);
numlib_vector *numlib_vector_calloc(int n);
int numlib_vector_div(numlib_vector *v1, numlib_vector *v2);
void  numlib_vector_free(numlib_vector *v);
double numlib_vector_get(numlib_vector *v,int i);
void numlib_vector_memcpy(numlib_vector *v1, numlib_vector *v2);
int numlib_vector_mul(numlib_vector *v1,numlib_vector *v2);
int numlib_vector_scale(numlib_vector *v, double c);
void numlib_vector_set(numlib_vector *v, int i, double x);
void numlib_vector_set_all(numlib_vector *v, double x);
int numlib_vector_sub(numlib_vector *v1, numlib_vector *v2);

int numlib_vector_int_get(numlib_vector_int *v,int i);

double numlib_sf_exp(double x);
double numlib_sf_log(double x);
double numlib_pow_2(double x);

#endif

#endif
