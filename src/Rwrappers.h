#ifndef CGHSEG_R_WRAPPERS_H
#define CGHSEG_R_WRAPPERS_H

#include <R.h>
#include <Rinternals.h>
#define Min(i, j)  ((i) < (j) ? (i) : (j))
#define Max(i, j)  ((i) > (j) ? (i) : (j))
extern "C" {
  SEXP sc_estep(SEXP logdensity, SEXP rowslogdensity, SEXP colslogdensity, SEXP phi);
  SEXP sc_hybrid(SEXP x, SEXP P, SEXP Kmax, SEXP lmin, SEXP lmax, SEXP vh, SEXP fast);
  SEXP sc_logdens(SEXP xk, SEXP P, SEXP phi); 
  SEXP sc_segmean(SEXP x, SEXP lmin, SEXP lmax, SEXP Kmax, SEXP vh);
  SEXP sc_segmixt(SEXP x, SEXP lmin, SEXP lmax, SEXP Kmax, SEXP phi, SEXP P);
  SEXP sc_segibp(SEXP ContrastR, SEXP KseqR, SEXP multiKmaxR);
  SEXP sc_EMalgo(SEXP xR, SEXP phiR, SEXP ruptR, SEXP KR, SEXP PR, SEXP vhR);
  SEXP sc_EMinit(SEXP xR, SEXP ruptR, SEXP KR, SEXP PR, SEXP vhR);
}

#endif
