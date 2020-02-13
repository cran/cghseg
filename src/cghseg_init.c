#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void ClassiSeg_cc(void *, void *, void *, void *, void *, void *, void *);
extern void ClassiSegProba(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void colibriR_cc(void *, void *, void *, void *, void *, void *, void *);
extern void LinProgDyn(void *, void *, void *, void *, void *);
extern void LinProgDynMelange(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void meanRuptR_cc(void *, void *, void *, void *);
extern void meansqRuptR_cc(void *, void *, void *, void *);

/* .Call calls */
extern SEXP sc_compactEMalgo(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sc_compactEMinit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sc_EMalgo(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sc_EMinit(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sc_estep(SEXP, SEXP, SEXP, SEXP);
extern SEXP sc_hybrid(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sc_logdens(SEXP, SEXP, SEXP);
extern SEXP sc_Lvinc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sc_segibp(SEXP, SEXP, SEXP);
extern SEXP sc_segmean(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sc_segmixt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"ClassiSeg_cc",      (DL_FUNC) &ClassiSeg_cc,      7},
  {"ClassiSegProba",    (DL_FUNC) &ClassiSegProba,    9},
  {"colibriR_cc",       (DL_FUNC) &colibriR_cc,       7},
  {"LinProgDyn",        (DL_FUNC) &LinProgDyn,        5},
  {"LinProgDynMelange", (DL_FUNC) &LinProgDynMelange, 9},
  {"meanRuptR_cc",      (DL_FUNC) &meanRuptR_cc,      4},
  {"meansqRuptR_cc",    (DL_FUNC) &meansqRuptR_cc,    4},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"sc_compactEMalgo", (DL_FUNC) &sc_compactEMalgo, 7},
  {"sc_compactEMinit", (DL_FUNC) &sc_compactEMinit, 7},
  {"sc_EMalgo",        (DL_FUNC) &sc_EMalgo,        6},
  {"sc_EMinit",        (DL_FUNC) &sc_EMinit,        5},
  {"sc_estep",         (DL_FUNC) &sc_estep,         4},
  {"sc_hybrid",        (DL_FUNC) &sc_hybrid,        7},
  {"sc_logdens",       (DL_FUNC) &sc_logdens,       3},
  {"sc_Lvinc",         (DL_FUNC) &sc_Lvinc,         7},
  {"sc_segibp",        (DL_FUNC) &sc_segibp,        3},
  {"sc_segmean",       (DL_FUNC) &sc_segmean,       5},
  {"sc_segmixt",       (DL_FUNC) &sc_segmixt,       6},
  {NULL, NULL, 0}
};

void R_init_cghseg(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
