#ifndef CGHSEG_COMPACTEM_ALGO_H
#define CGHSEG_COMPACTEM_ALGO_H

#include "numlib.h"

#include <iostream>
//#include <fstream>
#define Min(i, j)  ((i) < (j) ? (i) : (j))
#define Max(i, j)  ((i) > (j) ? (i) : (j))


namespace cghseg
{

  void compute_compactEM_algo(numlib_vector *xkbar,numlib_vector *x2kbar, numlib_matrix *nk, int P, numlib_vector *phistart, bool vh, numlib_vector **phiend, numlib_matrix **tau, double *lvinc_pointer, int *empty_pointer, int *dv_pointer);

  class compactEM_algo
  {
  public:
    double    _lengthx;        // size of the data
    int       _K;              // Number of segments 
    int       _P;              // Number of clusters
    double   *_phi;            // vector of parameters
    double    _lvinc ;         // Incomplete data loglikelihood
    int       _empty;          // =1 if empty clusters
    int       _dv;             // =2 if EM did not converge
    double  **_tau;            // matrix of posterior proba
    double   *_xkbar;          // empirical means
    double   *_x2kbar;         // empirical means of squares
    double   *_nk;             // length of segments
    double   *_wk;             // within variances for segments
    bool      _vh;             // variance homogeneity;
    
    compactEM_algo(int nbseg, int nbclust, bool varh);
    void Estep();
    void Mstep();
    void compactEM();
    void Init(double *Datak, double *Data2k,double *param, double *Datank);
    friend std::ostream & operator << (std::ostream &s, const compactEM_algo & compactEMa);
    ~compactEM_algo();
  };
}
#endif

  

