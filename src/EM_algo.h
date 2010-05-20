#ifndef CGHSEG_EM_ALGO_H
#define CGHSEG_EM_ALGO_H

#include "numlib.h"

#include <iostream>
//#include <fstream>
#define Min(i, j)  ((i) < (j) ? (i) : (j))
#define Max(i, j)  ((i) > (j) ? (i) : (j))


namespace cghseg
{

  void compute_EM_algo(numlib_vector *x, numlib_matrix *rupt, int P, numlib_vector *phistart, bool vh, numlib_vector **phiend, numlib_matrix **tau, double *lvinc_pointer, int *empty_pointer, int *dv_pointer);
  
  class EM_algo
  {
  public:
    int       _lengthx;        // size of the data
    int       _K;              // Number of segments 
    int       _P;              // Number of clusters
    double   *_phi;            // vector of parameters
    double    _lvinc ;         // Incomplete data loglikelihood
    int       _empty;          // =1 if empty clusters
    int       _dv;             // =2 if EM did not converge
    double  **_tau;            // matrix of posterior proba
    double   *_x;              // data
    double   *_xkbar;          // empirical means
    double   *_x2kbar;         // empirical means of squares
    double   *_wk;             // within variances for segments
    int      *_Breaks ;        // breakpoints vector of size K+1
    bool      _vh;             // variance homogeneity;
    
    EM_algo(int lengthdata, int nbseg, int nbclust, bool varh);
    void Estep();
    void Mstep();
    void EM();
    void Init(double *Data, double *param, int *RuptVect);
    friend std::ostream & operator << (std::ostream &s, const EM_algo & EMa);
    ~EM_algo();
  };
}
#endif

  

