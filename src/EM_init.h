#ifndef CGHSEG_EMINIT_H
#define CGHSEG_EMINIT_H

#include "numlib.h"

#include <iostream>
#include <fstream>

namespace cghseg
{
  void compute_EM_init(numlib_vector *x, numlib_matrix *rupt, int P , numlib_vector **phiend);

  class EM_init{
    
  public:
    int       _lengthx;        // size of the data
    int       _K;              // Number of segments 
    int       _P;              // Number of clusters
    double   *_phi;            // parameters
    double   *_x;              // data
    double   *_mk;             // empirical means
    double   *_vk;             // empirical var
    int      *_nk;             // length of segments
    double   *_mk0;             // empirical means
    double   *_vk0;             // empirical var
    int      *_nk0;             // length of segments
    int      *_Breaks ;        // breakpoints
    double   **_D;             // distance matrix
    double   *_Dtmp;           // distance vector of the merged groups
    double   **_tau;           // posterior
    
    EM_init(int lengthdata, int nbseg, int nbclust);
    void CAH();
    void compute_phi();
    void Init(double *Data, int *RuptVect);    
    friend std::ostream & operator << (std::ostream &s, const EM_init & EMi);
    ~EM_init();
  };
}
#endif



