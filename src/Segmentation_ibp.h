#ifndef CGHSEG_SEG_IBP_H
#define CGHSEG_SEG_IBP_H

#include "numlib.h"

#include <iostream>
//#include <fstream>
#define Min(i, j)  ((i) < (j) ? (i) : (j))
#define Max(i, j)  ((i) > (j) ? (i) : (j))




namespace cghseg
{

  void compute_segmentation_ibp(numlib_vector *Jvector, numlib_vector_int *Kvector, int multiK, numlib_vector **J_est, numlib_matrix **t_est);


  class Segmentation_ibp
  {
  public:
    int      _M;           // number of series
    int      _multiKmax;   // max number of segments for all series
    int     *_Kseq;        // number of segments per series
    double **_Contrast;   // contrasts
    double **_J;           // matrix of the lengths of the minimum distances
    int    **_Breaks;      // possible breakpoints
    int    **_BestBreaks;  // best breakpoints
    Segmentation_ibp(int nbseries, int *seqk, int mKmax);
    void DynProg(int m);
    void Terminate(int Ktot);
    void Init(double **cont);
    friend std::ostream & operator << (std::ostream &s, const Segmentation_ibp & Seg_ibp);
    ~Segmentation_ibp();
  };
}
#endif
