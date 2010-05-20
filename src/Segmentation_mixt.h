#ifndef CGHSEG_SEG_MIXT_H
#define CGHSEG_SEG_MIXT_H

#include "numlib.h"

#include <iostream>
/**#include <fstream>**/
#define Min(i, j)  ((i) < (j) ? (i) : (j))
#define Max(i, j)  ((i) > (j) ? (i) : (j))

namespace cghseg
{

  void compute_segmentation_mixt(numlib_vector *x, int Kmax, int lmin, int lmax, numlib_vector *phi, int P, numlib_vector **J_est, numlib_matrix **t_est) ;


class Segmentation_mixt
{
  public:
  int      _lengthx;     // length of the data
  int     _lmin;         // minimum size for a segment
  int     _lmax;         // maximum size for a segment
  double  *_x;           // data
  int      _Kmax;        // maximum number of segments
  int     _P;            // number of clusters 
  double *_phi;          // mixture parameters
  double **_Cost;        // Cost matrix
  double **_D;           // matrix of the lengths of the minimum distances
  int    **_Breaks;      // possible breakpoints
  int    **_BestBreaks;  // best breakpoints
  Segmentation_mixt(int n, int k, int minlength, int maxlength, int nbclust);
  void Init(double *Data, double *param);
  void DynProg(int k);
  void Terminate();
  friend std::ostream & operator << (std::ostream &s, const Segmentation_mixt &Seg_mixt);
  ~Segmentation_mixt();
};
}
#endif
