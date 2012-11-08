#ifndef CGHSEG_COMPACTEMINIT_H
#define CGHSEG_COMPACTEMINIT_H

#include "numlib.h"

#include <iostream>
#include <limits>
#include <fstream>
#include <cmath>

#define  Dist(i,j) _D[i-2][j-1]
#define  nk(i)     _nk[i-1]
#define  mk(i)     _mk[i-1]
#define  vk(i)     _vk[i-1]
#define  Dtmp(i)   _Dtmp[i-1]

namespace cghseg
{
  void compute_compactEM_init(numlib_vector *xk,numlib_vector *x2k, numlib_vector *nk, int P , numlib_vector **phiend, int OMP_NUM_THREADS = 1);
       
  class compactEM_init{
  private:
    inline void CAHmin(int k, int& imin, int& jmin, double& dmin);
    inline void CAHcore(int i, int ntmp, double mtmp, double vtmp);
    inline void CAHcopy(int k, int imin, int jmin);
    
  public:
    int       _lengthx;        // size of the data
    int       _K;              // Number of segments 
    int       _P;              // Number of clusters
    double   *_phi;            // parameters
    double   *_xk;              // data
    double   *_mk;             // empirical means
    double   *_m2k;             // empirical means
    double   *_vk;             // empirical var
    double   *_nk;             // length of segments
    double   *_mk0;             // empirical means
    double   *_vk0;             // empirical var
    int      *_nk0;             // length of segments
    double   **_D;             // distance matrix
    double   *_Dtmp;           // distance vector of the merged groups
    double   **_tau;           // posterior
    
    compactEM_init(int nbseg, int nbclust, int OMP_NUM_THREADS=1);
    void CAH();
    void compute_phi();
    void Init(double *Datak, double *Data2k,double *Datank);    
    friend std::ostream & operator << (std::ostream &s, const compactEM_init & compactEMi);
    ~compactEM_init();
  };

  void compactEM_init::CAHcore(int i, int ntmp, double mtmp, double vtmp){
    double ybar    =  (nk(i)*mk(i)+ntmp*mtmp)/(nk(i)+ntmp);
    double varpool =  (  ntmp*vtmp + nk(i)*vk(i) + ntmp*(mtmp-ybar)*(mtmp-ybar) + nk(i)*(mk(i)-ybar)*(mk(i)-ybar) ) / (ntmp+nk(i));
    Dtmp(i)        = -nk(i)*vk(i)-ntmp*vtmp+ (nk(i)+ntmp)*varpool;
  }

  void compactEM_init::CAHmin(int k, int& imin, int& jmin, double& dmin){
    int nbD = ((k-1)*k)/2; // 1 -> _K-1 elements 
    int iimin;    
    int iimin_shared;
    double dmin_p = std::numeric_limits<double>::max();
    double dmin_shared = std::numeric_limits<double>::max();

#pragma omp parallel if (nbD>1000) shared(dmin_shared, iimin_shared) private(iimin) firstprivate(dmin_p)
    {
#pragma omp for
      for(int ii=0; ii<nbD; ++ii)
	{
	  double distij = *(_D[0]+ii);
	  if(distij<=dmin_p){
	    dmin_p = distij;
	    iimin = ii;
	  }
	}
#pragma omp critical 
      {
	if(dmin_p<=dmin_shared){
	  dmin_shared = dmin_p;
	  iimin_shared = iimin;	
	}
      }
    } 
    dmin = dmin_shared;
    //std::cerr<<dmin<<" "<<iimin_shared<<std::endl;
    imin = int(floor(-0.5+sqrt(0.25+2*iimin_shared)));
    jmin = iimin_shared-(imin+1)*imin/2;
    imin += 2; jmin += 1;
    //std::cerr<<dmin<<" "<<imin<<" "<<jmin<<std::endl<<std::endl;    
  }

  void compactEM_init::CAHcopy(int k, int imin, int jmin){
    for (int i=1; i<(jmin-1+1);i++)
      Dist(jmin,i) = Dtmp(i);
    for (int i=jmin+1;i<(k-1+1);i++)
      Dist(i,jmin) = Dtmp(i);
    for (int j=1;j<(jmin-1+1);j++)
      Dist(imin,j) = Dist(k,j);    
    for (int j=jmin+1;j<(imin-1+1);j++)
      Dist(imin,j) = Dist(k,j);
    for (int i=imin+1;i<(k-1+1);i++)
      Dist(i,imin) = Dist(k,i);
  }

}
#endif



