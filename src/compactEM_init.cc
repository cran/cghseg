#include "compactEM_init.h"
#include <fstream>
#include <list>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <memory>
#include <vector>
#include <cstring>

#include <Rconfig.h>
#ifdef SUPPORT_OPENMP
  #include <omp.h>
#endif

//#define  Dist(i,j) _D[i-2][j-1]
//#define  nk(i)     _nk[i-1]
//#define  mk(i)     _mk[i-1]
//#define  vk(i)     _vk[i-1]
//#define  Dtmp(i)   _Dtmp[i-1]


using namespace std;

namespace cghseg
{

  void compactEM_init::compute_phi(){
    vector<int> index(_K, 0);  
    for (int k=0;k<_K;k++)
      index[k] = 0;
    for (int k=0;k<_K;k++){
      double tmp    = 0;
      double mintmp = _nk0[k]*(_vk0[k]+_mk0[k]*_mk0[k])-2*_nk0[k]*_mk0[k]*_mk[0] + _nk0[k]*_mk[0]*_mk[0];
      int h         = 0;
      for (int p=1;p<_P;p++){
	tmp = _nk0[k]*(_vk0[k]+_mk0[k]*_mk0[k])-2*_nk0[k]*_mk0[k]*_mk[p] + _nk0[k]*_mk[p]*_mk[p];
	if (tmp < mintmp){
	  mintmp = tmp;
	  h = p;
	}
      }    
      _tau[k][h] = 1;
      index[k] = h;
    }
    double tmp=0;
    for (int k=0;k<_K;k++){
      int j = index[k];
      tmp += _nk0[k]*(_vk0[k]+_mk0[k]*_mk0[k])-2*_nk0[k]*_mk0[k]*_mk[j] + _nk0[k]*_mk[j]*_mk[j];
    }    
    for (int p=0;p<_P;p++){
      _phi[p] = _mk[p];
      _phi[p+_P] = sqrt(tmp / _lengthx) ;
      _phi[p+2*_P] = double(_nk[p])/_lengthx;
    }
  }

  void compactEM_init::CAH(){
    for (int k=_K;k>=_P+1;k--){
      int imin    = 2; 
      int jmin    = 1; 
      double dmin = Dist(2,1);
      CAHmin(k, imin, jmin, dmin);

      int    ntmp = nk(imin)+nk(jmin);    
      double mtmp = (nk(imin)*mk(imin)+nk(jmin)*mk(jmin))/ntmp;
      double vtmp = ( nk(imin)*vk(imin) + nk(jmin)*vk(jmin) + nk(imin)*(mk(imin)-mtmp)*(mk(imin)-mtmp)+nk(jmin)*(mk(jmin)-mtmp)*(mk(jmin)-mtmp) ) / ntmp;      
      mk(jmin) = mtmp;
      vk(jmin) = vtmp;
      nk(jmin) = ntmp;
  
      int i1 = 0;
      int i2 = 0;
      int i3 = 0;

      vector<int> index;

      if (jmin>1){ // 1:(jmin-1)
	i1     = jmin-1;     
	for (int h = 1; h<i1+1 ;h++){
	    CAHcore(h, ntmp, mtmp, vtmp);
	}
      } else {
	i1 = 1;
        CAHcore(1, ntmp, mtmp, vtmp);
      }
      if (jmin+1<=imin-1){ // (jmin+1):(imin-1)
	i2     = imin-1-jmin-1+1;
	for (int h=1;h<(i2+1);h++)
          CAHcore(h+jmin, ntmp, mtmp, vtmp);
      } 
      if (imin+1<=(k)){// (imin+1):k
	i3     = k-imin-1+1;
	for (int h=1;h<(i3+1);h++)
          CAHcore(h+imin, ntmp, mtmp, vtmp);
      }  
    
      double auxm = mk(imin); mk(imin) = mk(k);  mk(k) = auxm;
      double auxv = vk(imin); vk(imin) = vk(k);  vk(k) = auxm;
      int    auxn = nk(imin); nk(imin) = nk(k);  nk(k) = auxn;
      double auxk = Dtmp(k);
      Dtmp(k)     = Dtmp(imin);
      Dtmp(imin)  = auxk; 
      CAHcopy(k, imin, jmin);
    
    } // end k
  } //end CAH



  void compactEM_init::Init(double *Datak, double *Data2k,double *Datank){

    memcpy(_mk,Datak,sizeof(double)*_K);
    memcpy(_m2k,Data2k,sizeof(double)*_K);
    memcpy(_nk,Datank,sizeof(double)*_K);
  
    for (int k = 0 ;  k < _K; k++){
	_vk[k] = _m2k[k]-_mk[k]*_mk[k];
	_lengthx += _nk[k];
    }

    for (int k = 0 ;  k < _K; k++){
      _mk0[k] = _mk[k];
      _vk0[k] = _vk[k];    
      _nk0[k] = _nk[k];
    }
  
    for (int k = 0 ;  k < _K-1; k++){
      for (int r = 0 ;  r < k+1; r++){
	double ybar    = (_nk[k+1]*_mk[k+1]+_nk[r]*_mk[r])/(_nk[k+1]+_nk[r]);
	double varpool =  (  _nk[r]*_vk[r] + _nk[k+1]*_vk[k+1] + _nk[r]*(_mk[r]-ybar)*(_mk[r]-ybar) + _nk[k+1]*(_mk[k+1]-ybar)*(_mk[k+1]-ybar) ) / (_nk[r]+_nk[k+1]);
	_D[k][r]        = (_nk[k+1]+_nk[r])*varpool -_nk[k+1]*_vk[k+1]-_nk[r]*_vk[r];
      }
    }
  } //end Init


  compactEM_init::compactEM_init(int nbsegments, int nbclusters, int OMP_NUM_THREADS){
#ifdef SUPPORT_OPENMP
    omp_set_num_threads(OMP_NUM_THREADS);
#endif

    _K       = nbsegments; 
    _P       = nbclusters;
    _D       = new double *[_K-1];
    _Dtmp    = new double [_K];
    _phi     = new double [3*_P];
    _tau     = new double *[_K];
    _mk      = new double[_K];
    _m2k      = new double[_K];
    _vk      = new double[_K];
    _nk      = new double[_K];
    _mk0     = new double[_K];
    _vk0     = new double[_K];
    _nk0     = new int[_K];
  
    for (int p = 0; p < 3*_P; p++)
      _phi[p] = 0;
  
    _tau[0] = new double[_K*_P];
    for (int k =1; k<_K; k++)
      _tau[k] =  _tau[k-1] + _P;    
    for (int k =0; k<_K; k++){
      for (int p =0; p<_P; p++){
	_tau[k][p] = 0;
      }
    }

    int nbD = ((_K-1)*_K)/2; // 1 -> _K-1 elements
    _D[0] = new double[nbD];
    for (int k =1; k<_K-1; k++){
      _D[k] = _D[k-1]+ k;
    }
    // for (int k =0; k<_K-1; k++){
    //   for (int r=0; r<k+1; r++){
    // 	_D[k][r] = 0;
    //   }
    // }

#pragma omp parallel if (nbD>1000)
    { // page placement by first touch
#pragma omp for
      for(int ii=0; ii<nbD; ++ii)
	{
	  *(_D[0]+ii) = 0.;
	}
    }
  
    for (int k =0; k<_K; k++){
      _Dtmp[k] = 0;
      _mk[k]   = 0;
      _m2k[k]   = 0;
      _vk[k]   = 0;
      _nk[k]   = 0;
      _mk0[k]  = 0;
      _vk0[k]  = 0;
      _nk0[k]  = 0;
    }
  
  } // end constructor


  // destructor

  compactEM_init::~compactEM_init(){
  
    delete[] _phi;
    delete[] _mk;
    delete[] _m2k;
    delete[] _vk;
    delete[] _nk;
    delete[] _mk0;
    delete[] _vk0;
    delete[] _nk0;
    delete[] _Dtmp;
  
    delete[] _tau[0];
    delete[] _tau;
    delete[] _D[0];
    delete[] _D; 
  } // end destructor

  void   
  compute_compactEM_init(numlib_vector *xk,numlib_vector *x2k, numlib_vector *nk, int P , numlib_vector **phiend, int OMP_NUM_THREADS)
  {
        
    double *datak=new double[xk->size];
    for (int i=0;i<xk->size;i++)
      datak[i]=numlib_vector_get(xk,i);

    double *data2k=new double[x2k->size];
    for (int i=0;i<x2k->size;i++)
      data2k[i]=numlib_vector_get(x2k,i);

    double *datank=new double[nk->size];
    for (int i=0;i<nk->size;i++)
      datank[i]=numlib_vector_get(nk,i);

    int K   = nk->size;

    compactEM_init compactEMi(K,P,OMP_NUM_THREADS);
    compactEMi.Init(datak,data2k,datank);
    compactEMi.CAH();
    compactEMi.compute_phi();
           
    (*phiend)=numlib_vector_calloc(3*P);
    for (int p=0;p<3*P;p++)
      numlib_vector_set((*phiend),p,compactEMi._phi[p]);
    
    delete[] datak;
    delete[] data2k;
    delete[] datank;

  }  // end compute_EM_init

} // end namespace cghseg






