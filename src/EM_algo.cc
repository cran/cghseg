#include "EM_algo.h"
//#include <fstream>
#include <list>
#include <cmath>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <memory>
#include <limits>
#include <cstring>
using namespace std;

namespace cghseg
{
  
  void EM_algo::Estep(){
    
    double pi = acos(-1);
    vector<double> max(_K, 0);  

    
    for (int k = 0 ;  k < _K; k++){
      int nk     = _Breaks[k+1]-_Breaks[k];
      double dkp = 0;
      vector<double> logfxk(_P, 0);  
      for (int p=0; p<_P; p++) {
	dkp = (_xkbar[k] - _phi[p])*(_xkbar[k] - _phi[p]);
	logfxk[p] = 0.5* nk * (-(dkp+_wk[k])/(_phi[_P+p]*_phi[_P+p]) -log(2*pi*_phi[_P+p]*_phi[_P+p]) )+ log(_phi[2*_P+p]);                          
      }
      double tmpmax = logfxk[0];
      for (int p=1; p<_P; p++){
	if (logfxk[p]>tmpmax) 
	  tmpmax = logfxk[p];
      }     
      max[k] = tmpmax;
      for (int p=0; p<_P; p++){
	_tau[k][p] = exp(logfxk[p]-max[k]);
      }
    }
    
    double tmp = 0 ;
    _lvinc = 0;
    for (int k = 0 ;  k < _K; k++){
      tmp = 0;
      for (int p=0; p<_P; p++){
	tmp +=_tau[k][p]; 
      }
      _lvinc += log(tmp) +max[k];
    }  
    
    for (int k = 0 ;  k < _K; k++){
      tmp = 0;
      for (int p=0; p<_P; p++){
	tmp +=_tau[k][p]; 
      }
      for (int p=0; p<_P; p++){
	_tau[k][p]/=tmp; 
      }
    }  
  } //end Estep
  
  void EM_algo::Mstep(){
    
    for (int p =0;p<3*_P;p++){
      _phi[p] = 0;
    }
    for (int p =0;p<_P;p++){
      double denom=0;     
      for (int k=0;k<_K;k++){       
	int nk        = _Breaks[k+1]-_Breaks[k];
	_phi[p]      += nk*_tau[k][p]*_xkbar[k];
	_phi[2*_P+p] += _tau[k][p] ;
	denom        += nk*_tau[k][p];       
      }
      _phi[p] /=denom;
      _phi[2*_P+p] /=_K;
    }
    
    if (_vh==false){
      for (int p =0;p<_P;p++){
	double denom=0;
	double dkp;
	double s2p=0;
	for (int k=0;k<_K;k++){
	  int nk        = _Breaks[k+1]-_Breaks[k];    
	  dkp           = (_xkbar[k] - _phi[p])*(_xkbar[k] - _phi[p]);
	  s2p          += nk*_tau[k][p]*(dkp+_wk[k]);
	  denom        += nk*_tau[k][p];       
	}
	s2p /= denom;
	_phi[_P+p] = sqrt(double(s2p));
      }

      
    } else if (_vh==true){
      double s2=0;
      for (int p =0;p<_P;p++){
	double dkp;
	for (int k=0;k<_K;k++){
	  dkp    = (_xkbar[k] - _phi[p])*(_xkbar[k] - _phi[p]);
	  int nk = _Breaks[k+1]-_Breaks[k];      
	  s2    += nk*_tau[k][p]*(dkp+_wk[k]);
	}
      }
      s2 /=_lengthx;
      for (int p =0;p<_P;p++){
	_phi[_P+p] = sqrt(double(s2));  
      }
    }  
  } // end Mstep
  
  
  void EM_algo::EM(){

    vector<double> delta(3*_P, 0);  
    double maxdelta = 1e-4;
    
    int iter        = 0;
    vector<double> np(_P, 0);  
    double eps      = 10e-10;
    vector<double> phi_tmp(3*_P, 0);  
    double min_np = _lengthx;
    
    for (int p=0; p<_P; p++){
      np[p] = 0;
    }   
    
    while ( (maxdelta>=1e-4) & (min_np>eps) & (iter<=5000) ){     
      iter       += 1;     
      for (int p=0;p<3*_P;p++){
	phi_tmp[p] = _phi[p];
      }     
      Estep();          
      Mstep();          
      for (int p=0;p<_P;p++){
	for (int k=0;k<_K;k++){
	  np[p] += _tau[k][p];
	}
      }  
      for (int p=0; p<_P; p++){
	if (np[p]<min_np) 
	  min_np = np[p];
      }     
      for (int p=0; p<_P; p++){
	np[p] = 0;
      }
      for (int p=0;p<3*_P;p++){
	delta[p]= fabs((phi_tmp[p]-_phi[p])/_phi[p]);
      }
      
      for (int p=0;p<3*_P;p++){
	if (delta[p]<maxdelta) 
	  maxdelta = delta[p];
      }     
      
    } // end while
    
    
    if (min_np<eps){
      _empty = 1;
      _lvinc = - numeric_limits<double>::infinity( ) ;
    }
    
    if (iter>5000){
      _dv    = 2;
      _lvinc = - numeric_limits<double>::infinity( ) ;
    }
    
  } // end EM()
  
  //fills vector _x and then _xkbar, _x2kbar and _phi
  void EM_algo::Init(double *Data, double *param, int *RuptVect)
  {
    memcpy(_x,Data,sizeof(double)*_lengthx);
    memcpy(_phi,param,sizeof(double)*3*_P);
    
    for (int k = 0 ;  k < _K+1; k++){
      _Breaks[k] = RuptVect[k];
    }
    
    for (int k = 0 ;  k < _K; k++){
      int start = _Breaks[k];
      int end   = _Breaks[k+1]-1;
      double nk = end-start+1  ;
      double *xt = &_x[start];
      for (int t = start ; t<end+1 ; t++, xt++){
	_xkbar[k]  += *xt;
	_x2kbar[k] += (*xt)*(*xt);
      }
      xt      = &_x[start];
      _xkbar[k]  /= nk;
      _x2kbar[k] /= nk;
      _wk[k]      = _x2kbar[k] - _xkbar[k]*_xkbar[k];
    }
    
  } // end Init
  
  
  
  EM_algo::EM_algo(int n, int nbsegments, int nbclusters, bool varh)
  {
    _lengthx    = n ;
    _K          = nbsegments; 
    _P          = nbclusters;
    _vh         = varh;
    _phi        = new double [3*_P];
    _tau        = new double *[_K];
    _x          = new double[_lengthx];
    _xkbar      = new double[_K];
    _x2kbar     = new double[_K];
    _wk         = new double[_K];
    _Breaks     = new int [(_K+1)];

    
    _lvinc    = 0;
    _empty = 0;
    _dv    = 0;
    
    for (int i = 0; i < _lengthx; i++)
      _x[i] = 0;
    for (int p = 0; p < 3*_P; p++)
      _phi[p] = 0;
    for (int k =0; k<_K; k++){
      _tau[k] = new double[_P];
    }
    for (int k =0; k<_K; k++){
      for (int p =0; p<_P; p++){
	_tau[k][p] = 0;
      }
    }
    for (int k =0; k<_K; k++){
      _xkbar[k]  = 0;
      _x2kbar[k] = 0;
      _wk[k]     = 0;
    }
    
    for (int k =0; k<_K+1; k++)
      _Breaks[k] = 0;
    
  } // end constructor
  
  
  // destructor
  
  EM_algo::~EM_algo()
  {
    
    delete[] _x;
    delete[] _phi;
    delete[] _xkbar;
    delete[] _x2kbar;
    delete[] _wk;
    delete[] _Breaks;
    
    for (int k =0; k<_K; k++)
      delete[] _tau[k];
    delete[] _tau;
    
  } // end destructor
  
  
  void   
  compute_EM_algo(numlib_vector *x, numlib_matrix *rupt, int P, numlib_vector *phistart, bool vh, numlib_vector **phiend, numlib_matrix **tau, double *lvinc_pointer, int *empty_pointer, int *dv_pointer)
  {
    
    
    double *data=new double[x->size];
    for (int i=0;i<x->size;i++)
      data[i]=numlib_vector_get(x,i);
    
    double *phi0=new double[phistart->size];
    for (int p=0;p<3*P;p++)
      phi0[p]=numlib_vector_get(phistart,p);
    

    // WARNING : rupt = matrix (begin,end)
    // Breaks = vector: rupt[1,] with t0,...,tK with t0=1 and tk=n+1
    // -- > of size K+1
    
    int K   = rupt->size1;
    int *Bp = new int[K+1];
    for (int k=0;k<K;k++)
      Bp[k]=(int)numlib_matrix_get(rupt,k,0)-1;
    Bp[K] = (x->size);
    
    EM_algo EMa(x->size,K,P,vh);  
    EMa.Init(data,phi0,Bp);
    EMa.EM();
    
    
    (*tau)=numlib_matrix_calloc(K,P);
    for (int k=0;k<K;k++)
      for (int p=0;p<P;p++)
	numlib_matrix_set((*tau),k,p,EMa._tau[k][p]);
    
    (*phiend)=numlib_vector_calloc(3*P);
    for (int p=0;p<3*P;p++)
      numlib_vector_set((*phiend),p,EMa._phi[p]);
    
    *lvinc_pointer = EMa._lvinc;
    *dv_pointer    = EMa._dv;
    *empty_pointer = EMa._empty;
    
    delete[] data;
    delete[] phi0;
    delete[] Bp;
  }  // end compute_EM_algo
  
} // end namespace cghseg
  

