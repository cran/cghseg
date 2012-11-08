#include "compactEM_algo.h"
//#include <fstream>
#include <list>
#include <cmath>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <memory>
#include <limits>
#include <cstring>
#include <iostream>
using namespace std;

namespace cghseg
{
  
  void compactEM_algo::Estep(){
    
    double pi = acos(-1);
    vector<double> max(_K, 0);  

    
    for (int k = 0 ;  k < _K; k++){
      double dkp = 0;
      vector<double> logfxk(_P, 0);  
      for (int p=0; p<_P; p++) {
	dkp = (_xkbar[k] - _phi[p])*(_xkbar[k] - _phi[p]);
	logfxk[p] = 0.5* _nk[k] * (-(dkp+_wk[k])/(_phi[_P+p]*_phi[_P+p]) -log(2*pi*_phi[_P+p]*_phi[_P+p]) )+ log(_phi[2*_P+p]);                          
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
  
  void compactEM_algo::Mstep(){
    
    for (int p =0;p<3*_P;p++){
      _phi[p] = 0;
    }
    for (int p =0;p<_P;p++){
      double denom=0;     
      for (int k=0;k<_K;k++){       

	_phi[p]      += _nk[k]*_tau[k][p]*_xkbar[k];
	_phi[2*_P+p] += _tau[k][p] ;
	denom        += _nk[k]*_tau[k][p];       
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
	  dkp           = (_xkbar[k] - _phi[p])*(_xkbar[k] - _phi[p]);
	  s2p          += _nk[k]*_tau[k][p]*(dkp+_wk[k]);
	  denom        += _nk[k]*_tau[k][p];       
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
	  s2    += _nk[k]*_tau[k][p]*(dkp+_wk[k]);
	}
      }
      s2 /=_lengthx;
      for (int p =0;p<_P;p++){
	_phi[_P+p] = sqrt(double(s2));  
      }
    }  
  } // end Mstep
  
  
  void compactEM_algo::compactEM(){

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
  
  //fills vectors xkbar, _x2kbar and _phi
  void compactEM_algo::Init(double *Datak, double *Data2k, double *param, double *Datank)
  {
    memcpy(_xkbar,Datak,sizeof(double)*_K);
    memcpy(_x2kbar,Data2k,sizeof(double)*_K);
    memcpy(_nk,Datank,sizeof(double)*_K);
    memcpy(_phi,param,sizeof(double)*3*_P);

    for (int k = 0 ;  k < _K; k++){
      _wk[k]   = _x2kbar[k] - _xkbar[k]*_xkbar[k];
      _lengthx +=_nk[k];
    }   
  } // end Init
  
  
  
  compactEM_algo::compactEM_algo(int nbsegments, int nbclusters, bool varh)
  {
    
    _K          = nbsegments; 
    _P          = nbclusters;
    _vh         = varh;
    _phi        = new double [3*_P];
    _tau        = new double *[_K];
    _xkbar      = new double[_K];
    _x2kbar     = new double[_K];
    _nk         = new double[_K];
    _wk         = new double[_K];

    _lengthx = 0;
    _lvinc   = 0;
    _empty   = 0;
    _dv      = 0;
    
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
      _nk[k]     = 0;
    }
    
        
  } // end constructor
  
  
  // destructor
  
  compactEM_algo::~compactEM_algo()
  {
    
    delete[] _phi;
    delete[] _xkbar;
    delete[] _x2kbar;
    delete[] _nk;
    delete[] _wk;
    
    for (int k =0; k<_K; k++)
      delete[] _tau[k];
    delete[] _tau;
    
  } // end destructor
  
  
  void   
  compute_compactEM_algo(numlib_vector *xkbar,numlib_vector *x2kbar, numlib_vector *nk, int P, numlib_vector *phistart, bool vh, numlib_vector **phiend, numlib_matrix **tau, double *lvinc_pointer, int *empty_pointer, int *dv_pointer)
  {    

    double *datak=new double[xkbar->size];
    for (int i=0;i<xkbar->size;i++)
      datak[i]=numlib_vector_get(xkbar,i);
    
    double *data2k=new double[x2kbar->size];
    for (int i=0;i<x2kbar->size;i++)
      data2k[i]=numlib_vector_get(x2kbar,i);

    double *datank=new double[nk->size];
    for (int i=0;i<nk->size;i++)
      datank[i]=numlib_vector_get(nk,i);

    double *phi0=new double[phistart->size];
    for (int p=0;p<3*P;p++)
      phi0[p]=numlib_vector_get(phistart,p);
        
    int K   = nk->size;

    compactEM_algo compactEMa(K,P,vh);  
    compactEMa.Init(datak,data2k,phi0,datank);
    compactEMa.compactEM();
    
    
    (*tau)=numlib_matrix_calloc(K,P);
    for (int k=0;k<K;k++)
      for (int p=0;p<P;p++)
	numlib_matrix_set((*tau),k,p,compactEMa._tau[k][p]);
    
    (*phiend)=numlib_vector_calloc(3*P);
    for (int p=0;p<3*P;p++)
      numlib_vector_set((*phiend),p,compactEMa._phi[p]);
    
    *lvinc_pointer = compactEMa._lvinc;
    *dv_pointer    = compactEMa._dv;
    *empty_pointer = compactEMa._empty;
    
    delete[] datak;
    delete[] data2k;
    delete[] datank;
    delete[] phi0;

  }  // end compute_compactEM_algo
  
} // end namespace cghseg
  

