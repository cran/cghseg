#include "EM_init.h"
#include <fstream>
#include <list>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <memory>
#include <vector>
#include <cstring>

#define  Dist(i,j) _D[i-2][j-1] 
#define  nk(i)     _nk[i-1]
#define  mk(i)     _mk[i-1]
#define  vk(i)     _vk[i-1]
#define  Dtmp(i)   _Dtmp[i-1]


using namespace std;

namespace cghseg
{

  void EM_init::compute_phi(){
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

  void EM_init::CAH(){
    for (int k=_K;k>=_P+1;k--){
      int imin    = 2; 
      int jmin    = 1; 
      double dmin = Dist(2,1);      
      for (int i=2; i<k+1 ; i++){
	for (int j=1; j<i ; j++){
	  if (Dist(i,j)<=dmin){
	    dmin = Dist(i,j); 
	    imin = i; 
	    jmin = j;
	  }
	}
      }

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
	  index.push_back(h);
	}
      } else {
	i1 = 1;
	index.push_back(1); 
      }
      if (jmin+1<=imin-1){ // (jmin+1):(imin-1)
	i2     = imin-1-jmin-1+1;
	for (int h=1;h<(i2+1);h++)
	  index.push_back(h+jmin);
      } 
      if (imin+1<=(k)){// (imin+1):k
	i3     = k-imin-1+1;
	for (int h=1;h<(i3+1);h++)
	  index.push_back(h+imin);
      }  
    
      int size_index= index.size();
    
      for (int h=1;h<size_index+1;h++){
	int i          = index[h-1];
	double ybar    =  (nk(i)*mk(i)+ntmp*mtmp)/(nk(i)+ntmp);
	double varpool =  (  ntmp*vtmp + nk(i)*vk(i) + ntmp*(mtmp-ybar)*(mtmp-ybar) + nk(i)*(mk(i)-ybar)*(mk(i)-ybar) ) / (ntmp+nk(i));
	Dtmp(i)        = -nk(i)*vk(i)-ntmp*vtmp+ (nk(i)+ntmp)*varpool;	
      }
      double auxm = mk(imin); mk(imin) = mk(k);  mk(k) = auxm;
      double auxv = vk(imin); vk(imin) = vk(k);  vk(k) = auxm;
      int    auxn = nk(imin); nk(imin) = nk(k);  nk(k) = auxn;
      double auxk = Dtmp(k);
      Dtmp(k)     = Dtmp(imin);
      Dtmp(imin)  = auxk;  
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
    
    } // end k
  } //end CAH


  void EM_init::Init(double *Data, int *RuptVect){
    memcpy(_x,Data,sizeof(double)*_lengthx);
  
    for (int k = 0 ;  k < _K+1; k++)
      _Breaks[k] = RuptVect[k];
  
    for (int k = 0 ;  k < _K; k++){
      int start  = _Breaks[k];
      int end    = _Breaks[k+1]-1;
      _nk[k]     = end-start+1  ;
      double *xt = &_x[start];
      for (int t = start ; t<end+1 ; t++, xt++){
	_mk[k]  += *xt;
	_vk[k]  += (*xt)*(*xt);
      }
      xt      = &_x[start];
      _mk[k] /= end-start+1 ;
      _vk[k] /= end-start+1 ;
      _vk[k] -= _mk[k]*_mk[k];
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


  EM_init::EM_init(int n, int nbsegments, int nbclusters){
    _lengthx = n ;
    _K       = nbsegments; 
    _P       = nbclusters;
    _D       = new double *[_K-1];
    _Dtmp    = new double [_K];
    _phi     = new double [3*_P];
    _tau     = new double *[_K];
    _x       = new double[_lengthx];
    _mk      = new double[_K];
    _vk      = new double[_K];
    _nk      = new int[_K];
    _mk0     = new double[_K];
    _vk0     = new double[_K];
    _nk0     = new int[_K];
    _Breaks  = new int [(_K+1)];
  
    for (int i = 0; i < _lengthx; i++)
      _x[i] = 0;
  
    for (int p = 0; p < 3*_P; p++)
      _phi[p] = 0;
  
    for (int k =0; k<_K; k++)
      _tau[k] = new double[_P];
    
    for (int k =0; k<_K; k++){
      for (int p =0; p<_P; p++){
	_tau[k][p] = 0;
      }
    }
  
    for (int k =0; k<_K-1; k++){
      _D[k] = new double[k+1];
      for (int r=0; r<k+1; r++){
	_D[k][r] = 0;
      }
    }
  
    for (int k =0; k<_K; k++){
      _Dtmp[k] = 0;
      _mk[k]   = 0;
      _vk[k]   = 0;
      _nk[k] = 0;
      _mk0[k]   = 0;
      _vk0[k]   = 0;
      _nk0[k] = 0;
    }
  
    for (int k =0; k<_K+1; k++)
      _Breaks[k] = 0;  
  } // end constructor


  // destructor

  EM_init::~EM_init(){
  
    delete[] _x;
    delete[] _phi;
    delete[] _mk;
    delete[] _vk;
    delete[] _nk;
    delete[] _mk0;
    delete[] _vk0;
    delete[] _nk0;
    delete[] _Breaks;
    delete[] _Dtmp;
  
    for (int k =0; k<_K; k++)
      delete[] _tau[k];
    delete[] _tau;
    for (int k =0; k<_K-1; k++)
      delete[] _D[k];
    delete[] _D;  
  } // end destructor


  void   
  compute_EM_init(numlib_vector *x, numlib_matrix *rupt, int P , numlib_vector **phiend)
  {
    
    
    double *data=new double[x->size];
    for (int i=0;i<x->size;i++)
      data[i]=numlib_vector_get(x,i);
    
    // WARNING : rupt = matrix (begin,end)
    // Breaks = vector: rupt[1,] with t0,...,tK with t0=1 and tk=n+1
    // -- > of size K+1
    
    int K   = rupt->size1;
    int *Bp = new int[K+1];
    for (int k=0;k<K;k++)
      Bp[k]=(int)numlib_matrix_get(rupt,k,0)-1;
    Bp[K] = (x->size);


    EM_init EMi(x->size,K,P);
    EMi.Init(data,Bp);
    EMi.CAH();
    EMi.compute_phi();
           
    (*phiend)=numlib_vector_calloc(3*P);
    for (int p=0;p<3*P;p++)
      numlib_vector_set((*phiend),p,EMi._phi[p]);
    
    delete[] data;
    delete[] Bp;
  }  // end compute_EM_init

} // end namespace cghseg







