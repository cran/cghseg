#include "Segmentation_mixt.h"
/**#include <fstream>**/
#include <list>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <math.h>
#include <cstring>
#include <limits>
#include <vector>


using namespace std;

namespace cghseg
{

  void Segmentation_mixt::DynProg(int k)
{
  int t = min((k+1)*_lmax, _lengthx);

  for (int i = (k+1)*_lmin - 1; i < t ; i++)
  {
    int threshold = _lmax-_lmin+2;

    if (i < threshold){
      double  *Ci  = _Cost[i-_lmin +1]; 
      int pos_min  = k;
      int       a  = k*_lmin;
      int       b  = min(k*_lmax+1,i-_lmin+2);
      double   min = numeric_limits<double>::max();

      for (int h = a; h < b; h++){
	if (_D[k - 1][(h - 1)-k*_lmin+1] + Ci[h] < min)
	  {
	    min     = _D[k - 1][(h - 1) - k*_lmin+1] + Ci[h];
	    pos_min = h;
	  }
      }    
      _D[k][i-(k+1)*_lmin+1] = min;   
      _Breaks[k-1][i-(k+1)*_lmin+1] = pos_min - 1;      

    } else {
      double  *Ci  = _Cost[i-_lmin +1]; 
      int pos_min  = k;
      int       a  = k*_lmin;
      int       b  = min(k*_lmax+1,i-_lmin+2);
      int       h2 = 0;
      double   min = numeric_limits<double>::max();

      for (int h = a; h < b; h++){
	int h2 = h - i + threshold - 1;       	  
	if ((h2>=0)){
	  if (_D[k - 1][(h - 1)-k*_lmin+1] + Ci[h2] < min)
	    {
	      min     = _D[k - 1][(h - 1) - k*_lmin+1] + Ci[h2];
	      pos_min = h;
	    }
	} 
      } /*for h*/ 
      _D[k][i-(k+1)*_lmin+1] = min;   
      _Breaks[k-1][i-(k+1)*_lmin+1] = pos_min - 1;  
    } /*end else*/

  }

}

// fills vector _x and the first column of matrix _D
void Segmentation_mixt::Init(double *Data, double *param)
{

  double aux,Cij;
  double pi = acos(-1);
  memcpy(_x,Data,sizeof(double)*_lengthx);
  memcpy(_phi,param,sizeof(double)*3*_P);

  for (int i = 0; i < _lmax - _lmin+1; i++)
    for (int j = 0; j < i + 1; j++)
    {
      double ykbar = 0 ;
      double y2kbar = 0; 
      double wk =0;
      Cij=_Cost[i][j];
      double ni   = i - j + _lmin ;	
      double *xt  = &_x[j];
      
      for (int t = j; t < i + _lmin ; t++, xt++)
      {
        ykbar  += *xt;
        y2kbar += (*xt)*(*xt);
      }
      xt      = &_x[j];
      ykbar  /= ni;
      y2kbar /= ni;
      wk      = y2kbar - ykbar*ykbar ;

     
      vector<double> dkp(_P, 0);  



      for (int p=0; p<_P; p++) 
        dkp[p] = (ykbar - _phi[p])*(ykbar - _phi[p]);

      vector<double> fyk(_P, 0);  

      for (int p=0; p<_P; p++) 
      {
        fyk[p] = 0.5* ni * (-(dkp[p]+wk)/(_phi[_P+p]*_phi[_P+p]) -log(2*pi*_phi[_P+p]*_phi[_P+p]) ) + log(_phi[2*_P+p]);
      }
      double max = fyk[0];
      for (int p=1; p<_P; p++)
      {
      if (fyk[p]>max) 
        max = fyk[p];
      }

      for (int p=0; p<_P; p++)
        fyk[p] = exp(fyk[p]-max);

      for (int p=0; p<_P; p++)
        Cij += fyk[p];

      Cij = -log(Cij)-max;
     _Cost[i][j]=Cij; 
    }
  if (_lmax<_lengthx)
    {
      for (int i = _lmax-_lmin+1; i < _lengthx -_lmin + 1; i++)
	for (int j = 0 ; j < _lmax -_lmin + 1; j++)
	  {
	    double ykbar = 0 ;
	    double y2kbar = 0; 
	    double wk =0;
	    Cij=_Cost[i][j];
	    
	    double ni   = _lmax - j;
	    int       h = j+i-_lmax+_lmin;
	    double *xt  = &_x[h];
	    
	    for (int t = h; t < i + _lmin ; t++, xt++)
	      {
		ykbar  += *xt;
		y2kbar += (*xt)*(*xt);
	      }
	    xt      = &_x[h];
	    ykbar  /= ni;
	    y2kbar /= ni;
	    wk      = y2kbar - ykbar*ykbar ;
	    
	    vector<double> dkp(_P, 0);  
	    
	    for (int p=0; p<_P; p++) 
	      dkp[p] = (ykbar - _phi[p])*(ykbar - _phi[p]);
	    
	    vector<double> fyk(_P, 0);  
	    for (int p=0; p<_P; p++) 
	      {
		fyk[p] = 0.5* ni * (-(dkp[p]+wk)/(_phi[_P+p]*_phi[_P+p]) -log(2*pi*_phi[_P+p]*_phi[_P+p]) ) + log(_phi[2*_P+p]);
	      }
	    double max = fyk[0];
	    for (int p=1; p<_P; p++)
	      {
		if (fyk[p]>max) 
		  max = fyk[p];
	      }
	    
	    for (int p=0; p<_P; p++)
	      fyk[p] = exp(fyk[p]-max);
	    
	    for (int p=0; p<_P; p++)
	      Cij += fyk[p];
	    
	    Cij = -log(Cij)-max;
	    _Cost[i][j]=Cij; 
	  }
    }

  for (int i = 0; i < _lmax-_lmin+1; i++)
    _D[0][i] = _Cost[i][0];
  
}








Segmentation_mixt::Segmentation_mixt(int n, int c, int minlength, int maxlength,  int nbclust)
{
  _lengthx    = n; 
  _Kmax       = c;
  _P          = nbclust;
  _lmin       = minlength;
  _lmax       = maxlength;
  _x          = new double[_lengthx];
  _phi        = new double[3*_P];
  _Cost       = new double *[_lengthx-_lmin+1];
  _D          = new double *[_Kmax];
  _Breaks     = new int  *[_Kmax-1];
  _BestBreaks = new int  *[_Kmax+1];

  for (int i = 0; i < _lengthx; i++)
    _x[i] = 0;

  for (int p = 0; p < 3*_P; p++)
    _phi[p] = 0;

  for (int i = 0; i < _lengthx-_lmin+1; i++)
    _Cost[i] = new double[i + 1];

  for (int i = 0; i < _lmax-_lmin+1; i++)
    for (int j = 0; j < i + 1; j++)
      _Cost[i][j] = 0;  

  if (_lmax<_lengthx){
    for (int i=_lmax-_lmin+1; i < _lengthx-_lmin+1 ; i++)
      for (int j = 0 ; j < _lmax- _lmin +1; j++)
	_Cost[i][j] = 0 ;
  } /*end if*/
  

  for (int k = 0; k< _Kmax;k++){
    int t = min((k+1)*_lmax,_lengthx);
/*    _D[k] = new double[_lengthx - (k+1)* _lmin +1];*/
      _D[k] = new double[t - (k+1)* _lmin +1];
    for (int i = 0; i < t - (k+1)* _lmin +1; i++)
      _D[k][i] = 0;  
  }

  for (int k = 0; k < _Kmax-1; k++){
    int t = min((k+2)*_lmax,_lengthx);
/*    _Breaks[k] = new int [_lengthx -(k+2)* _lmin +1];*/
    _Breaks[k] = new int [ t  -(k+2)* _lmin +1];
    for (int i = 0; i < t  - (k+2)* _lmin +1 ; i++)
      _Breaks[k][i] = 0;
  }

  for (int k = 0; k < _Kmax; k++)
    _BestBreaks[k] = new int [k+1];
  for (int k = 0; k < _Kmax; k++)
    for (int h = 0; h < k+1; h++)
      _BestBreaks[k][h] = 0;
}


// finds the Beast Beakpoints 
void Segmentation_mixt::Terminate()
{
  for (int k = 0; k < _Kmax ; k++)
    _BestBreaks[k][k] = min((k+1)*_lmax - 1, _lengthx - 1 );
  for (int k = 1; k < _Kmax; k++)
    for (int i = k - 1; i >= 0; i--)
    {
      _BestBreaks[k][i] = _Breaks[i][ _BestBreaks[k][i + 1] -(i+2)*_lmin+1];
    }

}


// destructor

Segmentation_mixt::~Segmentation_mixt()
{

 delete[] _x;
 for (int i = 0; i < _lengthx -_lmin+1; i++)
    delete[] _Cost[i];
 delete[] _Cost;
 for (int k = 0; k< _Kmax; k++)
    delete[] _D[k];
 delete[] _D;
 for (int k = 0; k<_Kmax-1; k++)
    delete[] _Breaks[k];
 delete[] _Breaks;

}


void
compute_segmentation_mixt(numlib_vector *x, int Kmax, int lmin, int lmax, numlib_vector *phi, int P, numlib_vector **J_est, numlib_matrix **t_est) {

	double *data=new double[x->size];
	for (int i=0;i<x->size;i++)
		data[i]=numlib_vector_get(x,i);

        int length_phi = phi->size;
	double *param  = new double[length_phi];
	for (int i=0;i<length_phi;i++)
		param[i] = numlib_vector_get(phi,i);

	Segmentation_mixt Seg_mixt(x->size,Kmax,lmin,lmax,P);
	Seg_mixt.Init(data,param);

	for (int k=1; k<Kmax; k++)
		Seg_mixt.DynProg(k);

	Seg_mixt.Terminate();

	if (J_est != NULL) {
		(*J_est)=numlib_vector_alloc(Kmax);
		for (int k=0; k<Kmax; k++)
		  {
		    int t = Min((x->size), (k+1)*lmax);
		    numlib_vector_set((*J_est),k,Seg_mixt._D[k][t-(k+1)*lmin]);
		  }
	}

	(*t_est)=numlib_matrix_calloc(Kmax,Kmax);
	for (int c=0;c<Kmax;c++)
		for (int r=0;r<Kmax;r++)
			if (r>=c) {
				numlib_matrix_set((*t_est),r,c,Seg_mixt._BestBreaks[r][c]+1);
			}

        delete[] param;
	delete[] data;

}


} /** end namespace segclust**/
