#include "Segmentation_mean.h"
//#include <fstream>
#include <list>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <cstring>
#include <limits>

using namespace std;

namespace cghseg
{
ostream & operator << (ostream &s, const Segmentation_mean &Seg_mean)
{
/**
  s << "nMax = " << Seg_mean._lengthx << endl;
  s << "Kmax = " << Seg_mean._Kmax << endl;
  s << "Affichage des donnees : " << endl;
  for (int i = 0; i < Seg_mean._lengthx; i++)
    s << Seg_mean._x[i] << " ";
  s << endl;

  s << "Affichage de la matrice _Cost : " << endl;
  for (int i = 0; i < Seg_mean._lengthx-Seg_mean._lmin+1; i++)
  {
    for (int j = 0; j < i + 1; j++)
      s << Seg_mean._Cost[i][j] << "\t";
    s << endl;
  }

  s << "Affichage de la matrice _D[k][i] : " << endl;
  for (int k = 0; k < Seg_mean._Kmax; k++)
  {
    for (int i = 0; i < Seg_mean._lengthx - (k+1)* Seg_mean._lmin +1; i++)
      s << Seg_mean._D[k][i] << "\t";
    s << endl;
  }

  s << "Affichage du contraste : " << endl;
  {
    for (int j = 0; j < Seg_mean._Kmax; j++)
      s << Seg_mean._D[j][Seg_mean._lengthx - (j+1)* Seg_mean._lmin ] << " ";
    s << endl;
  }

  s << "Affichage de la matrice des ruptures : " << endl;
  for (int k = 0; k < Seg_mean._Kmax-1; k++)
  {
    for (int i = 0; i < Seg_mean._lengthx - (k+2)* Seg_mean._lmin+1; i++)
      s << Seg_mean._Breaks[k][i] << " ";
    s << endl;
  }

  s << "Affichage de la matrice des meilleures ruptures : " << endl;
  for (int k = 0; k < Seg_mean._Kmax; k++)
  {
    for (int l = 0; l < k + 1; l++)
      s << Seg_mean._BestBreaks[k][l] << " ";
    s << endl;
  }

  return s;
**/
}

// Reccurrence : calculates level k, assuming that level (k - 1) has already been done
  void Segmentation_mean::DynProg(int k)
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
void Segmentation_mean::Init(double *Data, bool vh)
{
  double mean   = 0;
  double ykbar  = 0;
  double y2kbar = 0;
  double aux;
  double Cij;
  memcpy(_x,Data,sizeof(double)*_lengthx);

  if (vh==true){      
    for (int i = 0; i < _lmax - _lmin+1; i++)
      for (int j = 0; j < i + 1; j++)
	{
	  mean        = 0;
	  Cij         =_Cost[i][j];
	  double ni   = i - j + _lmin ;	
	  double *xt  = &_x[j];
	  for (int t = j; t < i + _lmin ; t++, xt++)
	    mean += *xt;
	  mean /= ni;
	  xt=&_x[j];
	  for (int t = j; t < i + _lmin; t++, xt++)
	    {
	      aux = *xt - mean;
	      Cij += aux * aux;
	    }
	  _Cost[i][j]=Cij;
	}
    if (_lmax<_lengthx)
      {
	for (int i = _lmax-_lmin+1; i < _lengthx -_lmin + 1; i++)
	  for (int j = 0 ; j < _lmax -_lmin + 1; j++)
	    {
	      mean        = 0;
	      Cij         =_Cost[i][j];
	      double ni   = _lmax - j;	
	      int    h    = j+i-_lmax+_lmin;
	      double *xt  = &_x[h];
	      for (int t = h ; t < i+_lmin  ; t++, xt++)
		mean += *xt;
	      mean /= ni;
	      xt=&_x[h];
	      for (int t = h; t < i+_lmin ; t++, xt++)
		{
		  aux = *xt - mean;
		  Cij += aux * aux;
		}
	      _Cost[i][j]=Cij;
	    }
      } /*end if _lmax*/
  } /*end if vh==true*/



  else /*--> if (vh==FALSE) case with heterogeneous variances */
    {
      for (int i = 0; i < _lmax - _lmin+1; i++)
	for (int j = 0; j < i + 1; j++)
	  {
	    ykbar       = 0;
	    y2kbar      = 0;
	    Cij         =_Cost[i][j];
	    double ni   = i - j + _lmin ;	
	    double *xt  = &_x[j];
	    for (int t = j; t < i + _lmin ; t++, xt++)
	      {
		ykbar  += *xt;
		y2kbar += (*xt) * (*xt);
	      }
	    ykbar      /= ni;
	    y2kbar     /= ni;
	    Cij         = ni * log (y2kbar-ykbar*ykbar);
	    xt          = &_x[j];
	    _Cost[i][j] = Cij;
	  }
      if (_lmax<_lengthx)
	{
	  for (int i = _lmax-_lmin+1; i < _lengthx -_lmin + 1; i++)
	    for (int j = 0 ; j < _lmax -_lmin + 1; j++)
	      {
		ykbar       = 0;
		y2kbar      = 0;
		Cij         =_Cost[i][j];
		double ni   = _lmax - j;
		int       h = j+i-_lmax+_lmin;
		double *xt  = &_x[h];
		for (int t = h; t < i + _lmin ; t++, xt++)
		  {
		    ykbar  += *xt;
		    y2kbar += (*xt) * (*xt);
		  }
		ykbar      /= ni;
		y2kbar     /= ni;
		Cij         = ni * log (y2kbar-ykbar*ykbar);
		xt          = &_x[h];
		_Cost[i][j]=Cij;
	      }
	}  
    } 

  
  for (int i = 0; i < _lmax-_lmin+1; i++)
    _D[0][i] = _Cost[i][0];

}

 Segmentation_mean::Segmentation_mean(int n, int c, int minlength, int maxlength, bool varh)
{
  _lengthx    = n; 
  _Kmax       = c;
  _lmin       = minlength;
  _lmax       = maxlength;
  _vh         = varh;
  _x          = new double[_lengthx];
  _Cost       = new double *[_lengthx-_lmin+1];
  _D          = new double *[_Kmax];
  _Breaks     = new int  *[_Kmax-1];
  _BestBreaks = new int  *[_Kmax+1];

  for (int i = 0; i < _lengthx; i++)
    _x[i] = 0;


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
      _D[k] = new double[t - (k+1)* _lmin +1];
    for (int i = 0; i < t - (k+1)* _lmin +1; i++)
      _D[k][i] = 0;  
  }

  for (int k = 0; k < _Kmax-1; k++){
    int t = min((k+2)*_lmax,_lengthx);
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
void Segmentation_mean::Terminate()
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

Segmentation_mean::~Segmentation_mean()
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
compute_segmentation_mean(numlib_vector *x, int Kmax, int lmin, int lmax, bool vh, numlib_vector **J_est, numlib_matrix **t_est) {

	double *data=new double[x->size];
	for (int i=0;i<x->size;i++)
		data[i]=numlib_vector_get(x,i);

	Segmentation_mean Seg_mean(x->size,Kmax,lmin,lmax,vh);
	Seg_mean.Init(data,vh);

	for (int k=1; k<Kmax; k++)
		Seg_mean.DynProg(k);

	Seg_mean.Terminate();

	if (J_est != NULL) {
		(*J_est)=numlib_vector_alloc(Kmax);
		for (int k=0; k<Kmax; k++){	        
			int t = Min((x->size), (k+1)*lmax);
			numlib_vector_set((*J_est),k,Seg_mean._D[k][t-(k+1)*lmin]);
	   	}
	}

	(*t_est)=numlib_matrix_calloc(Kmax,Kmax);
	for (int c=0;c<Kmax;c++)
		for (int r=0;r<Kmax;r++)
			if (r>=c) {
				numlib_matrix_set((*t_est),r,c,Seg_mean._BestBreaks[r][c]+1);
			}

	delete[] data;


}

} // end namespace cghseg
