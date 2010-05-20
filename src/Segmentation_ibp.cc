#include "Segmentation_ibp.h"
#include <list>
//#include <fstream>
#include <cmath>
#include <cstdlib>
#include <memory>


using namespace std;


namespace cghseg
{

  ostream & operator << (ostream &s, const Segmentation_ibp &Seg_ibp)
  {
  }
  
  
  // Reccurrence : calculates level m, assuming that level (m - 1) has already been done
  void Segmentation_ibp::DynProg(int m)
  {
    int Km   = _Kseq[m];
    int Kmm1 = 0;
    for (int i = 0; i < m ; i++)
      Kmm1 += _Kseq[i];
    
    if (Kmm1>_multiKmax){
      Kmm1 = _multiKmax;
    }

    _J[m][0] = _J[m-1][0] + _Contrast[m][0];
  
    for (int k = m ; k < Km+m; k++)
      { double aux   = 0.0 ;
	double lemin = _J[m][0];
	int posmin   = m-1;      
	for (int  h = 0; h< (k+1)-(m+1)+1; h++)
	  {int i = k-h-1;
	    if (i-m+2<=Kmm1-2)
	      { aux   = _J[m-1][i-m+1] + _Contrast[m][h]; 
		if (aux<lemin) {lemin = aux;posmin = i;} 
	      } 
	  } 
	_J[m][k-m]  = lemin ; 
	_Breaks[m][(k-m)] = posmin;        
      }/*end for k*/

    int i = 0;
    for (int k = Kmm1+Km-1 ;  k > Km+m-1; k--)
      { 
	i += 1;
	double aux = 0.0;
	double lemin = _J[m-1][k-Km-m+1]+_Contrast[m][Km-1];
	int posmin = Kmm1-i;
	for (int h = -1 ; h < Km-1; h++)
	  { 
	    int j = Kmm1+Km-(k+1+h+1);
	    if (j>=0){
	      aux = _J[m-1][Kmm1-j-m]+_Contrast[m][Km-h-2] ;
	      if (aux <lemin){lemin = aux; posmin = Kmm1-j-1;}
	    } 
	  }  
	_J[m][(k-m)]      = lemin;
	_Breaks[m][(k-m)] = posmin;
      } 

  }

  // fills _Contrast and the first row of matrix _J
  void Segmentation_ibp::Init(double **Cont)
  {

    for (int m = 0; m < _M; m++)
      for (int i = 0; i < _Kseq[m]; i++)
	_Contrast[m][i] = Cont[m][i];
   
    for (int i = 0; i < _Kseq[0]; i++)
      _J[0][i] = _Contrast[0][i];


  }

  Segmentation_ibp::Segmentation_ibp(int nbseries, int *seqk, int mKmax)
  {
    _M          = nbseries; 
    _multiKmax  = mKmax;
    _Kseq       = seqk;
    _Contrast   = new double *[_M];
    _J          = new double *[_M];
    _Breaks     = new int    *[_M];
    _BestBreaks = new int    *[_M];

    for (int m = 0; m < _M; m++)
      _Contrast[m] = new double[_Kseq[m]];
    for (int m = 0; m < _M; m++)
      for (int i = 0; i < _Kseq[m]; i++)
	_Contrast[m][i] = 0;  

    _J[0] = new double[_Kseq[0]];
    for (int i=0; i<_Kseq[0] ; i++)
      _J[0][i] = 0;

    for (int m = 0; m < _M; m++)
      {

	int Km   = _Kseq[m];
	int Kmm1 = 0;
	for (int i = 0; i < m ; i++)
	  Kmm1 += _Kseq[i];
        
        if (Kmm1>_multiKmax){
          Kmm1 = _multiKmax;
        }
	_J[m]      = new double[Kmm1+Km-(m+1)+1];    
	_Breaks[m] = new int[Kmm1+Km-(m+1)+1];

	for (int i=0;i< Kmm1+Km-(m+1)+1; i++)
	  {
	    _J[m][i]      = 0;
	    _Breaks[m][i] = 0;
	  }

      }
  
    
    for (int m = 0; m< _M ; m++)
      _BestBreaks[m] = new int[_multiKmax-_M+1];
    for (int m = 0 ; m < _M; m++)
      for (int i = 0; i < _multiKmax-_M+1; i++)
	_BestBreaks[m][i] = -5 ;
   
  }

  // finds the Beast Breakpoints 
  void Segmentation_ibp::Terminate(int Ktot)
  {

    for (int k = Ktot-1; k > _M - 2; k--)
      {   
	int tmp1 = _Breaks[_M-1][k-_M +1];
	_BestBreaks[_M-1][k -_M +1] = k - tmp1 - 1;
	for (int m = _M-2; m > 0; m--)
	  {
	    int tmp2 = _Breaks[m][tmp1-m];
	    _BestBreaks[m][k-_M+1] = tmp1 - tmp2 - 1;
	    tmp1 = tmp2;
	  }
      }
    for (int k=0 ; k < Ktot-_M+1; k++)
      { 
	int Ktmp = 0;
	for (int m=1; m<_M; m++)
	  {
	    Ktmp += _BestBreaks[m][k];
	  }
	_BestBreaks[0][k] = k - Ktmp; 
      }
   
    
  }

  // destructor

  Segmentation_ibp::~Segmentation_ibp()
  {

    for (int m = 0; m < _M; m++)
      delete[] _Contrast[m] ;
    delete[] _Contrast;

    for (int m = 0; m < _M; m++)
      delete[] _J[m] ;
    delete[] _J;

    for (int m = 0; m < _M; m++)
      delete[] _Breaks[m] ;
    delete[] _Breaks;

    for (int m = 0; m < _M; m++)
      delete[] _BestBreaks[m] ;
    delete[] _BestBreaks;
    

  }


  void compute_segmentation_ibp(numlib_vector *Jvector, numlib_vector_int *Kvector, int multiK, numlib_vector **J_est, numlib_matrix **t_est)
  {
    double *ContrastVector=new double[Jvector->size];
    for (int i=0;i<Jvector->size;i++)
      ContrastVector[i]=numlib_vector_get(Jvector,i);

    int *Kseq =new int[Kvector->size];
    for (int i=0;i<Kvector->size;i++)
      Kseq[i]=numlib_vector_int_get(Kvector,i);


    int M               = Kvector->size;
    double ** Contrast  = new double*[M];
    int a               = 0;
    int b               = 0;

    for (int m = 0; m< M ; m++)
      {
	Contrast[m] = new double[Kseq[m]];
	b = a + Kseq[m];
	for (int i = a ; i<b; i++)
	  Contrast[m][i-a] = ContrastVector[i];
	a = b;
      }
    
    int multiKmax = multiK;

    Segmentation_ibp Seg_ibp(M,Kseq,multiKmax);
    Seg_ibp.Init(Contrast);
    for (int m = 1; m < M; m++)
      Seg_ibp.DynProg(m);
    Seg_ibp.Terminate(multiKmax);

    if (J_est != NULL) {
      (*J_est)=numlib_vector_alloc(multiKmax-M+1);
      for (int k=0; k<multiKmax-M+1; k++){	        
	numlib_vector_set((*J_est),k,Seg_ibp._J[M][k]);
      }
    }
    
    (*t_est)=numlib_matrix_calloc(M,multiKmax-M+1);
    for (int m=0;m<M;m++)
       for (int k=0;k<multiKmax-M+1;k++)
	numlib_matrix_set((*t_est),k,m,Seg_ibp._BestBreaks[m][k]+1);
    
    delete[] Kseq;
    delete[] ContrastVector;
    for (int m = 0; m < M; m++)
      delete[] Contrast[m] ;
    delete[] Contrast;

  }
     





} // end namespace seglm //
