
#include <cstdio>
//using namespace std;


#include "numlib.h"

#include "Estep.h"
#include "EM_algo.h"
#include "EM_init.h"
#include "hybrid.h"
#include "logdens.h"
#include "Segmentation_mean.h"
#include "Segmentation_mixt.h"
#include "Segmentation_ibp.h"
/**
 * Include Rwrappers.h after other includes to avoid
 * name clashes of R macros with stdlib components (i.e the length macro).
 **/
#include "Rwrappers.h"
#include "Rdefines.h"

using namespace cghseg;


SEXP sc_estep(SEXP logdensityR, SEXP rowslogdensityR, SEXP colslogdensityR, SEXP phiR)
{

  int rows_logdensity=*(INTEGER(rowslogdensityR));
  int cols_logdensity=*(INTEGER(colslogdensityR));

  numlib_matrix *logdensity=numlib_matrix_alloc(rows_logdensity,cols_logdensity);
  for (int index=0,col=0;col<cols_logdensity;col++)
	for (int row=0;row<rows_logdensity;row++,index++)
		numlib_matrix_set(logdensity,row,col,REAL(logdensityR)[index]);

  int length_phi=length(phiR);
  numlib_vector *phi=numlib_vector_alloc(length_phi);
  for (int i=0;i<length_phi;i++)
    numlib_vector_set(phi,i,REAL(phiR)[i]);

  numlib_matrix *tau=NULL;
  double lvinc;
  Estep(logdensity, phi,&tau, &lvinc);

  SEXP res;
  PROTECT(res=allocMatrix(REALSXP,tau->size1,tau->size2));

 for (int index=0,c=0;c<tau->size2;c++)
    for(int r=0;r<tau->size1;r++,index++)
	REAL(res)[index]=numlib_matrix_get(tau,r,c);

  UNPROTECT(1);

  numlib_vector_free(phi);
  numlib_matrix_free(logdensity);

  return res;
}


SEXP sc_logdens(SEXP xkR, SEXP PR, SEXP phiR)
{
  SEXP res;
  int length_xk=length(xkR);
  
  numlib_vector *xk=numlib_vector_calloc(length_xk);
  for (int i=0;i<xk->size;i++)
    numlib_vector_set(xk,i,REAL(xkR)[i]);

  int P=*(INTEGER(PR));

  int length_phi=length(phiR);

  numlib_vector *phi=numlib_vector_calloc(length_phi);
  for (int i=0;i<phi->size;i++)
    numlib_vector_set(phi,i,REAL(phiR)[i]);

  numlib_matrix *ld=logdens(xk,P,phi);

  PROTECT(res=allocMatrix(REALSXP,ld->size1,ld->size2));

 for (int index=0,c=0;c<ld->size2;c++)
    for(int r=0;r<ld->size1;r++,index++)
	REAL(res)[index]=numlib_matrix_get(ld,r,c);

  UNPROTECT(1);

  numlib_matrix_free(ld);
  numlib_vector_free(phi);
  numlib_vector_free(xk);

  return res;

  
}

SEXP sc_hybrid(SEXP xR, SEXP PR, SEXP KmaxR, SEXP lminR, SEXP lmaxR, SEXP vhR, SEXP fastR)
{
  SEXP res;
    
  int length_x=length(xR);
  numlib_vector *x=numlib_vector_calloc(length_x);
  for (int i=0;i<x->size;i++)
    numlib_vector_set(x,i,REAL(xR)[i]);


  int P=*(INTEGER(PR));
  int Kmax=*(INTEGER(KmaxR));
  int lmin=*(INTEGER(lminR));
  int lmax=*(INTEGER(lmaxR));
  int vhInt=*(LOGICAL(vhR));
  bool vh=true;
  if (!vhInt)
    vh=false;
  int fastInt = *(LOGICAL(fastR));
  bool fast=false;
  if (!fastInt)
    fast=true;

  numlib_matrix *Linc_res=0;
  param_struct *param_res=0;
  long param_res_size=0;

  hybrid(x,P,Kmax,lmin,lmax, vh,fast, &Linc_res,&param_res,&param_res_size);

  numlib_vector_free(x);
  

  PROTECT(res=allocVector(VECSXP,2));

  SEXP Linc_resR;
  PROTECT(Linc_resR=allocMatrix(REALSXP,Linc_res->size1,Linc_res->size2));
  for (int c=0;c<Linc_res->size2;c++)
    for (int r=0;r<Linc_res->size1;r++)
      REAL(Linc_resR)[r+c*Linc_res->size1]=numlib_matrix_get(Linc_res,r,c);
  numlib_matrix_free(Linc_res);
  SET_VECTOR_ELT(res,0,Linc_resR);

  SEXP param_resR;
  PROTECT(param_resR=allocVector(VECSXP,param_res_size));
  SET_VECTOR_ELT(res,1,param_resR);

  SEXP param_i;
  SEXP param_i_phi;
  SEXP param_i_rupt;
  SEXP param_i_names;
  int nprotect=0;
  for (int i=0;i<param_res_size;i++) {
    PROTECT(param_i=allocVector(VECSXP,2));
    nprotect++;
    SET_VECTOR_ELT(param_resR,i,param_i);

    PROTECT(param_i_names=allocVector(STRSXP,2));
    nprotect++;
    SET_STRING_ELT(param_i_names,0,mkChar("phi"));
    SET_STRING_ELT(param_i_names,1,mkChar("rupt"));
    setAttrib(param_i,R_NamesSymbol,param_i_names);

    if (param_res[i].phi) {
      PROTECT(param_i_phi=allocVector(REALSXP,param_res[i].phi->size));
      nprotect++;
      SET_VECTOR_ELT(param_i,0,param_i_phi);
      for (int j=0;j<param_res[i].phi->size;j++)
	REAL(param_i_phi)[j]=numlib_vector_get(param_res[i].phi,j);
      numlib_vector_free(param_res[i].phi);
    }

    if (param_res[i].rupt) {
      PROTECT(param_i_rupt=allocMatrix(REALSXP,param_res[i].rupt->size1,
				       param_res[i].rupt->size2));
      nprotect++;
      for (int c=0;c<param_res[i].rupt->size2;c++)
	for (int r=0;r<param_res[i].rupt->size1;r++)
	  REAL(param_i_rupt)[r+c*param_res[i].rupt->size1]=numlib_matrix_get(param_res[i].rupt,r,c);
      numlib_matrix_free(param_res[i].rupt);
      SET_VECTOR_ELT(param_i,1,param_i_rupt);
    }
  }
  
  delete[] param_res;

  SEXP res_names;
  PROTECT(res_names=allocVector(STRSXP,2));
  SET_STRING_ELT(res_names,0,mkChar("Linc"));
  SET_STRING_ELT(res_names,1,mkChar("param"));
  setAttrib(res,R_NamesSymbol,res_names);

  UNPROTECT(4+nprotect);

  return res;

}


SEXP sc_segmean(SEXP xR, SEXP lminR, SEXP lmaxR, SEXP KmaxR, SEXP vhR)
{ 
  SEXP res;
  int lengthx = length(xR);
  double *x = new double[lengthx];

  for (int i=0;i<lengthx;i++)
    x[i] = REAL(xR)[i];

  int lmin  = *(INTEGER(lminR));
  int lmax  = *(INTEGER(lmaxR));
  int Kmax  = *(INTEGER(KmaxR)); 
  bool vh   = true;
  int vhInt = *(LOGICAL(vhR));
  if (!vhInt)
    vh=false;

  Segmentation_mean Seg_mean(lengthx, Kmax, lmin, lmax, vh);
  Seg_mean.Init(x, vh);

  for (int k = 1; k < Kmax; k++)
    Seg_mean.DynProg(k);

  Seg_mean.Terminate();

  double *J_est = new double[Kmax];

  for (int k=0; k<Kmax;k++)
  {
    int t = Min(lengthx,(k+1)*lmax);
    J_est[k] = Seg_mean._D[k][t-(k+1)*lmin];
  }

  SEXP J_estR;
  SEXP t_estR;

  PROTECT(J_estR=allocVector(REALSXP,Kmax));

  for (int i=0;i<Kmax;i++)
    REAL(J_estR)[i]= J_est[i];

  PROTECT(t_estR=allocMatrix(REALSXP,Kmax,Kmax));

 for (int index=0,c=0;c<Kmax;c++)
    for(int r=0;r<Kmax;r++,index++)
      if (r>=c){
	REAL(t_estR)[index]= Seg_mean._BestBreaks[r][c]+1;
      } else {
	REAL(t_estR)[index]=0;
      }

  PROTECT(res=allocVector(VECSXP,2));
  SET_VECTOR_ELT(res,0,J_estR);
  SET_VECTOR_ELT(res,1,t_estR);

  SEXP res_names;
  PROTECT(res_names=allocVector(STRSXP,2));
  SET_STRING_ELT(res_names,0,mkChar("J.est"));
  SET_STRING_ELT(res_names,1,mkChar("t.est"));
  setAttrib(res,R_NamesSymbol,res_names);
   
  UNPROTECT(4);
  
  delete[] x;
  delete[] J_est;

  return res;
}




SEXP sc_segmixt(SEXP xR, SEXP lminR, SEXP lmaxR, SEXP KmaxR, SEXP phiR, SEXP PR)
{ 
  SEXP res;
  int lengthx = length(xR);
  double *x = new double[lengthx];

  for (int i=0;i<lengthx;i++)
    x[i] = REAL(xR)[i];

  int lmin  = *(INTEGER(lminR));
  int lmax  = *(INTEGER(lmaxR));
  int Kmax  = *(INTEGER(KmaxR)); 

  int length_phi=length(phiR);

  double *phi = new double[length_phi];
  for (int i=0;i<length_phi;i++)
    phi[i] = REAL(phiR)[i];
  
  int P=*(INTEGER(PR));

  Segmentation_mixt Seg_mixt(lengthx, Kmax, lmin, lmax, P);
  Seg_mixt.Init(x,phi);

  for (int k = 1; k < Kmax; k++)
    Seg_mixt.DynProg(k);

  Seg_mixt.Terminate();

  double *J_est = new double[Kmax];

  for (int k=0; k<Kmax;k++)
  {
    int t = Min(lengthx,(k+1)*lmax);
    J_est[k] = Seg_mixt._D[k][t-(k+1)*lmin];
  }


  SEXP J_estR;
  SEXP t_estR;

  PROTECT(J_estR=allocVector(REALSXP,Kmax));

  for (int i=0;i<Kmax;i++)
    REAL(J_estR)[i]= J_est[i];

  PROTECT(t_estR=allocMatrix(REALSXP,Kmax,Kmax));

  for (int index=0,c=0;c<Kmax;c++)
    for(int r=0;r<Kmax;r++,index++)
      if (r>=c){
	REAL(t_estR)[index]= Seg_mixt._BestBreaks[r][c]+1;
      } else {
	REAL(t_estR)[index]=0;
      }

  PROTECT(res=allocVector(VECSXP,2));
  SET_VECTOR_ELT(res,0,J_estR);
  SET_VECTOR_ELT(res,1,t_estR);

  SEXP res_names;
  PROTECT(res_names=allocVector(STRSXP,2));
  SET_STRING_ELT(res_names,0,mkChar("J.est"));
  SET_STRING_ELT(res_names,1,mkChar("t.est"));
  setAttrib(res,R_NamesSymbol,res_names);
   
  UNPROTECT(4);

  delete[] J_est;
  delete[] phi;  
  delete[] x;
  return res;
}


SEXP sc_segibp(SEXP ContrastVectR, SEXP KseqR, SEXP multiKmaxR)
{ 
  SEXP res;
  
  int     lengthContrastVect = length(ContrastVectR);
  int     M                  = length(KseqR);
  double *ContrastVect       = new double[lengthContrastVect];
  double ** Contrast         = new double*[M];
  int    *Kseq               = new int[M];
  int    multiKmax           = *(INTEGER(multiKmaxR));
  
  
  for (int i=0;i<lengthContrastVect;i++)
    ContrastVect[i] = REAL(ContrastVectR)[i];
  
  for (int m=0;m<M;m++)
    Kseq[m] = INTEGER(KseqR)[m];
  
  int a = 0;
  int b = 0;
  for (int m = 0; m < M ; m++)
    {
      Contrast[m] = new double[Kseq[m]];
      b = a + Kseq[m];
      
      for (int i = a ; i<b; i++)
        Contrast[m][i-a] = ContrastVect[i];
      a = b;
    }
  
  
  Segmentation_ibp Seg_ibp(M,Kseq,multiKmax);
  Seg_ibp.Init(Contrast);
  for (int m = 1; m < M; m++)
    Seg_ibp.DynProg(m);    
  Seg_ibp.Terminate(multiKmax);
  
 
  // dimensions verif ok
  double *J_est = new double[multiKmax-M+1];
  for (int i=0; i<(multiKmax-M+1); i++)
    J_est[i] = Seg_ibp._J[M-1][i];

  SEXP J_estR;
  SEXP t_estR;
  
  PROTECT(J_estR=allocVector(REALSXP,multiKmax-M+1));

  for (int i=0;i<multiKmax-M+1;i++)
    REAL(J_estR)[i]= J_est[i];
 
  PROTECT(t_estR=allocMatrix(REALSXP,M,multiKmax-M+1));

  // dimension verif ok
  for (int index=0,k=0;k<multiKmax-M+1;k++)
    for(int m=0;m<M;m++,index++)
	REAL(t_estR)[index]= Seg_ibp._BestBreaks[m][k]+1;
    
  PROTECT(res=allocVector(VECSXP,2));
  SET_VECTOR_ELT(res,0,J_estR);
  SET_VECTOR_ELT(res,1,t_estR);
  

  SEXP res_names;
  PROTECT(res_names=allocVector(STRSXP,2));
  SET_STRING_ELT(res_names,0,mkChar("J.est"));
  SET_STRING_ELT(res_names,1,mkChar("t.est"));
  setAttrib(res,R_NamesSymbol,res_names);
  
  UNPROTECT(4);
  
  delete[] ContrastVect;
  delete[] Kseq;
  for (int m = 0; m < M; m++)
    delete[] Contrast[m] ;
  delete[] Contrast;
  delete[] J_est;  


  return res;
}

SEXP sc_EMalgo(SEXP xR, SEXP phiR, SEXP ruptR, SEXP KR, SEXP PR, SEXP vhR)
{
  
  SEXP res;
  int lengthx = length(xR);
  long K      = *(INTEGER(KR));
  long P      = *(INTEGER(PR));
  bool vh     = true;
  int vhInt   = *(LOGICAL(vhR));
  if (!vhInt)
    vh=false;
  double lvinc = 0;
  int    empty = 0;
  int    dv    = 0;
  
  
  double *data = new double[lengthx];
  for (int i=0;i<length(xR);i++)
    data[i] = REAL(xR)[i];
  
  double *phi0 = new double[length(phiR)];
  for (int p=0;p<length(phiR);p++)
    phi0[p] = REAL(phiR)[p];

  numlib_matrix  *rupt=numlib_matrix_alloc(K,2);
  for (int index=0,col=0;col<2;col++)
    for (int row=0;row<K;row++,index++)
      numlib_matrix_set(rupt,row,col,REAL(ruptR)[index]);

  int *Bp = new int[K+1];
  for (int k=0;k<K;k++){
    Bp[k] = (int)numlib_matrix_get(rupt,k,0)-1;
  }
  Bp[K] = length(xR);
  
  EM_algo EMa(lengthx,K,P,vh);  
  EMa.Init(data,phi0,Bp);
  EMa.EM();


  SEXP phi_out;
  PROTECT(phi_out = allocVector(REALSXP,3*P));
  for (int p=0;p<3*P;p++)
    REAL(phi_out)[p] = EMa._phi[p];
  
  SEXP tau_out;
  PROTECT(tau_out = allocMatrix(REALSXP,K,P));
  for (int index=0,c=0;c<P;c++)
    for(int r=0;r<K;r++,index++)
      REAL(tau_out)[index]= EMa._tau[r][c];
  
  SEXP lv_out;
  SEXP empty_out;
  SEXP dv_out;
  
  PROTECT(lv_out=NEW_NUMERIC(1));
  NUMERIC_POINTER(lv_out)[0] = EMa._lvinc;
  PROTECT(empty_out=NEW_INTEGER(1));
  INTEGER_POINTER(empty_out)[0] = EMa._empty;
  PROTECT(dv_out=NEW_INTEGER(1));  
  INTEGER_POINTER(dv_out)[0] = EMa._dv;

  PROTECT(res=allocVector(VECSXP,5));
  SET_VECTOR_ELT(res,0,phi_out);
  SET_VECTOR_ELT(res,1,tau_out);
  SET_VECTOR_ELT(res,2,lv_out);
  SET_VECTOR_ELT(res,3,empty_out);
  SET_VECTOR_ELT(res,4,dv_out);

  SEXP res_names;
  PROTECT(res_names=allocVector(STRSXP,5));
  SET_STRING_ELT(res_names,0,mkChar("phi"));
  SET_STRING_ELT(res_names,1,mkChar("tau"));
  SET_STRING_ELT(res_names,2,mkChar("lvinc"));
  SET_STRING_ELT(res_names,3,mkChar("empty"));
  SET_STRING_ELT(res_names,4,mkChar("dv"));
  setAttrib(res,R_NamesSymbol,res_names);
    
  UNPROTECT(7);


  numlib_matrix_free(rupt);


   delete[] data;
   delete[] phi0;
   delete[] Bp;


  return res;
}

SEXP sc_EMinit(SEXP xR, SEXP ruptR, SEXP KR, SEXP PR, SEXP vhR)
{
  
  SEXP res;
  int lengthx = length(xR);
  long K      = *(INTEGER(KR));
  long P      = *(INTEGER(PR));
  bool vh     = true;
  int vhInt   = *(LOGICAL(vhR));
  if (!vhInt)
    vh=false;
  
  double *data = new double[lengthx];
  for (int i=0;i<length(xR);i++)
    data[i] = REAL(xR)[i];
  
  numlib_matrix  *rupt=numlib_matrix_alloc(K,2);
  for (int index=0,col=0;col<2;col++)
    for (int row=0;row<K;row++,index++)
      numlib_matrix_set(rupt,row,col,REAL(ruptR)[index]);

  int *Bp = new int[K+1];
  for (int k=0;k<K;k++){
    Bp[k] = (int)numlib_matrix_get(rupt,k,0)-1;
  }
  Bp[K] = length(xR);
  
  EM_init EMi(lengthx,K,P);  
  EMi.Init(data,Bp);
  EMi.CAH();
  EMi.compute_phi();
  
  PROTECT(res = allocVector(REALSXP,3*P));
  for (int p=0;p<3*P;p++)
    REAL(res)[p] = EMi._phi[p];
  
  numlib_matrix_free(rupt);
  UNPROTECT(1);
  
  delete[] data;
  delete[] Bp;
  return res;
}




