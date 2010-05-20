#include "neighbors.h"

#include <cstdio>
#include <iostream>
using namespace std;

#include "EM_algo.h"
#include "Segmentation_mixt.h" 
#include "helperfuncs.h"

namespace cghseg {

void
neighbors(numlib_vector *x, long Kmax,
	  numlib_vector *L, long k, 
	  param_struct param[], long param_size,
	  long P, long lmin, long lmax, bool vh,
	  numlib_vector **L_res, 
	  param_struct *param_res[], long *param_res_size)
{

  if (*param_res) {
    for (int index=0;index<*param_res_size;index++) {
      if ((*param_res)[index].phi)
	numlib_vector_free((*param_res)[index].phi);
      if ((*param_res)[index].rupt) {
	numlib_matrix_free((*param_res)[index].rupt);
      }
    }
    delete[] *param_res;
  }

  *param_res=new param_struct[param_size];
  *param_res_size=param_size;
  for (int index=0;index<param_size;index++) {
    numlib_vector *phi=0;
    if (param[index].phi) {
      phi=numlib_vector_calloc(param[index].phi->size);
      numlib_vector_memcpy(phi,param[index].phi);
    }
    (*param_res)[index].phi=phi;

    numlib_matrix *rupt=0;
    if (param[index].rupt) {
      rupt=numlib_matrix_calloc(param[index].rupt->size1,
			     param[index].rupt->size2);
      numlib_matrix_memcpy(rupt,param[index].rupt);
    }
    (*param_res)[index].rupt=rupt;
  }

  if (*L_res)
    numlib_vector_free(*L_res);
  *L_res=numlib_vector_calloc(L->size);
  numlib_vector_memcpy(*L_res,L);



  //# left neighbor
  //a  = max(L[1:(k-1)])
  //K1 = which.max(L[1:(k-1)])
  double a=NUMLIB_NAN;
  int K1=-1;
  //  fprintf(stdout,"k=%d\n",k);
  if (1<=k-1) {
    numlib_vector *L_range=copy_range_numlib_vector(L,1-1,k-1-1);
    a=checked_max_numlib_vector(L_range);
    //    fprintf(stdout,"L_range\n");
    //    numlib_vector_fprintf(stdout,L_range,"%lf");
    K1=checked_max_index_numlib_vector(L_range)+1;
    numlib_vector_free(L_range);
  }

  //  fprintf(stdout,"a=%lf\nK1=%d\n",a,K1);
  numlib_vector *phi1=0;

  numlib_matrix *rupt1=0;
  numlib_vector *out_EM1_phi=0;
  numlib_matrix *out_EM1_tau=0;
  double out_EM1_lvinc=0.0;
  int out_EM1_empty=0;
  int out_EM1_dv=0;

  //if ( (is.null(K1)) || (a==-Inf) || (is.na(a) ) ){
  if ((K1<=0) || (numlib_isinf(a)==-1) || numlib_isnan(a)) {

    //  phi1          = rep(-Inf,3*P)
    phi1=numlib_vector_calloc(3*P);
    numlib_vector_set_all(phi1,NUMLIB_NEGINF);

    //  out.EM1$lvinc = list(lvinc=- Inf)
    out_EM1_lvinc=NUMLIB_NEGINF;

  //}    else { 
  } else {

    //phi1                     = param[[K1]]$phi     
    if (phi1)
      numlib_vector_free(phi1);
    phi1=numlib_vector_calloc(param[K1-1].phi->size);
    numlib_vector_memcpy(phi1,param[K1-1].phi);
    
    //G                        = Gmixt(x,lmin=lmin,phi1,P)
    numlib_matrix *t_est=0;

    compute_segmentation_mixt(x,k,lmin,lmax,phi1,P,NULL,&t_est);

    //rupt1      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
    if (rupt1) {
      numlib_matrix_free(rupt1);
    }
    rupt1=numlib_matrix_calloc(t_est->size2,2);
    numlib_matrix_set(rupt1,0,0,1.0);
    for (int index=1;index<=k-1;index++)
      numlib_matrix_set(rupt1,index,0,numlib_matrix_get(t_est,k-1,index-1)+1.0);
    numlib_vector *rowvec=numlib_vector_calloc(t_est->size2);
    numlib_matrix_get_row(rowvec,t_est,k-1);
    numlib_matrix_free(t_est);
    numlib_matrix_set_col(rupt1,1,rowvec);
    numlib_vector_free(rowvec);

    //out.EM1    = EM.algo(x,rupt1,P,phi1,vh)
    //EM_algo(x,rupt1,P,phi1,vh,&out_EM1_phi,&out_EM1_tau,&out_EM1_lvinc,&out_EM1_empty,&out_EM1_dv);

    compute_EM_algo(x,rupt1,P,phi1,vh,&out_EM1_phi,&out_EM1_tau,&out_EM1_lvinc,&out_EM1_empty,&out_EM1_dv);

    numlib_matrix_free(out_EM1_tau);

    // } #end else
  }
  numlib_vector_free(phi1);
  phi1=0;

  //# right neighbor 
  //a  = max(L[(k+1):Kmax])  
  //K2 = which.max(L[(k+1):Kmax])
  //K2 = K2 + k
  a=NUMLIB_NAN;
  int K2=-1;
  if (k+1<=Kmax) {
    numlib_vector *L_range=copy_range_numlib_vector(L,k+1-1,Kmax-1);
    a=checked_max_numlib_vector(L_range);

    K2=checked_max_index_numlib_vector(L_range)+1;
    numlib_vector_free(L_range);

    K2=K2+k;
  }

  numlib_vector *phi2=0;
  numlib_matrix *rupt2=0;
  numlib_vector *out_EM2_phi=0;
  numlib_matrix *out_EM2_tau=0;
  double out_EM2_lvinc=0.0;
  int out_EM2_empty=0;
  int out_EM2_dv=0;


  //if ( (K2==0) || (a==-Inf) || (is.na(a)) ){
  if ((K2<=0) || (numlib_isinf(a)==-1) || (numlib_isnan(a))) {




    //  phi2          = rep(-Inf,3*P)
    phi2=numlib_vector_calloc(3*P);
    numlib_vector_set_all(phi2,NUMLIB_NEGINF);

    //  out.EM2       = list(lvinc=-Inf)
    out_EM2_lvinc=NUMLIB_NEGINF;

    //}   else {
  } else {

    //phi2                     = param[[K2]]$phi
    phi2=numlib_vector_calloc(param[K2-1].phi->size);
    numlib_vector_memcpy(phi2,param[K2-1].phi);

    //G                        = Gmixt(x,lmin=lmin,phi2,P)
    numlib_matrix *t_est=0;

    compute_segmentation_mixt(x,k,lmin,lmax,phi2,P,NULL,&t_est);


    //rupt2      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
    if (rupt2)
      numlib_matrix_free(rupt2);
    rupt2=numlib_matrix_calloc(t_est->size2,2);
    numlib_matrix_set(rupt2,0,0,1.0);
    for (int index=1;index<=k-1;index++)
      numlib_matrix_set(rupt2,index,0,numlib_matrix_get(t_est,k-1,index-1)+1.0);
    numlib_vector *rowvec=numlib_vector_calloc(t_est->size2);
    numlib_matrix_get_row(rowvec,t_est,k-1);
    numlib_matrix_free(t_est);
    numlib_matrix_set_col(rupt2,1,rowvec);
    numlib_vector_free(rowvec);
    
    //out.EM2    = EM.algo(x,rupt2,P,phi2,vh)
//    EM_algo(x,rupt2,P,phi2,vh,&out_EM2_phi,&out_EM2_tau,&out_EM2_lvinc,&out_EM2_empty,&out_EM2_dv);
    
    compute_EM_algo(x,rupt2,P,phi2,vh,&out_EM2_phi,&out_EM2_tau,&out_EM2_lvinc,&out_EM2_empty,&out_EM2_dv);
    numlib_matrix_free(out_EM2_tau);
    out_EM2_tau=0;
  
    //} #end else
  }

  numlib_vector_free(phi2);
  phi2=0;
  

    
  //# choice of the best likelihood
  //a = which.max( c(out.EM1$lvinc, out.EM2$lvinc,  L[k]) ) 
  numlib_vector *maxvec=numlib_vector_calloc(3);
  numlib_vector_set(maxvec,0,out_EM1_lvinc);
  numlib_vector_set(maxvec,1,out_EM2_lvinc);
  numlib_vector_set(maxvec,2,numlib_vector_get(L,k-1));
  int aindex=checked_max_index_numlib_vector(maxvec)+1;
  //  fprintf(stdout,"maxvec=\n");
  //  numlib_vector_fprintf(stdout,maxvec,"%lf");
  //  fprintf(stdout,"max=%d\n",aindex);
    numlib_vector_free(maxvec);
  

  //# parameter update
  //if (length(a)==0) {
  if (aindex==-1) {

    //      L[k] = L
    numlib_vector_set(*L_res,k-1,numlib_vector_get(*L_res,0));

    //     param[[k]] = param
    if ((*param_res)[k-1].phi)
      numlib_vector_free((*param_res)[k-1].phi);
    (*param_res)[k-1].phi=0;
    if ((*param_res)[0].phi) {
      (*param_res)[k-1].phi=numlib_vector_calloc((*param_res)[0].phi->size);
      numlib_vector_memcpy((*param_res)[k-1].phi,(*param_res)[0].phi);
    }

    if ((*param_res)[k-1].rupt)
      numlib_matrix_free((*param_res)[k-1].rupt);
    (*param_res)[k-1].rupt=0;
    if ((*param_res)[0].rupt) {
      (*param_res)[k-1].rupt=numlib_matrix_calloc((*param_res)[0].rupt->size1,
					       (*param_res)[0].rupt->size2);
      numlib_matrix_memcpy((*param_res)[k-1].rupt,(*param_res)[0].rupt);
    }

    if (out_EM1_phi)
      numlib_vector_free(out_EM1_phi);
    if (out_EM2_phi)
      numlib_vector_free(out_EM2_phi);
    if (rupt1)
      numlib_matrix_free(rupt1);
    rupt1=0;
    if (rupt2)
      numlib_matrix_free(rupt2);
    rupt2=0;

  //} else {
  } else {
  //if (a==1){
    if (aindex==1) {

      //  param[[k]] = list(phi = out.EM1$phi, rupt = rupt1)
      if ((*param_res)[k-1].phi)
	numlib_vector_free((*param_res)[k-1].phi);
      (*param_res)[k-1].phi=out_EM1_phi;
      if (out_EM2_phi)
	numlib_vector_free(out_EM2_phi);
      out_EM2_phi=0;

      if ((*param_res)[k-1].rupt)
	numlib_matrix_free((*param_res)[k-1].rupt);
      (*param_res)[k-1].rupt=rupt1;
      if (rupt2)
	numlib_matrix_free(rupt2);
      rupt2=0;
      
      //  L[k] = out.EM1$lvinc}
      numlib_vector_set(*L_res,k-1,out_EM1_lvinc);
      
      
    }
    //if (a==2) {
    if (aindex==2) {
      //  param[[k]] = list(phi = out.EM2$phi,rupt = rupt2)
      if ((*param_res)[k-1].phi)
	numlib_vector_free((*param_res)[k-1].phi);
      (*param_res)[k-1].phi=out_EM2_phi;
      if (out_EM1_phi)
	numlib_vector_free(out_EM1_phi);
      out_EM1_phi=0;

      if ((*param_res)[k-1].rupt)
	numlib_matrix_free((*param_res)[k-1].rupt);
      (*param_res)[k-1].rupt=rupt2;
      if(rupt1)
	numlib_matrix_free(rupt1);
      rupt1=0;
      //  L[k] = out.EM2$lvinc}
      numlib_vector_set(*L_res,k-1,out_EM2_lvinc);
    }

    if (aindex==3) {
      if (out_EM1_phi)
	numlib_vector_free(out_EM1_phi);
      if (out_EM2_phi)
	numlib_vector_free(out_EM2_phi);
      if (rupt1)
	numlib_matrix_free(rupt1);
      rupt1=0;
      if (rupt2)
	numlib_matrix_free(rupt2);
      rupt2=0;
    }

  }
// invisible(list(L=L,param=param))  
    
}

}
