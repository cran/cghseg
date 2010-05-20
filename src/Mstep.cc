#include "Mstep.h"

#include <cmath>
#include <iostream>
using namespace std;

/**
#include <gsl/numlib_blas.h>
#include <gsl/numlib_math.h>
#include <gsl/numlib_permutation.h>
#include <gsl/numlib_sort_vector.h>
**/

#include "helperfuncs.h"
#include "repmat.h"

namespace cghseg {

numlib_vector *
Mstep(numlib_vector *x, numlib_matrix *rupt, numlib_matrix *tau, numlib_vector *phi,
      bool vh)
{
  //K = nrow(tau)
  int K=tau->size1;

  //P = ncol(tau)
  int P=tau->size2;

  // m = matrix(nrow=P,ncol=1)
  numlib_vector *m=numlib_vector_calloc(P);
  
  //s = matrix(nrow=P,ncol=1)
  numlib_vector *s=numlib_vector_calloc(P);

  //prop = matrix(nrow=P,ncol=1)
  numlib_vector *prop=numlib_vector_calloc(P);

  //Yk = apply(rupt,1,FUN=function(y) sum(x[y[1]:y[2]]))
  numlib_vector *Yk=numlib_vector_calloc(rupt->size1);
  for (int r=0;r<rupt->size1;r++) {
    int y1=(int)numlib_matrix_get(rupt,r,1-1);
    int y2=(int)numlib_matrix_get(rupt,r,2-1);
    if (y1>y2) {
      int tmp_y=y1;
      y1=y2;
      y2=tmp_y;
    }
    double sum=0.0;
    for (int index=y1; index<=y2;index++)
      sum+=numlib_vector_get(x,index-1);
    numlib_vector_set(Yk,r,sum);
  }

  

  //nk = rupt[,2]-rupt[,1]+1
  numlib_vector *nk=numlib_vector_calloc(rupt->size1);

  numlib_vector *rupt2=numlib_vector_calloc(rupt->size1);
  numlib_matrix_get_col(rupt2,rupt,2-1);

  numlib_vector *rupt1=numlib_vector_calloc(rupt->size1);
  numlib_matrix_get_col(rupt1,rupt,1-1);

  numlib_vector_memcpy(nk,rupt2);
  numlib_vector_sub(nk,rupt1);
  numlib_vector_add_constant(nk,+1.0);

  numlib_vector_free(rupt1);
  numlib_vector_free(rupt2);

  //n  = sum(nk)
  double n=sum_numlib_vector(nk);

  //#homogeneous variances

  if (vh==true) {

    //for (p in (1:P)){
    for (int p=1;p<=P;p++) {
      
      // np    = nk %*% tau[,p]
      numlib_vector *taup=numlib_vector_calloc(tau->size1);
      numlib_matrix_get_col(taup,tau,p-1);
      double np=0.0;
      numlib_blas_ddot(nk,taup,&np);
      
      //m[p]  = Yk %*% tau[,p]/np
      numlib_vector *taup_div_np=numlib_vector_calloc(taup->size);
      numlib_vector_memcpy(taup_div_np,taup);
      numlib_vector_scale(taup_div_np,1.0/np);
      double prod=0.0;
      numlib_blas_ddot(Yk,taup_div_np,&prod);
      numlib_vector_free(taup_div_np);
      
      numlib_vector_set(m,p-1,prod);
      
      //tmp   = tau[,p] * ( apply(rupt,1,FUN=function(y) sum(    (x[y[1]:y[2]]-m[p])^2   )))  
      numlib_vector *res_apply=numlib_vector_calloc(taup->size);
      for (int r=0;r<rupt->size1;r++) {
	int y1=(int)numlib_matrix_get(rupt,r,1-1);
	int y2=(int)numlib_matrix_get(rupt,r,2-1);
	if (y1>y2) {
	  int tmp_y=y1;
	  y1=y2;
	  y2=tmp_y;
	}
	numlib_vector *diff_vec=copy_range_numlib_vector(x,y1-1,y2-1);
	numlib_vector_add_constant(diff_vec,-1.0*numlib_vector_get(m,p-1));
	numlib_vector_mul(diff_vec,diff_vec);
	double sum=sum_numlib_vector(diff_vec);
	numlib_vector_free(diff_vec);
	numlib_vector_set(res_apply,r,sum);
      }
      numlib_vector *tmp=numlib_vector_calloc(taup->size);
      numlib_vector_memcpy(tmp,taup);
      numlib_vector_mul(tmp,res_apply);
      numlib_vector_free(res_apply);
      
      //s[p]=sum(tmp)
      double sum=sum_numlib_vector(tmp);
      numlib_vector_free(tmp);
      numlib_vector_set(s,p-1,sum);
      
      numlib_vector_free(taup);
      //} for p
    }
    
    numlib_vector_free(Yk);
    numlib_vector_free(nk);
    
    
    //s = repmat(sum(s)/n,P,1)
    double sum=sum_numlib_vector(s)/(1.0*n);
    for (int index=0;index<s->size;index++)
      numlib_vector_set(s,index,sum);

  } else {
    
    numlib_vector *taup=0;
    numlib_vector *taup_div_np=0;
    numlib_vector *tmp=0;
    numlib_vector *ruptk1_k2=0;
    //for (p in (1:P)){
    for (int p=1;p<=P;p++) {
      
      // np    = nk %*% tau[,p]
      if (taup)
	numlib_vector_free(taup);
      taup=numlib_vector_calloc(tau->size1);
      numlib_matrix_get_col(taup,tau,p-1);
      double np=0.0;
      numlib_blas_ddot(nk,taup,&np);
      
      //m[p]  = Yk %*% tau[,p]/np
      if (taup_div_np)
	numlib_vector_free(taup_div_np);
      taup_div_np=numlib_vector_calloc(taup->size);
      numlib_vector_memcpy(taup_div_np,taup);
      numlib_vector_scale(taup_div_np,1.0/np);
      double prod=0.0;
      numlib_blas_ddot(Yk,taup_div_np,&prod);
      numlib_vector_free(taup_div_np);
      taup_div_np=0;

      numlib_vector_set(m,p-1,prod);

      //tmp = rep(Inf,K)
      if (tmp)
	numlib_vector_free(tmp);
      tmp=numlib_vector_calloc(K);
      numlib_vector_set_all(tmp,NUMLIB_POSINF);

      //for (k in (1:K)) {
      for (int k=1;k<=K;k++) {

	// tmp[k] = tau[k,p] *  sum ( (x[(rupt[k,2]:rupt[k,1])]-m[p])^2)
	int ruptk2=(int)numlib_matrix_get(rupt,k-1,2-1);
	int ruptk1=(int)numlib_matrix_get(rupt,k-1,1-1);
	if (ruptk1_k2)
	  numlib_vector_free(ruptk1_k2);
	ruptk1_k2=copy_range_numlib_vector(x,ruptk2-1,ruptk1-1);
	numlib_vector_add_constant(ruptk1_k2,-1.0*numlib_vector_get(m,p-1));
	numlib_vector_mul(ruptk1_k2,ruptk1_k2);
	double sum=sum_numlib_vector(ruptk1_k2);
	numlib_vector_free(ruptk1_k2);
	ruptk1_k2=0;
	numlib_vector_set(tmp,k-1,numlib_matrix_get(tau,k-1,p-1)*sum);

      //}
      }

      // s[p] = sum(tmp)/np
      numlib_vector_set(s,p-1,sum_numlib_vector(tmp)/np);
      numlib_vector_free(tmp);
      tmp=0;
      //}
      
    }
      
  }





  //s = sqrt(s)
  numlib_vector *tmp_s=apply_basicfunc_numlib_vector(s,sqrt);
  numlib_vector_free(s);
  s=tmp_s;

  //prop = apply(tau,2,sum)/K
  numlib_vector *colsum=colsum_numlib_matrix(tau);
  numlib_vector_scale(colsum,1.0/K);
  numlib_vector_memcpy(prop,colsum);
  numlib_vector_free(colsum);

  /**
   * !!! WARNING !!!
   * Due to numerical precision limits, the sort order of m
   * in the C++ code can be different form the order in the R code.
  **/

  //b    = order(m)
  numlib_permutation *b=numlib_permutation_calloc(m->size);
  numlib_sort_vector_index(b,m);

  //m    = sort(m)
  numlib_sort_vector(m);

  //s    = s[b]
  tmp_s=numlib_vector_calloc(s->size);
  for (int index=0; index<s->size;index++)
    numlib_vector_set(tmp_s,index,numlib_vector_get(s,numlib_permutation_get(b,index)));
  numlib_vector_free(s);
  s=tmp_s;

  //prop = prop[b]
  numlib_vector *tmp_prop=numlib_vector_alloc(prop->size);
  for (int index=0;index<prop->size;index++)
    numlib_vector_set(tmp_prop,index,numlib_vector_get(prop,numlib_permutation_get(b,index)));
  numlib_vector_free(prop);
  prop=tmp_prop;

  numlib_permutation_free(b);

  // phi  = c(m,s,prop)
  numlib_vector *phi_res=numlib_vector_calloc(m->size+s->size+prop->size);
  for (int index=0;index<m->size;index++)
    numlib_vector_set(phi_res,index,numlib_vector_get(m,index));
  for (int index=0;index<s->size;index++)
    numlib_vector_set(phi_res,m->size+index,numlib_vector_get(s,index));
  for (int index=0;index<prop->size;index++)
    numlib_vector_set(phi_res,m->size+s->size+index,numlib_vector_get(prop,index));

  numlib_vector_free(prop);
  numlib_vector_free(s);
  numlib_vector_free(m);

  return phi_res;
}

}
