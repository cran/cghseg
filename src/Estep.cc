#include "Estep.h"

#include <cmath>
#include <iostream>
using namespace std;

#include "helperfuncs.h"
#include "numlib.h"
#include "repmat.h"

namespace cghseg {

void
Estep(numlib_matrix *logdensity, numlib_vector *phi,
      numlib_matrix **tau_res, double *lvinc_res)
{

  //K = nrow(logdensity)
  int K=logdensity->size1;

  //P = ncol(logdensity)
  int P=logdensity->size2;

  //tau     = repmat( t (log( phi[(2*P+1):(3*P)] )),K,1)+logdensity
  int cmin=2*P+1;
  int cmax=3*P;
  numlib_vector *logphi=copy_range_numlib_vector(phi,cmin-1,cmax-1);
  numlib_vector *logphi_tmp=apply_basicfunc_numlib_vector(logphi,log);
  numlib_vector_free(logphi);
  logphi=logphi_tmp;

  numlib_matrix *logphi_mat=numlib_matrix_alloc(1,logphi->size);
  numlib_matrix_set_row(logphi_mat,0,logphi);
  numlib_vector_free(logphi);

  numlib_matrix *logphi_repmat=repmat(logphi_mat,K,1);
  numlib_matrix_free(logphi_mat);

  numlib_matrix *tau=numlib_matrix_calloc(logphi_repmat->size1,logphi_repmat->size2);
  numlib_matrix_memcpy(tau,logphi_repmat);
  numlib_matrix_free(logphi_repmat);
  numlib_matrix_add(tau,logdensity);


  //tau_max = apply(tau,1,max)
  numlib_vector *tau_max=rowmax_numlib_matrix(tau);

  //tau     = exp(tau-repmat(tau_max,1,P))
  numlib_matrix *tau_max_mat=numlib_matrix_calloc(tau_max->size,1);
  numlib_matrix_set_col(tau_max_mat,0,tau_max);
  numlib_matrix *tau_max_repmat=repmat(tau_max_mat,1,P);
  numlib_matrix_free(tau_max_mat);
  numlib_matrix_sub(tau,tau_max_repmat);
  numlib_matrix_free(tau_max_repmat);
  numlib_matrix *tau_tmp=apply_basicfunc_numlib_matrix(tau,exp);
  numlib_matrix_free(tau);
  tau=tau_tmp;

  //lvinc  = sum(log( apply(tau,1,sum)) + tau_max)
  numlib_vector *tau_rowsum=rowsum_numlib_matrix(tau);
  numlib_vector *tau_rowsumtmp=apply_basicfunc_numlib_vector(tau_rowsum,log);
  numlib_vector_free(tau_rowsum);
  tau_rowsum=tau_rowsumtmp;
  numlib_vector_add(tau_rowsum,tau_max);
  double lvinc=sum_numlib_vector(tau_rowsum);
  numlib_vector_free(tau_max);
  numlib_vector_free(tau_rowsum);

  //tau     = tau / repmat( apply(tau,1,sum) ,1,P)
  tau_rowsum=rowsum_numlib_matrix(tau);
  numlib_matrix *tau_rowsum_mat=numlib_matrix_calloc(tau_rowsum->size,1);
  numlib_matrix_set_col(tau_rowsum_mat,0,tau_rowsum);
  numlib_vector_free(tau_rowsum);
  numlib_matrix *tau_rowsum_repmat=repmat(tau_rowsum_mat,1,P);
  numlib_matrix_free(tau_rowsum_mat);
  numlib_matrix_div_elements(tau,tau_rowsum_repmat);
  numlib_matrix_free(tau_rowsum_repmat);
  
  //invisible(list(tau,lvinc))
  if (*tau_res)
    numlib_matrix_free(*tau_res);
  *tau_res=tau;

  *lvinc_res=lvinc;
}

}
