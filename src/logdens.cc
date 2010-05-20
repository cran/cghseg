#include <cmath>
using namespace std;

#include "logdens.h"

namespace cghseg {

numlib_matrix *
logdens(numlib_vector *xk, int P, numlib_vector *phi) {

  static double sqrt_2_pi=sqrt(2*M_PI);

  int i,p;
  double ld,sp,mp,sum,sumterm;

  int m_index=0;
  int s_index=P;
  int nk=xk->size;

  numlib_matrix *tmp=numlib_matrix_alloc(1,P);
  
  for (p=0;p<P;p++) {
    sp=numlib_vector_get(phi,s_index+p);
    mp=numlib_vector_get(phi,m_index+p);

    sum=0.0;
    for (i=0;i<xk->size;i++) {
      sumterm=numlib_vector_get(xk,i)-mp;
      sum+=sumterm*sumterm;
    }
    ld=-nk*log(sqrt_2_pi*sp)-0.5*sum/(sp*sp);

    numlib_matrix_set(tmp,0,p,ld);
  }
  
  return tmp;
}

}
