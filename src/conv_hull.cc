#include "conv_hull.h"

#include "helperfuncs.h"

namespace cghseg {

numlib_vector *
conv_hull(numlib_vector *J, numlib_vector *pen)
{


  //K = length(J)
  long K=J->size;

  //k = 1
  long k=1;

  //kv = c()
  numlib_vector *kv=numlib_vector_calloc(K);
  long kv_next=0;

  //pv = c()

  //  numlib_vector *pv=numlib_vector_calloc(K);
  //long pv_next=0;

  //while (k<K){
  while (k<K) {
    //pk = (J[(k+1):K]-J[k]) / (pen[k]-pen[(k+1):K])

    numlib_vector *Jrange=copy_range_numlib_vector(J,k+1-1,K-1);
    numlib_vector_add_constant(Jrange,-1.0*numlib_vector_get(J,k-1));

    numlib_vector *penrange=copy_range_numlib_vector(pen,k+1-1,K-1);
    numlib_vector_scale(penrange,-1.0);
    numlib_vector_add_constant(penrange,numlib_vector_get(pen,k-1));

    numlib_vector *pk=numlib_vector_alloc(Jrange->size);
    numlib_vector_memcpy(pk,Jrange);
    numlib_vector_free(Jrange);
    numlib_vector_div(pk,penrange);
    numlib_vector_free(penrange);

    //dm = which.max(pk)
    long dm=checked_max_index_numlib_vector(pk)+1;
    numlib_vector_free(pk);

    //kv = c(kv,k)
    numlib_vector_set(kv,kv_next,k);
    kv_next++;

    //k  = k + dm 
    k=k+dm;

    //}
  }
 
  //kv = c(kv,K)
  numlib_vector_set(kv,kv_next,K);
  kv_next++;

  //invisible(kv=kv)
  numlib_vector *kv_ret=numlib_vector_alloc(kv_next);
  for (int i=0;i<kv_next;i++)
    numlib_vector_set(kv_ret,i,numlib_vector_get(kv,i));

  numlib_vector_free(kv);

  return kv_ret;

}

}
