#include "repmat.h"

namespace cghseg {

numlib_matrix *
repmat(numlib_matrix *a, int n, int m)
{
  numlib_matrix *res=numlib_matrix_alloc(n*a->size1,m*a->size2);
  int i,j,r,c;

  for (i=0;i<n;i++) 
    for (j=0;j<m;j++)
      for (r=0;r<a->size1;r++)
	for (c=0;c<a->size2;c++)
	  numlib_matrix_set(res,i*a->size1+r,j*a->size2+c,numlib_matrix_get(a,r,c));

  return res;
}

}
