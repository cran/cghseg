#include "helperfuncs.h"

#include <cmath>
#include <cassert>
using namespace std;

#include "numlib.h"

namespace cghseg {

numlib_matrix * 
apply_basicfunc_numlib_matrix(numlib_matrix *mat, basicfunc bf) {

  numlib_matrix *res=numlib_matrix_calloc(mat->size1,mat->size2);
  for (int r=0;r<mat->size1;r++)
    for (int c=0;c<mat->size2;c++)
      numlib_matrix_set(res,r,c,bf(numlib_matrix_get(mat,r,c)));

  return res;
  
}

numlib_vector * 
apply_basicfunc_numlib_vector(numlib_vector *v, basicfunc bf) {

  numlib_vector *res=numlib_vector_calloc(v->size);
  for (int i=0;i<v->size;i++)
      numlib_vector_set(res,i,bf(numlib_vector_get(v,i)));

  return res;
  
}


numlib_vector *
copy_range_numlib_vector(numlib_vector *v, int min, int max)
{
  if (min>max) {
    int tmp=min;
    min=max;
    max=min;
  }

  numlib_vector *res=numlib_vector_calloc((max-min)+1);

  for (int i=min;i<=max;i++)
    numlib_vector_set(res,i-min,numlib_vector_get(v,i));

  return res;
}


numlib_vector *
build_range_numlib_vector(int min, int max)
{
  int inc=1;
  int size=(max-min)+1;

  if (min>max) {
    size=(min-max)+1;
    inc=-1;
  }

  numlib_vector *res=numlib_vector_calloc(max-min+1);
  for (int index=0,val=min;index<size;index++,val+=inc) {
    numlib_vector_set(res,index,val);
  }

  return res;
}

numlib_vector *
rowsum_numlib_matrix(numlib_matrix *mat)
{
  numlib_vector *vec=numlib_vector_calloc(mat->size1);

  for (int r=0;r<mat->size1;r++)
    for (int c=0;c<mat->size2;c++)
      numlib_vector_set(vec,r,numlib_vector_get(vec,r)+numlib_matrix_get(mat,r,c));

  return vec;
}

numlib_vector *
colsum_numlib_matrix(numlib_matrix *mat)
{
  numlib_vector *vec=numlib_vector_calloc(mat->size2);

    for (int c=0;c<mat->size2;c++)
      for (int r=0;r<mat->size1;r++)
	numlib_vector_set(vec,c,numlib_vector_get(vec,c)+numlib_matrix_get(mat,r,c));

  return vec;
}

double
sum_numlib_matrix(numlib_matrix *mat)
{
  double sum=0.0;

  for (int r=0;r<mat->size1;r++)
    for (int c=0;c<mat->size2;c++)
      sum+=numlib_matrix_get(mat,r,c);

  return sum;
}

double
sum_numlib_vector(numlib_vector *v)
{
  double sum=0.0;

  for (int i=0;i<v->size;i++)
    sum+=numlib_vector_get(v,i);

  return sum;
}

numlib_vector *
cumsum_numlib_matrix(numlib_matrix *m)
{
  numlib_vector *v=numlib_vector_calloc(m->size1*m->size2);

  int index=0;
  for (int r=0;r<m->size1;r++)
    for (int c=0;c<m->size2;c++) {
      if (r == 0 && c == 0) {
	numlib_vector_set(v,index,numlib_matrix_get(m,r,c));
      } else {
	numlib_vector_set(v,index,numlib_vector_get(v,index-1)+numlib_matrix_get(m,r,c));
      }
      index++;
    }
  return v;
}

numlib_vector *
cumsum_numlib_vector(numlib_vector *v)
{
  numlib_vector *res=numlib_vector_calloc(v->size);
  numlib_vector_set(res,0,numlib_vector_get(v,0));
  for (int index=1;index<v->size;index++)
    numlib_vector_set(res,index,numlib_vector_get(res,index-1)+numlib_vector_get(v,index));
  return res;
  
}

numlib_vector *
colmax_numlib_matrix(numlib_matrix *m)
{
  numlib_vector *res=numlib_vector_calloc(m->size2);
  for (int c=0;c<m->size2;c++) {
    double colmax=numlib_matrix_get(m,0,c);
    for (int r=1;r<m->size1;r++) {
      colmax=NUMLIB_MAX(colmax,numlib_matrix_get(m,r,c));
    }
    numlib_vector_set(res,c,colmax);
  }

  return res;
}

numlib_vector *
rowmax_numlib_matrix(numlib_matrix *m)
{
  numlib_vector *res=numlib_vector_calloc(m->size1);
  for (int r=0;r<m->size1;r++) {
    double rowmax=numlib_matrix_get(m,r,0);
    for (int c=1;c<m->size2;c++) {
      rowmax=NUMLIB_MAX(rowmax,numlib_matrix_get(m,r,c));
    }
    numlib_vector_set(res,r,rowmax);
  }

  return res;
}

numlib_vector *
colmin_numlib_matrix(numlib_matrix *m)
{
  numlib_vector *res=numlib_vector_calloc(m->size2);
  for (int c=0;c<m->size2;c++) {
    double colmin=numlib_matrix_get(m,0,c);
    for (int r=1;r<m->size1;r++) {
      colmin=NUMLIB_MIN(colmin,numlib_matrix_get(m,r,c));
    }
    numlib_vector_set(res,c,colmin);
  }

  return res;
}

numlib_vector *
rowmin_numlib_matrix(numlib_matrix *m)
{
  numlib_vector *res=numlib_vector_calloc(m->size1);
  for (int r=0;r<m->size1;r++) {
    double rowmin=numlib_matrix_get(m,r,0);
    for (int c=1;c<m->size2;c++) {
      rowmin=NUMLIB_MIN(rowmin,numlib_matrix_get(m,r,c));
    }
    numlib_vector_set(res,r,rowmin);
  }

  return res;
}

double
mean_numlib_vector(numlib_vector *v) {
  double mean=NUMLIB_NAN;
  if (v->size) {
    mean=0.0;
    for (int i=0;i<v->size;i++)
      mean+=numlib_vector_get(v,i);
    mean/=(1.0*v->size);
  }
  return mean;
  
}

double
var_numlib_vector(numlib_vector *v) {
  double var=NUMLIB_NAN;
  if (v->size) {
    var=0.0;
    double mean=mean_numlib_vector(v);
    for (int i=0;i<v->size;i++) {
      var+=(numlib_vector_get(v,i)*numlib_vector_get(v,i))-mean*mean;
    }
    if (v->size>1)
      var/=(v->size-1.0);
  }

  return var;

}

double
stddev_numlib_vector(numlib_vector *v) {

  double stddev=NUMLIB_NAN;
  if (v->size)
    stddev=sqrt(var_numlib_vector(v));
  
  return stddev;

}

numlib_vector *
diff_numlib_vector(numlib_vector *v) {

  numlib_vector *res=numlib_vector_calloc(v->size-1);
  for (int i=1; i<v->size;i++)
    numlib_vector_set(res,i-1,numlib_vector_get(v,i)-numlib_vector_get(v,i-1));

  return res;
}

int count_diffs_numlib_vectors(numlib_vector *v1, numlib_vector *v2)
{
  numlib_vector *shortest=v1;
  numlib_vector *longest=v2;
  
  if (shortest->size > longest->size) {
    numlib_vector *tmp=shortest;
    shortest=longest;
    longest=tmp;
  }

  int n=longest->size-shortest->size;

  for (int i=0;i<shortest->size;i++)
    if (numlib_vector_get(shortest,i) != numlib_vector_get(longest,i))
      n++;

  return n;
}

int
sum_diffs_numlib_matrices(numlib_matrix *m1, numlib_matrix *m2)
{
  assert(m1->size1 == m2->size1);
  assert(m1->size2 == m2->size2);

  int n=0;

  for (int r=0;r<m1->size1;r++)
    for (int c=0;c<m1->size2;c++)
      if (numlib_matrix_get(m1,r,c) != numlib_matrix_get(m2,r,c))
	n++;

  return n;
  
}

double checked_max_numlib_vector(numlib_vector *v)
{
  double maxval=NUMLIB_NEGINF;
  for (int index=0;index<v->size;index++)
    if (!isnan(numlib_vector_get(v,index)) && numlib_vector_get(v,index)>maxval)
      maxval=numlib_vector_get(v,index);


  return maxval;
}

int checked_max_index_numlib_vector(numlib_vector *v)
{
  double maxval=NUMLIB_NEGINF;
  int index=-1;

  for (int i=0;i<v->size;i++) {
    if ((!isnan(numlib_vector_get(v,i))) && numlib_vector_get(v,i)>maxval) {
      maxval=numlib_vector_get(v,i);
      index=i;
    }
  }

  return index;
}

double checked_min_numlib_vector(numlib_vector *v)
{
  double minval=NUMLIB_POSINF;
  for (int index=0;index<v->size;index++)
    if (!isnan(numlib_vector_get(v,index)) && numlib_vector_get(v,index)<minval)
      minval=numlib_vector_get(v,index);

  return minval;
}

int checked_min_index_numlib_vector(numlib_vector *v)
{
  double minval=NUMLIB_POSINF;
  int index=-1;

  for (int i=0;i<v->size;i++)
    if ((!isnan(numlib_vector_get(v,i))) && numlib_vector_get(v,i)<minval) {
      minval=numlib_vector_get(v,i);
      index=i;
    }

  return index;
}

}
