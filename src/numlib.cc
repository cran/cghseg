#include "numlib.h"

#ifndef USE_GSL

#include <cstdlib>
#include <cstring>

void 
numlib_blas_ddot(numlib_vector *v1,numlib_vector *v2,double *r)
{
  *r=0.0;
  double *d1=v1->data;
  double *d2=v2->data;
  for (int index=0;index<v1->size;index++,d1++,d2++)
    *r+=(*d1)*(*d2);
}



int
numlib_matrix_add(numlib_matrix *m1, numlib_matrix *m2)
{

  double *e1=m1->data;
  double *e2=m2->data;
  int etot=m1->size1*m2->size2;
  for (int i=0;i<etot;i++,e1++,e2++)
    *e1=*e1+*e2;
  return 0;
  
}

numlib_matrix *numlib_matrix_alloc(int size1,int size2)
{
  numlib_matrix *res=(numlib_matrix *)malloc(sizeof(numlib_matrix));
  res->size1=size1;
  res->size2=size2;
  res->data=(double *)malloc(size1*size2*sizeof(double));
  return res;
}

numlib_matrix *
numlib_matrix_calloc(int size1, int size2) 
{
  numlib_matrix *res=(numlib_matrix *)malloc(sizeof(numlib_matrix));
  res->size1=size1;
  res->size2=size2;
  res->data=(double *)calloc(size1*size2,sizeof(double));
  return res;
}

int
numlib_matrix_div_elements(numlib_matrix *m1, numlib_matrix *m2) 
{
  double *e1=m1->data;
  double *e2=m2->data;
  int etot=m1->size1*m2->size2;
  for (int i=0;i<etot;i++,e1++,e2++)
    *e1=*e1/(*e2);
  return 0;
}

void
numlib_matrix_free(numlib_matrix *m) 
{
  free(m->data);
  free(m);
}

double
numlib_matrix_get(numlib_matrix *m,int i,int j) 
{
  return m->data[i*m->size2+j];
}

void
numlib_matrix_get_col(numlib_vector *c,numlib_matrix *m,int i)
{
  double *dest=c->data;
  double *src=m->data+i;
  for (int index=0;index<m->size1;index++,src+=m->size2,dest++)
    *dest=*src;
  
}

void 
numlib_matrix_get_row(numlib_vector *r,numlib_matrix *m,int i)
{
  double *dest=r->data;
  double *src=m->data+i*m->size2;
  memcpy(dest,src,sizeof(double)*m->size2);
}

void
numlib_matrix_memcpy(numlib_matrix *m1,numlib_matrix *m2)
{
  memcpy(m1->data,m2->data,sizeof(double)*m2->size1*m2->size2);
}

void
numlib_matrix_set(numlib_matrix *m,int i1,int i2,double v)
{
  m->data[i1*m->size2+i2]=v;
}

void
numlib_matrix_set_all(numlib_matrix *m,double v)
{
  int nvalues=m->size1*m->size2;
  double *dest=m->data;
  for (int index=0;index<nvalues;index++,dest++)
    *dest=v;
}

void
numlib_matrix_set_col(numlib_matrix *m,int i, numlib_vector * v)
{
  double *src=v->data;
  double *dest=m->data+i;
  for (int index=0;index<m->size1;index++,dest+=m->size2,src++)
    *dest=*src;
}

void
numlib_matrix_set_row(numlib_matrix *m,int i,numlib_vector * v)
{
  double *src=v->data;
  double *dest=m->data+i*m->size2;
  memcpy(dest,src,sizeof(double)*m->size2);

}

int
numlib_matrix_sub(numlib_matrix *m1,numlib_matrix *m2)
{

  double *e1=m1->data;
  double *e2=m2->data;
  int etot=m1->size1*m2->size2;
  for (int i=0;i<etot;i++,e1++,e2++)
    *e1=*e1-*e2;
  return 0;

}

numlib_matrix_view
numlib_matrix_view_array(double *a,int size1,int size2)
{
  numlib_matrix_view v;
  v.matrix.size1=size1;
  v.matrix.size2=size2;
  v.matrix.data=a;
  return v;
}

numlib_permutation *
numlib_permutation_calloc(int n)
{
  numlib_permutation *res=(numlib_permutation *)calloc(1,sizeof(numlib_permutation));
  res->size=n;
  res->data=(int *)malloc(n*sizeof(int));
  for (int i=0;i<n;i++)
    res->data[i]=i;
  return res;
}

void
numlib_permutation_free(numlib_permutation *p)
{
  free(p->data);
  free(p);
}

int
numlib_permutation_get(numlib_permutation *p,int i)
{
  return p->data[i];
}


typedef struct {
  int index;
  double val;
} indexed_double;


static int cmpidxdbl(const void *v1, const void *v2) {
  double dbl1=((indexed_double *)v1)->val;
  double dbl2=((indexed_double *)v2)->val;
  int res=0;
  res=(dbl1<dbl2)?-1:res;
  res=(dbl1>dbl2)?1:res;
  return res;
}

void
numlib_sort_vector_index(numlib_permutation *p,numlib_vector *v)
{
  indexed_double *tab=(indexed_double *)malloc(v->size*sizeof(indexed_double));
  for (int index=0;index<v->size;index++) {
    tab[index].index=index;
    tab[index].val=v->data[index];
  }

  qsort(tab,v->size,sizeof(indexed_double),cmpidxdbl);

  for (int index=0;index<v->size;index++) {
    p->data[index]=tab[index].index;
  }

  free(tab);
}

static int cmpdbl(const void *v1, const void *v2) {
  double dbl1=*(double *)v1;
  double dbl2=*(double *)v2;
  int res=0;
  res=(dbl1<dbl2)?-1:res;
  res=(dbl1>dbl2)?1:res;
  return res;
}

void
numlib_sort_vector(numlib_vector *v)
{
  qsort(v->data,v->size,sizeof(double),cmpdbl);
}

int
numlib_vector_add(numlib_vector *v1,numlib_vector *v2)
{
  double *data1=v1->data;
  double *data2=v2->data;
  for (int index=0;index<v1->size;index++,data1++,data2++)
    *data1+=*data2;
  return 0;
}

int
numlib_vector_add_constant(numlib_vector *v,double c)
{
  double *data=v->data;
  for (int index=0;index<v->size;index++,data++)
    *data+=c;
  return 0;
}

numlib_vector *
numlib_vector_alloc(int n)
{

  numlib_vector *res=(numlib_vector *)malloc(sizeof(numlib_vector));
  res->size=n;
  res->data=(double *)malloc(n*sizeof(double));
  return res;
}

numlib_vector *
numlib_vector_calloc(int n)
{

  numlib_vector *res=(numlib_vector *)calloc(1,sizeof(numlib_vector));
  res->size=n;
  res->data=(double *)calloc(n,sizeof(double));
  return res;
 
}

int
numlib_vector_div(numlib_vector *v1, numlib_vector *v2)
{
  double *data1=v1->data;
  double *data2=v2->data;
  for (int index=0;index<v1->size;index++,data1++,data2++)
    *data1/=*data2;
  return 0;

}

void
numlib_vector_free(numlib_vector *v)
{
  free(v->data);
  free(v);
}

double
numlib_vector_get(numlib_vector *v,int i)
{
  return v->data[i];
}

void 
numlib_vector_memcpy(numlib_vector *v1, numlib_vector *v2)
{
  memcpy(v1->data,v2->data,sizeof(double)*v1->size);

}

int
numlib_vector_mul(numlib_vector *v1,numlib_vector *v2)
{
  double *data1=v1->data;
  double *data2=v2->data;
  for (int index=0;index<v1->size;index++,data1++,data2++)
    *data1*=*data2;
  return 0;
}

int
numlib_vector_scale(numlib_vector *v, double c)
{
  double *data=v->data;
  for (int index=0;index<v->size;index++,data++)
    *data*=c;
  return 0;
}

void 
numlib_vector_set(numlib_vector *v, int i, double x)
{
  v->data[i]=x;
}

void
numlib_vector_set_all(numlib_vector *v, double x)
{
    double *data=v->data;
  for (int index=0;index<v->size;index++,data++)
    *data=x;
}

int
numlib_vector_sub(numlib_vector *v1, numlib_vector *v2)
{
  double *data1=v1->data;
  double *data2=v2->data;
  for (int index=0;index<v1->size;index++,data1++,data2++)
    *data1-=*data2;
  return 0;
}

int
numlib_vector_int_get(numlib_vector_int *v,int i)
{
  return v->data[i];
}

double
numlib_sf_exp(double x)
{
  return exp(x);
}

double
numlib_sf_log(double x)
{
  return log(x);
}

double
numlib_pow_2(double x)
{
  return pow(x,2.0);
}

#endif
