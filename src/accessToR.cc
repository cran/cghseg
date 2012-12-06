#include "colibri.h"
#include<R_ext/Arith.h>

extern "C" {

  void colibriR_c (double *profil, int *nbi, int *Kmaxi, double *mini, double *maxi, int *origine,
      double *cout_n){
    //char *test="azeaz";
    //cout_n[3]=3.0;
    //int i=0;
    //while(i < *Kmaxi){
    //printf("%f,", cout_n[i]);
    //i=i+1;
    //}
    //printf("\n", cout_n[i]);
    colibri_c (profil, nbi, Kmaxi, mini, maxi, origine, cout_n);
  }

  void meanRuptR_c(double * data, int* position, int*k, double* res){
    unsigned int j=0;
    for (unsigned int i=0; i<*k; i++){
        double sum = 0.;
        unsigned int count = 0;
        while(j<=(position[i]-1)){ // R to C++ -> -1
            //if (!isnan(data[j])){
            if (!R_IsNA(data[j])){
                count++;
                sum += data[j];
            }
            j++;
        }
        res[i] = sum/double(count);
    }
  }
  void meansqRuptR_c(double * data, int* position, int*k, double* res){
    unsigned int j=0;
    for (unsigned int i=0; i<*k; i++){
        double sum = 0.;
        unsigned int count = 0;
        while(j<=(position[i]-1)){ // R to C++ -> -1
            //if (!isnan(data[j])){
            if (!R_IsNA(data[j])){
                count++;
                sum += data[j]*data[j];
            }
            j++;
        }
        res[i] = sum/double(count);
    }
  }

}

