#include "colibri.h"

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
}

