/**
#include <R.h>
#include <Rinternals.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
**/

#include <cstdlib>

#include "numlib.h"

extern "C" {
  void LinProgDyn(double *, int *, int *, double *, int *);
}


void LinProgDyn(double *sequence, int *lgSeq, int *nStep, double *res1, int *res2){
	/* Lecture des données */
	/* On stocke les données comme un liste de vecteurs */
	//printf("ligne 11\n");
	char c = 13;
	numlib_matrix_view matResult1 = numlib_matrix_view_array(res1, *nStep, *lgSeq);
	numlib_matrix_set_all(&matResult1.matrix, NUMLIB_POSINF);
	int i, j, k, l;
	
	/* on calcule le carré de chaque point */
	double *normSeq = (double*) malloc(* lgSeq * sizeof(double));
	i=0;
	while(i < *lgSeq){
		//fprintf(stderr, "Boucle ligne 21 : %d \n", i);
		normSeq[i]= sequence[i] * sequence[i];
		i++;
	}
	
	//printf("ligne 25\n");

	/* Initialisation chemin allant de 0 à i, i de 1 à n-1*/
	double SommeCarre, Somme, Poids;
	SommeCarre=normSeq[0];
	Somme=sequence[0];
	Poids=0.0;
	
	i=0;
	while(i < *lgSeq-1){
		//fprintf(stderr, "Boucle ligne 31 : %d \n", i);
		numlib_matrix_set(&matResult1.matrix, 0, i, Poids);
		res2[i] = 0;
		
		SommeCarre = SommeCarre+normSeq[i+1];
		Somme = Somme + sequence[i+1];
		Poids=SommeCarre - (Somme*Somme)/(i+2);
		i++;
	}
	numlib_matrix_set(&matResult1.matrix, 0, i, Poids);
	res2[i] = 0;

	/* Suite */
	int minim;
	double coutTraj;
	i=1;
	/* i : point de depart, j arrivée, k nombre de pas*/
	while(i < *lgSeq){
		if( i+1 <= *nStep){ minim=i+1; } 
		else{ minim=*nStep; }
		//printf("Sommet %d et Mini %d\n", i, minim);

		SommeCarre=normSeq[i];
		Somme=sequence[i];
		Poids=0.0;
		j=i+1;
		while(j < *lgSeq){
			
			k=1;
			while(k < minim){
				
				coutTraj = Poids + numlib_matrix_get(&matResult1.matrix, k-1, i-1);
				//printf("I= %d, K=%d, J=%d, Poids: %f, %f\n", i, k, j, Poids, coutTraj);
				if( coutTraj  < numlib_matrix_get(&matResult1.matrix, k, j-1)){
					numlib_matrix_set(&matResult1.matrix, k, j-1, coutTraj);
					res2[(*lgSeq)*k+j-1] = i;
				}
				k++;
			}
		
			SommeCarre = SommeCarre+normSeq[j];
			Somme = Somme + sequence[j];
			Poids=SommeCarre - (Somme*Somme)/(j-i+1);

				
	
		j++;
		}
		
		k=1;
		while(k < minim){
				
			coutTraj = Poids + numlib_matrix_get(&matResult1.matrix, k-1, i-1);
			//printf("I= %d, K=%d, J=%d, Poids: %f, %f\n", i, k, j, Poids, coutTraj);
			if( coutTraj  < numlib_matrix_get(&matResult1.matrix, k, j-1)){
				numlib_matrix_set(&matResult1.matrix, k, j-1, coutTraj);
				res2[(*lgSeq)*k+j-1] = i;
			}
			k++;
		}
		
		i++;
	}
	free(normSeq);

}
