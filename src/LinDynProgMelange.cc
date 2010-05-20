#include <cmath>
#include <cstdlib>

#include "numlib.h"
using namespace std;

extern "C" {

  void LinProgDynMelange(double *, int *, int *, double *, int *, int *, double *, double *, double *);
}


void LogValue(double *logValue, double *Value, double constante, int *nbClasse){
	int i=0;
	while( i < *nbClasse){
		logValue[i] = constante * numlib_sf_log (Value[i]);
		i=i+1;
	}
}


void InverseValue(double *InverseValue, double *Value, double constante, int *nbClasse){
	int i=0;
	while( i < *nbClasse){
		InverseValue[i] = constante / Value[i]; 
		i=i+1;
	}
}

void MiseAJour(double y, double *currentLogClasse, double *logVarianceSur2, double *inverseVarianceSur2, double *moyennes, int *nbClasse){
	int i=0;
	
	while( i < *nbClasse){
		//printf("logCurrent : %f, logVar : %f, inverseVar : %f, moyenne %f \n", currentLogClasse[i], 
		//logVarianceSur2[i], inverseVarianceSur2[i], moyennes[i]);
		currentLogClasse[i] = currentLogClasse[i]	- 0.5 * M_LNPI -0.5 * M_LN2 - logVarianceSur2[i] 
												  - inverseVarianceSur2[i] * numlib_pow_2(y - moyennes[i]) ;
		i=i+1;
	}

}

double getMax(double *Values, int *nbClasse){
	int i=0;
	double max=NUMLIB_NEGINF;
	while( i < *nbClasse){
		if(Values[i] > max){
			max= Values[i];		
		}
		i=i+1;
	}
	return(max);
}

double ComputeCost(double *currentLogClasse, int *nbClasse){
	double Max = getMax(currentLogClasse, nbClasse);
	double Somme=0.0;
	int i=0;
	
	while( i < *nbClasse){
		//Somme = Somme + numlib_sf_exp (currentLogClasse[i] - Max);
		Somme = Somme + exp(currentLogClasse[i] - Max); // Add
		//fprintf(stderr,"%f, %f, %f \n", currentLogClasse[i], Max, Somme);
		i=i+1;
	}
	
	//return(numlib_sf_log (Somme) + Max);
	return(log (Somme) + Max); // Add
}


void LinProgDynMelange(double *sequence, int *lgSeq, int *nStep, double *res1, int *res2, int *nbClasse, 
		double *moyennes, double *variances, double *proportions){
	/* Lecture des données */
	/* On stocke les données comme un liste de vecteurs */
	//printf("ligne 11\n");
	char c = 13;
	numlib_matrix_view matResult1 = numlib_matrix_view_array(res1, *nStep, *lgSeq);
	numlib_matrix_set_all(&matResult1.matrix, NUMLIB_POSINF);
	int i, j, k, l, iClasse;

	
	//printf("ligne 25\n");

	/* Initialisation chemin allant de 0 à i, i de 1 à n-1*/
	double *inverseVarianceSur2 = (double*) malloc(* nbClasse * sizeof(double));
	InverseValue(inverseVarianceSur2, variances, 0.5, nbClasse);
	double *logVarianceSur2 = (double*) malloc(* nbClasse * sizeof(double));
	LogValue(logVarianceSur2, variances, 0.5, nbClasse);

	double *currentLogClasse = (double*) malloc(* nbClasse * sizeof(double));
	LogValue(currentLogClasse, proportions, 1.0, nbClasse);
	MiseAJour(sequence[0], currentLogClasse, logVarianceSur2, inverseVarianceSur2, moyennes, nbClasse);
	
	double Poids;
	i=0;
	while(i < *lgSeq-1){
		//fprintf(stderr, "Tour : %d \n", i);
		numlib_matrix_set(&matResult1.matrix, 0, i, -ComputeCost(currentLogClasse, nbClasse));
		res2[i] = 0;
		MiseAJour(sequence[i+1], currentLogClasse, logVarianceSur2, inverseVarianceSur2, moyennes, nbClasse);
		
		i++;
	}
	numlib_matrix_set(&matResult1.matrix, 0, i, -ComputeCost(currentLogClasse, nbClasse));
	res2[i] = 0;

	/* Suite */
	int minim;
	double coutTraj;
	i=1;
	/* i : point de depart, j arrivée, k nombre de pas*/
	 while(i < *lgSeq){
		//fprintf(stderr, "%c Node :   %d  / %d  ", c, i, *lgSeq);
		if( i+1 <= *nStep){ minim=i+1; } 
		else{ minim=*nStep; }
		//printf("Sommet %d et Mini %d\n", i, minim);

		LogValue(currentLogClasse, proportions, 1.0, nbClasse);
		MiseAJour(sequence[i], currentLogClasse, logVarianceSur2, inverseVarianceSur2, moyennes, nbClasse);
		//printf("Tour : %d\n", i);
		Poids=-ComputeCost(currentLogClasse, nbClasse);
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
		
			MiseAJour(sequence[j], currentLogClasse, logVarianceSur2, inverseVarianceSur2, moyennes, nbClasse);
			Poids=-ComputeCost(currentLogClasse, nbClasse);

				
	
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
	free(inverseVarianceSur2 );
	free(logVarianceSur2);
	free(currentLogClasse);
//	printf("\n End of Dynamic Programming\n");

}
