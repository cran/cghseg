#include "ApproxSeg.h"

void ClassiSeg(double *sequence, int *lgSeq, int *nStep, double *res1, int *res2, int *nbClasse, double *moyennes){
  /* Compteurs et autres */
  int i, j, k, l, indice;
  /* Variable temporaires */
  double * vTmp;
  vTmp = (double *) malloc( *nbClasse * sizeof(double));
  for(i =0; i < *nbClasse; i++) vTmp[i]=0;

  int * whichCome;
  whichCome = (int *) malloc( *nbClasse * sizeof(int));
  for(i =0; i < *nbClasse; i++) whichCome[i]=-1;

  double min, tMin;
  int whichMin;

  for(i = 0; i < *lgSeq; i++)
   for(k= 0; k < *nStep; k++) 
     res1[(*lgSeq)*k+i] = A_POSINF;

  /* intialisation */
  for(i =0; i < *lgSeq; i++)
  {
    min=A_POSINF;
    for(l =0; l < *nbClasse; l++)
    {		
      vTmp[l] = vTmp[l] + (sequence[i] - moyennes[l])*(sequence[i] - moyennes[l]);	   
      if(min > vTmp[l]) min  = vTmp[l];	
    }
    res1[i] = min;
    res2[i] = 0;  	
  }

  for(k =1; k < *nStep; k++)
  {
    for(l =0; l < *nbClasse; l++) vTmp[l]= A_POSINF;
    for(l =0; l < *nbClasse; l++) whichCome[l]=0;
	
    for(i=k; i < *lgSeq; i++)
    {
      min=A_POSINF;
      for(l =0; l < *nbClasse; l++)
      {
        indice = (*lgSeq)*(k-1)+i-1;
        if( vTmp[l] > res1[indice] ) /* on change de segment */
        {
	  vTmp[l] = res1[indice] + (sequence[i] - moyennes[l])*(sequence[i] - moyennes[l]);
          whichCome[l]=i;	
	} else /* pas de changement de segment */
        {
         vTmp[l] = vTmp[l] + (sequence[i] - moyennes[l])*(sequence[i] - moyennes[l]);
        }

	if(vTmp[l] < min)
        {
	  min=vTmp[l];
	  whichMin=l;
	}
      }	
      indice = (*lgSeq)*k+i;			
      res1[indice] = min;
      res2[indice] = whichCome[whichMin];
			
     }
   }
	
   free(vTmp);
   free(whichCome);
}

