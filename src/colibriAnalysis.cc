/* Copyright 2008 Guillem Rigaill <guillem.rigaill@curie.fr> 

   This file is part of colibri design for a fast segmentation in the mean

   Colibri is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   Colibri is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with Colibri; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "colibri.h"



/* argv 1 filein, 2 fileout, 3 fileoutint, 4 nb, 5 Kmax, 6 outPathFile, 7 fichierReal, 8 fileoutNb */

void colibri_Analysis (double *profil, int *nbi, int *Kmaxi, double *mini, double *maxi, char *fileOutInt)
{
	int nb=*nbi;
	int Kmax=*Kmaxi;
	double min=*mini;
	double max=*maxi;
	//double *minCostBefore = (double *) malloc(nb * sizeof(double));
	double *minCostBefore = new double[nb];
	//double *minCostCurrent = (double *) malloc(nb * sizeof(double));
	double *minCostCurrent = new double[nb];
	double *tmp;//1
	int minPosition;
	double minCurrent;
	//int * origine = (int *) malloc(nb * sizeof(int));
	int * origine = new int[nb];
    int * nombreCurrent = new int[nb];

	int i =0;
	double somme = 0;
	while(i < nb)
	{
		somme = somme + profil[i];
		minCostBefore[i] = - pow(somme, 2) / (i+1);
		origine[i]=0;
        nombreCurrent[i]=0;
		i++;
	}
	/* Save */
	//writeFileD(fileOut, minCostBefore[nb-1]);
	//writeFileI(fileOutInt, origine, nb, 0);
	
	/* Next */
	Polynome2 * p1;
	Liste * l1;  
	Polynome2 * pTest;
	int turn=1;
	char c = 13;
	//Polynome2 **stock= (Polynome2**)malloc(nb *sizeof(Polynome2*));//A
	Polynome2 **stock= new Polynome2* [nb];
	i=0;
	while(i < nb)
	{
		stock[i]=new Polynome2();
		i++;	
	}
	
	while( turn < Kmax)
	{
	  /* Print turn / Kmax */
	  //fprintf(stderr, "%c Turn :   %d  / %d  ", c, turn, Kmax);
	  /* initalisation */
	  i=turn;
	  stock[i]->reset(1.0, -2*profil[i], minCostBefore[turn -1],  turn);
	  stock[i]->setStatus(2);
	  l1 = new Liste(max, min, stock[i]);
	  /* Min */
	  l1->computeMinOrMax(&minCurrent, &minPosition);
	  minCostCurrent[i]=minCurrent;
	  origine[i] = i;
	  nombreCurrent[i]= l1->compte();
	  i=i+1;
	  while(i < nb)
		{
		 /* Slide 1 and Prune */
		 l1->computeRoots(minCostBefore[i-1]);
		 stock[i]->reset(0.0, 0.0, minCostBefore[i-1],  i);
		 l1->resetAllBorders(stock[i]);
		 l1->checkForDoublon();
		 l1->add(1.0, -2*profil[i], 0.0);

		 /* Compute Min */
		 l1->computeMinOrMax(&minCurrent, &minPosition);
		 minCostCurrent[i]=minCurrent;
         nombreCurrent[i]= l1->compte();
		 //origine[i] = minPosition;
		
		 /* iterate */
		 i++;	
	  	}

	  /* Save */
	  //writeFileD(fileOut, minCostCurrent[nb-1]);
	  writeFileI(fileOutInt, nombreCurrent, nb, turn);

	  /* */
	  tmp=minCostCurrent;
	  minCostCurrent=minCostBefore;
	  minCostBefore=tmp;
	
	
	  //delete(l1);
	  /* iterate */
	  turn++;	
	}
	
	/* Free All */
	/* free stock */
	i=0;
	while(i < nb)
	{
	    delete(stock[i]);	
		i++;
	}
	delete(stock);
	delete(minCostBefore);
	delete(minCostCurrent);
	delete(origine);
    delete(nombreCurrent);
	//std::cout << std::endl;

	/* Create matrix with Breakpoints positions for 0, ..., Kmax Breakpoints */
	

    //return 0;
}



















