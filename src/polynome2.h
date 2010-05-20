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

#include <iostream>
#include <math.h>
#include <float.h>
#include <stdlib.h>


class Polynome2 {
 private:
     double a2, a1, a0;
/* racine du polynome - A*/
     double rac1, rac2;
/* status of the polynome 0 = roots not computed, 1 = computed*/
     int status;
     int origine;
 public:
     /* constructors and destructors */
     Polynome2();
     ~Polynome2();
     Polynome2(double, double, double, int);
     /* a few operations */
	/* reset */
     void reset(double, double, double, int);

     /* getter and setter */
     double geta2();
     double geta1();
     double geta0();
     double getRacine1();
     double getRacine2();

     void seta2(double);
     void seta1(double);
     void seta0(double);
     void setRacine1(double);
     void setRacine2(double);
     void setStatus(int);
     int getStatus();
     int getOrigine();

     /* Delta and others */
     double eval(double);
     double delta();
	 /* Delta  of the Polynome - double */
     double delta(double);

     void roots();
	  /* Roots  of the Polynome - double */
     void roots(double);
	 void add(double, double, double);
     void minOrMax(double *, double *, int *);
     /* print and others */
     void show();
};

