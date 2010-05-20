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

#include "polynome2.h"
#include "numlib.h"


class Liste {
 private:
     double max, min;
	 Polynome2 *poly;
	 Liste *next;
 public:
    /* constructors and destructors */
    Liste();
	Liste(double max_, double min_);
	Liste(double max_, double min_, Polynome2 *poly_);
	Liste(Polynome2 *poly_);
	~Liste();
	/* fonction setter and getter */
    double getMax();
    void setMax(double max_);
	double getMin();
    void setMin(double min_);

	void setPolynome(Polynome2 * poly_);
	Polynome2 *getPolynome();

	Liste * getNext();
	void setNext(Liste * next_);

	/* Useful */
	void setToNull();
	void insert(Liste * maillon_);
    int compte();
	Liste * removeDoublon();
	void checkForDoublon();

	/* show and others */
	void show();
	void showAllNext();
	
	/* */
	void computeRoots(double);
	void add(double, double, double);
	void computeMinOrMax(double*, int*);
	void resetMaillonBorders(Polynome2*);
	void resetAllBorders(Polynome2*);
	
};

