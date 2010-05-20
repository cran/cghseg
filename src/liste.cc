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

#include "liste.h"


Liste::Liste()
{
  this->max=0.0;
  this->min=0.0;
  this->next=NULL;
}

Liste::Liste(double max_, double min_)
{
  this->max=max_;
  this->min=min_;
  this->next=NULL;
}

Liste::Liste(double max_, double min_, Polynome2 *poly_)
{

  this->max=max_;
  this->min=min_;
  this->next=NULL;
  this->poly= poly_;

}

Liste::Liste( Polynome2 *poly_)
{

  this->max=0.0;
  this->min=0.0;
  this->next=NULL;
  this->poly= poly_;

}
Liste::~Liste()
{
	 delete next;
	 delete poly;
}
/* Setter and Getter */
/* */  
double Liste::getMax()
{
	return(this->max);
}

void Liste::setMax(double max_)
{
	this->max = max_;
}

double Liste::getMin()
{
	return(this->min);
}

void Liste::setMin(double min_)
{
	this->min = min_;
}

void Liste::setPolynome(Polynome2 * poly_)
{
		this->poly=poly_;	
}
Polynome2* Liste::getPolynome()
{
		return(this->poly	);
}
/* */  
Liste * Liste::getNext()
{
	return(this->next);
}
void Liste::setNext(Liste * next_)
{
	this->next = next_;
}


/* Usefull */

void Liste::setToNull()
{
	max=NULL;
	min=NULL;
	poly=NULL;
	next=NULL;
}
void Liste::insert(Liste * maillon_)
{
	maillon_->setNext(this->getNext());
	this->setNext(maillon_);
}

int Liste::compte(){
	Liste *l;
	int tmp = 0;
	l=this;
	while(l != NULL){
		tmp= tmp+1;
		l=l->getNext();	
	}
	return(tmp);
}

Liste * Liste::removeDoublon()
{
	Liste *next = this->getNext();
	if(next != NULL)
	{
		if(next->getPolynome() == this->getPolynome())
		{
			this->setMin(next->getMin());
			this->setNext(next->getNext());
			next->setToNull();
			delete next;
			return(this);
		} else 
		{ 
			return(next);
		}
	} else 
	{
		return(NULL);
	}
}


void Liste::checkForDoublon()
{
	Liste *l = this;
	while(l != NULL)
	{
		l=l->removeDoublon();
	}
}
/* Show and others */  
void Liste::show()
{
	std::cout << "Max : " << this->getMax() << ", Min : " << this->getMin()<< std::endl;
	this->poly->show();
}

void Liste::showAllNext()
{
	Liste *l;
	l = this;
	while(l != NULL)
	{
		l->show();
		l=l->getNext();	
	}
	
}

void Liste::computeRoots(double a0_)
{
	Liste *l;
	l=this;
	while(l != NULL)
	{
		l->getPolynome()->roots(a0_);
		l=l->getNext();	
	}
}

void Liste::add(double a2_, double a1_, double a0_)
{
	Liste *l;
	l=this;
	while(l != NULL)
	{
		l->getPolynome()->add(a2_, a1_, a0_);
		l=l->getNext();	
	}
}

void Liste::computeMinOrMax(double * min, int * which)
{
	Liste *l;
	double tmp = NUMLIB_POSINF;
	*min = NUMLIB_POSINF;
	* which=-1;
	l=this;
	while(l != NULL)
	{
		l->getPolynome()->minOrMax(min, &tmp, which);
		l=l->getNext();	
	}
}

void Liste::resetMaillonBorders(Polynome2 *poly_)
{
	//if(this->getPolynome()->getRacine2() == NAN)
	if(this->getPolynome()->getRacine2() == NULL)
	//if( isnan(this->getPolynome()->getRacine2()) )
	{
		this->setPolynome(poly_);
	} else if(this->getPolynome()->getRacine1() >= this->getMax())
	{
		if(this->getPolynome()->getRacine2() >= this->getMax())
		{
			this->setPolynome(poly_);
		} else if(this->getPolynome()->getRacine2() > this->getMin()) 
			{
				Liste * maillon= new Liste(this->getPolynome()->getRacine2(), this->getMin(), poly_);
				this->insert(maillon);
				this->setMin(this->getPolynome()->getRacine2());
			} else 
			{
			}
		
	
	} else if( this->getPolynome()->getRacine1() > this->getMin())
		{
		
			if(this->getPolynome()->getRacine2() > this->getMin())
			{
				Liste *maillon3 = new Liste(this->getPolynome()->getRacine2(), this->getMin(), poly_);
				Liste *maillon2 = new Liste(this->getPolynome()->getRacine1(), this->getPolynome()->getRacine2(), this->getPolynome());
				this->setMin(this->getPolynome()->getRacine1());
				this->setPolynome(poly_);
				this->insert(maillon3);
				this->insert(maillon2);
			} else
			{
				Liste *maillon2 = new Liste(this->getPolynome()->getRacine1(), this->getMin(), this->getPolynome());
				this->setMin(this->getPolynome()->getRacine1());
				this->setPolynome(poly_);
				this->insert(maillon2);
			}
		
		} else 
		{
			this->setPolynome(poly_);
		}
}

void Liste::resetAllBorders(Polynome2 *poly_)
{

	Liste *lCurrent, *lNext;
	lCurrent = this;
	lNext = this->getNext();
	lCurrent->resetMaillonBorders(poly_);
	lCurrent=lNext;
	while(lCurrent != NULL)
	{
		lNext=lCurrent->getNext();	
		lCurrent->resetMaillonBorders(poly_);
		lCurrent=lNext;
	}
}

