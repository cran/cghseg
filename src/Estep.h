#ifndef CGHSEG_ESTEP_H
#define CGHSEG_ESTEP_H

#include "numlib.h"

namespace cghseg {

void Estep(numlib_matrix *logdensity, numlib_vector *phi,
	   numlib_matrix **tau, double *lvinc);

}
#endif

