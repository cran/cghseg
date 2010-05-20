#ifndef CGHSEG_MSTEP_H
#define CGHSEG_MSTEP_H

#include "numlib.h"

namespace cghseg {

numlib_vector * Mstep(numlib_vector *x, numlib_matrix *rupt, numlib_matrix *tau,
		   numlib_vector *phi, bool vh);

}

#endif
