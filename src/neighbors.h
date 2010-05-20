#ifndef CGHSEG_NEIGHBORS_H
#define CGHSEG_NEIGHBORS_H

#include "numlib.h"

#include "cghseg_types.h"

namespace cghseg {

void neighbors(numlib_vector *x, long Kmax,
	       numlib_vector *L, long k, 
	       param_struct param[], long param_size,
	       long P, long lmin, long lmax, bool vh,
	       numlib_vector **L_res, 
	       param_struct *param_res[],long *param_res_size);

}
#endif
