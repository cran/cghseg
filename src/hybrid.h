#ifndef CGHSEG_HYBRID_H
#define CGHSEG_HYBRID_H

#include "numlib.h"

#include "cghseg_types.h"

namespace cghseg {

  void hybrid(numlib_vector *x, long P, long Kmax, long lmin, long lmax,
	      bool vh, bool fast, numlib_matrix **Linc_res, param_struct *param_res[],
	    long *param_res_size);
}
#endif
