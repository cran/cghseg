#ifndef CGHSEG_LOGDENS_H
#define CGHSEG_LOGDENS_H

#include "numlib.h"

namespace cghseg {

numlib_matrix *
logdens(numlib_vector *xk, int P, numlib_vector *phi);

}
#endif
