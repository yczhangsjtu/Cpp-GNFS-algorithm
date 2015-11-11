#ifndef __LATTICE_UTIL_H__
#define __LATTICE_UTIL_H__

#include <flint/fmpz.h>
#include "mypair.h"

void gaussianLatticeReduce(MyPair &u, MyPair &v);
void getBound(slong &E1, slong &E2, slong &F1, slong &F2, slong C, slong D, MyPair s, MyPair t);

#endif
