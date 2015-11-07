#ifndef __FACTOR_BASE_H__
#define __FACTOR_BASE_H__

#include "mypair.h"

ulong boundForSmoothness(ulong d, const fmpz_t n);
void prepareRationalBase(ulong *RB, double *lRB, ulong &nRB, ulong bound);
void prepareAlgebraicBase(MyPair *AB, double *lAB, ulong &nAB, ulong size, fmpz_poly_t f);
void prepareQuadraticBase(MyPair *QB, ulong &nQB, ulong min, ulong max, fmpz_poly_t &f);

#endif
