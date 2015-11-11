#ifndef __FACTOR_BASE_H__
#define __FACTOR_BASE_H__

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include "mypair.h"

ulong boundForSmoothness(slong d, const fmpz_t n);
void prepareRationalBase(ulong *RB, ulong &nRB, ulong bound);
void prepareAlgebraicBase(MyPair *AB, ulong &nAB, ulong size, fmpz_poly_t f);
void prepareQuadraticBase(MyPair *QB, ulong &nQB, ulong min, ulong max, fmpz_poly_t &f);

#endif
