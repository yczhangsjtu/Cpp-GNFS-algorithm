#ifndef __POLY_SELECT_H__
#define __POLY_SELECT_H__

#include <iostream>
#include <cmath>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

bool selectPolynomial(const fmpz_t n, fmpz_poly_t f, fmpz_t m, ulong &d);

#endif
