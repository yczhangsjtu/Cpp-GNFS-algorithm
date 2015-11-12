#ifndef __ALGEBRAIC_H__
#define __ALGEBRAIC_H__

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include "mypair.h"

void algebraicSieve(double *nf_sieve_array, const fmpz_poly_t f, const MyPair *AB,
	const double *lAB, ulong nAB, ulong sieve_len, slong A, ulong b);

#endif
