#ifndef __ALGEBRAIC_H__
#define __ALGEBRAIC_H__

#include "mypair.h"

void algebraicSieve(double *nf_sieve_array, const fmpz_poly_t f, const MyPair *AB, const double *lAB, ulong nAB,
					ulong sieve_len, slong A, ulong b);

#endif
