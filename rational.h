#ifndef __RATIONAL_H__
#define __RATIONAL_H__

#include <flint/fmpz.h>
#include "mypair.h"

void rationalSieve(double *sieve_array, ulong sieve_len, const ulong *RB,
	const double *lRB, ulong nRB, slong A, const fmpz_t bm);

#endif
