#ifndef __SIEVE_H__
#define __SIEVE_H__

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <cmath>
#include "mypair.h"
#include "util.h"

void sieve(fmpz_poly_t f, const ulong *RB, const double *lRB, ulong nRB,
		   const MyPair *AB, const double *lAB, ulong nAB, MyPair *abPairs,
		   ulong num, slong N, fmpz_t m,
		   int mynode, int totalnodes, MPI_Status *status);

#endif
