#ifndef __LATTICE_SIEVE_H__
#define __LATTICE_SIEVE_H__

#include "mypair.h"
#include "HashTable.h"

void latticeRationalSieve(HashTable &pairs, const ulong *RB, const double *lRB,
	ulong iRB, ulong nRB, MyPair u, MyPair v, slong C, slong D, double logq, const fmpz_t m);
void latticeAlgebraicSieve(HashTable &abPairs, ulong &loc, slong num, HashTable &pairs, const fmpz_poly_t f,
	const MyPair *AB, const double *lAB, ulong iAB, ulong nAB, MyPair u, MyPair v, slong C, slong D);
void latticeSieve(const fmpz_poly_t f, const ulong *RB, const double *lRB, ulong nRB,
		   const MyPair *AB, const double *lAB, ulong nAB, MyPair *abPairs, ulong num, slong A, slong B, fmpz_t m);

#endif
