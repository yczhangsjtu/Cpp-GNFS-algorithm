#ifndef __LATTICE_SIEVE_H__
#define __LATTICE_SIEVE_H__

#include "mypair.h"

double **createCDTable(slong C, slong D);
int **createMarkTable(slong C, slong D);
void freeCDTable(double **table, slong C);
void freeMarkTable(int **table, slong C);
void latticeRationalSieve(double **cdTable, const ulong *RB, const double *lRB,
	ulong iRB, ulong nRB, MyPair u, MyPair v, slong C, slong D, double logq, const fmpz_t m);
void latticeAlgebraicSieve(double **cdTable, ulong &loc, slong num, const fmpz_poly_t f,
	const MyPair *AB, const double *lAB, ulong iAB, ulong nAB, MyPair u, MyPair v, slong C, slong D);
void latticeSieve(const fmpz_poly_t f, const ulong *RB, const double *lRB, ulong nRB,
		   const MyPair *AB, const double *lAB, ulong nAB, MyPair *abPairs, ulong num, slong A, slong B, fmpz_t m);

#endif
