#ifndef __LINEAR_H__
#define __LINEAR_H__

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod_mat.h>
#include "mypair.h"

void formMatrix(nmod_mat_t mat, ulong I, ulong J, const fmpz_t m, const fmpz_poly_t f, const MyPair *abPairs,
				const ulong *RB, ulong nRB, const MyPair* AB, ulong nAB, const MyPair* QB, ulong nQB);
void solveMatrix(nmod_mat_t mat, ulong I, ulong J, int *vec);

#endif
