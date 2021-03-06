#ifndef __UTIL_H__
#define __UTIL_H__

#include <flint/fmpz.h>
#include "mypair.h"

void productMod(fmpz_t res, const ulong *a, ulong k, const fmpz_t n);
void freeArrayOfFmpz(fmpz_t *array, ulong size);
double doublesum(const double array[], ulong l, ulong r);
void printListOfNumbers(FILE *f, slong *A, ulong s, ulong N);
void printListOfNumbers(FILE *f, ulong *A, ulong s, ulong N);
void printListOfPairs(FILE *f, MyPair *A, ulong s, ulong N);
void select(MyPair *pairs, int *vec, ulong I, ulong &n);

#endif
