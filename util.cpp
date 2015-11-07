#include "GNFS.h"
/**
 *	Product of a set of numbers modular n
 */
void productMod(fmpz_t res, ulong *a, ulong k, const fmpz_t n)
{
	fmpz_one(res);
	for(ulong i = 0; i < k; i++)
	{
		fmpz_mul_ui(res,res,a[i]);
		fmpz_mod(res,res,n);
	}
}

/**
 *	Clear an array of fmpz_t.
 */
void freeArrayOfFmpz(fmpz_t *array, ulong size)
{
	for(ulong i = 0; i < size; i++)
		fmpz_clear(array[i]);
}

/**
 *	Get the sum of an array of doubles.
 *
 *	This is tricky and recursive method is used instead of adding them one by
 *	one because of the properties of floating point computation in computer:
 *	the precision is badly hurt when summing two floating point numbers significantly
 *	different in order.
 */
double doublesum(double array[], ulong l, ulong r)
{
	if(r-l == 0) return 0;
	if(r-l == 1) return array[l];
	return doublesum(array,l,(l+r)/2) + doublesum(array,(l+r)/2,r);
}

/**
 *	Some print procedures.
 */
void printListOfNumbers(FILE *file, slong *A, ulong s, ulong N)
{
	for(ulong i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) fprintf(file,"\n");
		fprintf(file,"%10ld",A[i]);
	}
	fprintf(file,"\n");
}

void printListOfNumbers(FILE *file, ulong *A, ulong s, ulong N)
{
	for(ulong i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) fprintf(file,"\n");
		fprintf(file,"%10ld",A[i]);
	}
	fprintf(file,"\n");
}

void printListOfPairs(FILE *file, MyPair *A, ulong s, ulong N)
{
	for(ulong i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) fprintf(file,"\n");
		fprintf(file,"%10ld%10ld",A[i].r,A[i].p);
	}
	fprintf(file,"\n");
}

/**
 *	Select part of an array corresponding to the 1's in a 0-1 vector,
 *	discard the rest and keep only the selected ones, update the size of array
 */
void select(MyPair *pairs, int *vec, ulong I, ulong &n)
{
	ulong loc = 0;
	for(ulong i = 0; i < I; i++)
	{
		if(vec[i]) pairs[loc++] = pairs[i];
	}
	n = loc;
	assert(n > 0);
}
