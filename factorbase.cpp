#include "GNFS.h"
#include "poly.h"
#include "factorbase.h"
#include "util.h"

int MaxPrime = DefaultMaxPrime;
/**
 *	Get the bound for smoothness (that upperbound for the rational factor base)
 */
ulong boundForSmoothness(slong d, const fmpz_t n)
{
#ifdef DEBUG
	assert(d > 0);
	assert(n > 0);
#endif
	double dlogd = d*log(d);
	double temp = 1.0/d * fmpz_dlog(n);
#ifdef DEBUG
	assert(temp > 0);
#endif
	double e = dlogd + sqrt(dlogd*dlogd + 4*temp*log(temp));
	return 0.5*exp(0.5*e);
}

/**
 *	Prepare the bases
 *
 *	1. rational base: a set of prime numbers
 *	2. algebraic base: a set of pairs (p,r), where p is prime number,
 *		r is a root of f modular p.
 *	3. quadratic base: a set of pairs (q,s), where q is a prime number,
 *		s is a root of f modular p. It is required that (q,s) should never
 *		divide the pairs (a,b) sieved out by the (p,r) pairs, so the smallest
 *		q is usually larger than the largest p in algebraic base.
 */
void prepareRationalBase(ulong *RB, ulong &nRB, ulong bound)
{
	ulong p;
	n_primes_t iter;
	n_primes_init(iter);
	for(nRB = 0; (p = n_primes_next(iter)) <= bound; nRB++)
	{
		RB[nRB] = p;
	}
	//slong r = rand() % (nRB/10);
	//nRB -= r;
	assert(nRB <= MaxPrimeBufSize);
	n_primes_clear(iter);
}

void prepareAlgebraicBase(MyPair *AB, ulong &nAB, ulong size, fmpz_poly_t f)
{
	n_primes_t iter;
	n_primes_init(iter);
	for(ulong i = 0; i < size; i++)
	{
		ulong p = n_primes_next(iter);
		ulong roots[16]; ulong nroot;
		rootsMod(f,p,roots,nroot);
		for(ulong j = 0; j < nroot; j++)
		{
			AB[nAB] = MyPair(roots[j],p);
			nAB++;
		}
	}
	n_primes_clear(iter);
	assert(nAB <= MaxPrimeBufSize);
}

void prepareQuadraticBase(MyPair *QB, ulong &nQB, ulong min, ulong max, fmpz_poly_t &f)
{
	ulong p;
	n_primes_t iter;
	n_primes_init(iter);
	n_primes_jump_after(iter,min);
	for(ulong i = 0; (p=n_primes_next(iter)) <= max; i++)
	{
		ulong roots[16]; ulong nroot;
		rootsMod(f,p,roots,nroot);
		for(ulong j = 0; j < nroot; j++)
		{
			QB[nQB] = MyPair(roots[j],p);
			nQB++;
		}
	}
	assert(nQB <= MaxPrimeBufSize);
	n_primes_clear(iter);
}

using namespace std;
int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		cerr << "Usage: factorbase inputfile outputfile" << endl;
		exit(-1);
	}
	FILE *input = fopen(argv[1],"r");
	if(!input) perror(argv[1]);
	FILE *output = fopen(argv[2],"w");
	if(!output) perror(argv[2]);

	fmpz_t n, m;
	fmpz_poly_t f;

	fmpz_init(n);
	fmpz_init(m);
	fmpz_poly_init(f);

	fmpz_fread(input,n);
	fmpz_fread(input,m);
	fmpz_poly_fread(input,f);
	slong d = fmpz_poly_degree(f);

	ulong smoothBound = boundForSmoothness(d,n);
	ulong RB[MaxPrimeBufSize], nRB = 0, nAB = 0, nQB = 0;
	MyPair AB[MaxPrimeBufSize], QB[MaxPrimeBufSize];
	prepareRationalBase(RB,nRB,smoothBound);
	prepareAlgebraicBase(AB,nAB,nRB,f);
	prepareQuadraticBase(QB,nQB,smoothBound,smoothBound+20*log(smoothBound),f);
	fmpz_fprint(output,n); fprintf(output,"\n");
	fmpz_fprint(output,m); fprintf(output,"\n");
	fprintf(output,"%lu\n",smoothBound);
	fmpz_poly_fprint(output,f); fprintf(output,"\n");
	fprintf(output,"%lu\n",nRB);
	printListOfNumbers(output,RB,nRB,10);
	fprintf(output,"%lu\n",nAB);
	printListOfPairs(output,AB,nAB,5);
	fprintf(output,"%lu\n",nQB);
	printListOfPairs(output,QB,nQB,5);

	fmpz_clear(n);
	fmpz_clear(m);
	fmpz_poly_clear(f);

	fclose(input);
	fclose(output);
	return 0;
}
