#include "GNFS.h"
/**
 *	Get the bound for smoothness (that upperbound for the rational factor base)
 */
ulong boundForSmoothness(ulong d, const fmpz_t n)
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
void prepareRationalBase(ulong *RB, double *lRB, ulong &nRB, ulong bound)
{
	ulong p;
	n_primes_t iter;
	n_primes_init(iter);
	for(nRB = 0; (p = n_primes_next(iter)) <= bound; nRB++)
	{
		RB[nRB] = p;
		lRB[nRB] = log(RB[nRB]);
	}
	//slong r = rand() % (nRB/10);
	//nRB -= r;
	assert(nRB <= MaxPrimeBufSize);
	n_primes_clear(iter);
}

void prepareAlgebraicBase(MyPair *AB, double *lAB, ulong &nAB, ulong size, fmpz_poly_t f)
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
	for(ulong i = 0; i < nAB; i++)
		lAB[i] = log(AB[i].p);
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


