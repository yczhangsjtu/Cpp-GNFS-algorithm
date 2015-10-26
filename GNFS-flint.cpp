#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <list>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <limits.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>

//#define DEBUG 1
#define PRINT_PROCESS 1
#define ONLY_RESULT 0
#define PRINT_SIEVE_PROCESS 0

#if(PRINT_PROCESS)
#define PRINT_MDF
//#define PRINT_SMOOTH_BOUND
//#define PRINT_RATIONAL_BASE
//#define PRINT_ALGEBRAIC_BASE
//#define PRINT_QUADRATIC_BASE
//#define PRINT_SELECTED_ABPAIRS
//#define PRINT_MATRIX
//#define PRINT_SELECTED_SQUARE_ABPAIRS
//#define PRINT_PROD_ABM
//#define PRINT_UPPDER_BOUND
//#define PRINT_PRIMES
//#define PRINT_XI
#endif

using namespace std;

/*******************************************************************************
 *	Global Constants
 ******************************************************************************/
const int MaxPrimeBufSize = 30000;
const int MaxPrime = 20000000;
const int MaxSelectedPrimes = 10000;
const int MaxB = 10240;

/*******************************************************************************
 *	Structure definitions
 ******************************************************************************/

/**
 *	Pair of integer
 */
typedef struct MyPair
{
	slong r;
	slong p;
	MyPair(slong a, slong b){r=a; p=b;}
	MyPair(){r=0;p=0;}
} MyPair;

/*******************************************************************************
 * Tools used in the GNFS
 ******************************************************************************/

/* Tools about polynomials and norms. *****************************************/

/**
 *	Get the 'norm' of a polynomial f at pair (a,b)
 *	The 'norm' of a polynomial f is defined by
 *	Norm[a,b](f) = (-b)^d f(-a/b)
 *		= a^d - c1 a^(d-1) b + c2 a^(d-2) b^2 - ... + (-1)^(d-1) cd-1 a b^(d-1)
 *			+ (-1)^d cd b^d
 *	where d is the degree of f, and ci are coefficients of f.
 *	f(x) = x^d + c1 x^(d-1) + ... + cd-1 x + cd
 */
void norm(fmpz_t nm, const fmpz_poly_t f, const fmpz_t a, const fmpz_t b)
{
	fmpz_t poa,pob,mb,c,ab;
	fmpz_init(poa); /*Power of a*/
	fmpz_init(pob); /*Power of -b*/
	fmpz_init(mb); /*-b*/
	fmpz_init(c);
	fmpz_init(ab);
	fmpz_one(poa);
	fmpz_one(pob);
	fmpz_zero(nm);
	fmpz_neg(mb,b);
	ulong d = fmpz_poly_degree(f);

	/* If a = 0, then the norm is simply (-b)^d */
	if(fmpz_is_zero(a))
	{
		for(ulong i = 0; i < d; i++)
			fmpz_mul(pob,pob,mb);
		fmpz_poly_get_coeff_fmpz(c,f,0);
		fmpz_mul(nm,pob,c);
	}
	else
	{
		/* First raise a to power d. */
		for(ulong i = 0; i < d; i++)
			fmpz_mul(poa,poa,a);
		/* In each step multiply power of a and power of -b, add to norm
		 * then multiply pob by -b, divide poa by a */
		for(slong i = d; i >= 0; i--)
		{
			fmpz_poly_get_coeff_fmpz(c,f,i);
			fmpz_mul(ab,poa,pob);
			fmpz_mul(ab,ab,c);
			fmpz_add(nm,nm,ab);
			fmpz_mul(pob,pob,mb);
			fmpz_fdiv_q(poa,poa,a);
		}
	}

	fmpz_clear(poa);
	fmpz_clear(pob);
	fmpz_clear(mb);
	fmpz_clear(c);
	fmpz_clear(ab);
}

/**
 *	Given a number num, and a set of primes in ps[], check if num can be factored
 *	over this set of primes.
 */
bool isSmooth(const fmpz_t num, const ulong *ps, ulong np)
{
	fmpz_t k;
	fmpz_init_set(k,num);
	for(ulong i = 0; i < np; i++)
	{
		ulong p = ps[i];
		while(!fmpz_is_zero(k) && fmpz_divisible_si(k,p))
			fmpz_fdiv_q_ui(k,k,p);
	}
	fmpz_clear(k);
	return fmpz_is_pm1(k);
}

bool isSmooth(const fmpz_t num, const MyPair *ps, ulong np)
{
	fmpz_t k;
	fmpz_init_set(k,num);
	for(ulong i = 0; i < np; i++)
	{
		ulong p = ps[i].p;
		while(!fmpz_is_zero(k) && fmpz_divisible_si(k,p))
			fmpz_fdiv_q_ui(k,k,p);
	}
	fmpz_clear(k);
	return fmpz_is_pm1(k);
}

/**
 *	Find the root of a polynomial modular prime p.
 *
 *	The flint library does not implement this function directly,
 *	but polynomial factorization mod p is provided.
 *	The Algorithm here is simply factoring f mod p first and collect
 *	the unaries.
 */
void rootsMod(const fmpz_poly_t f, ulong p, ulong *roots, ulong &nroot)
{
#ifdef DEBUG
	assert(p > 1);
#endif
	ulong d = fmpz_poly_degree(f);
	fmpz_t c;
	nmod_poly_factor_t fac;
	nmod_poly_t g;
	fmpz_init(c);
	nmod_poly_factor_init(fac);
	nmod_poly_factor_fit_length(fac,d);
	nmod_poly_init(g,p);
	for(int i = 0; i <= d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,f,i);
		nmod_poly_set_coeff_ui(g,i,fmpz_mod_ui(c,c,p));
	}
	nmod_poly_factor(fac,g);

	ulong n = fac->num;
	nroot = 0;
	for(long i = 0; i < n; i++)
	{
		nmod_poly_struct* h = fac->p+i;
		if(nmod_poly_degree(h) == 1)
		{
			ulong a = nmod_poly_get_coeff_ui(h,0);
			roots[nroot] = (p-a)%p;
			nroot++;
		}
	}
	fmpz_clear(c);
	nmod_poly_clear(g);
	nmod_poly_factor_clear(fac);
}

/**
 *	Legender symbol (in fact jacobi symbol is the generalization of it)
 *
 *	Jacobi symbol is implemented by flint library, however, it is required
 *	that a < p, so we have to pack it in our own function to call it with
 *	a >= p
 */
int Leg(slong a, ulong p)
{
	fmpz_t c,q;
	fmpz_init_set_ui(q,p);
	fmpz_init_set_si(c,a);
	fmpz_mod(c,c,q);
	int l = fmpz_jacobi(c,q);
	fmpz_clear(c);
	fmpz_clear(q);
	return l;
}

/**
 *	Return true if fx is irreducible mod p.
 */
bool irreducibleMod(const fmpz_poly_t fx, ulong p)
{
	ulong d = fmpz_poly_degree(fx);
	fmpz_t c;
	nmod_poly_t f;
	fmpz_init(c);
	nmod_poly_init(f,p);
	for(slong i = 0; i <= d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,fx,i);
		nmod_poly_set_coeff_ui(f,i,fmpz_mod_ui(c,c,p));
	}
	int r = nmod_poly_is_irreducible(f);
	fmpz_clear(c);
	nmod_poly_clear(f);
	return r;
}

/**
 *	Select a non-quadratic-residual in the field F_p[X]/<fx>, which is easy because there are
 *	half of them in the field.
 *
 *	When fx is irreducible, an elment g in F_p[x]/<fx> is nonresidual iff g to power (q-1)/2==-1
 *	Where q = p^d, and d is the degree of fx.
 */
void selectNonResidual(nmod_poly_t px, const nmod_poly_t fx, ulong p, const fmpz_t e, ulong d)
{
	mpz_t ee;
	mpz_init(ee);
	fmpz_get_mpz(ee,e);
	flint_rand_t frt;
	flint_randinit(frt);
	nmod_poly_t mp;
	nmod_poly_init(mp,p);
	while(true)
	{
		nmod_poly_randtest_monic(px,frt,d);
		nmod_poly_powmod_mpz_binexp(mp,px,ee,fx);
#ifdef DEBUG
		assert(nmod_poly_degree(mp)==0);
#endif
		if(nmod_poly_get_coeff_ui(mp,0) == p-1) break;
	}
	flint_randclear(frt);
	nmod_poly_clear(mp);
	mpz_clear(ee);
}

/**
 *	It is assume that the order of px is 2^r, this function returns the r.
 *
 *	That is, px^2^r mod fx = 1
 */
ulong computeOrder(const nmod_poly_t px, const nmod_poly_t fx)
{
	mp_limb_t p = nmod_poly_modulus(px);
	nmod_poly_t g;
	nmod_poly_init(g,p);
	nmod_poly_set(g,px);
	ulong i;
	for(i = 0; i < 10000; i++)
	{
		if(nmod_poly_is_one(g)) break;
		nmod_poly_powmod_ui_binexp(g,g,2,fx);
	}
	nmod_poly_clear(g);
	if(i >= 10000) return -1;
	return i;
}

/* Something else**************************************************************/

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
void printListOfNumbers(slong *A, ulong s, ulong N)
{
	for(ulong i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) cout << endl;
		cout << setw(10) << A[i];
	}
	cout << endl;
}

void printListOfNumbers(ulong *A, ulong s, ulong N)
{
	for(ulong i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) cout << endl;
		cout << setw(10) << A[i];
	}
	cout << endl;
}

void printListOfPairs(MyPair *A, ulong s, ulong N)
{
	for(ulong i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) cout << endl;
		cout << "(" << A[i].r << "," << A[i].p << ") ";
	}
	cout << endl;
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

/*******************************************************************************
 * Procedures of GNFS
 ******************************************************************************/

/**
 * Given integer n, return f, m, d, such that f is a degree d polynomial,
 * and f(m) = n.
 *
 * f: Integer array representing the coefficients of polynomial
 * m: the number selected
 * d: degree of the polynomial
 */
void selectPolynomial(const fmpz_t n, fmpz_poly_t f, fmpz_t m, ulong &d)
{
	fmpz_t lead,tail,N,Nmm;
#ifdef DEBUG
	assert(fmpz_cmp_ui(n,0) > 0);
	assert(fmpz_dlog(n) > 0);
#endif
	fmpz_poly_init(f);
	fmpz_init(m);
	fmpz_init(N);
	fmpz_init(lead);
	fmpz_init(tail);
	fmpz_init(Nmm);
	/* For moderate size of n, d is usually 3. When n is very large, d could be 5. */
	d = pow(3*fmpz_dlog(n)/log(fmpz_dlog(n)),0.333);
	/* It is required that d is odd. */
	if(d % 2 == 0) d++;
	/* Loop until proper f,m are found. */
	while(true)
	{
		do /* Select m. Set m to root of d first, then subtract a random number from it. */
		{
			fmpz_root(m,n,d); /* Set m to the d'th root of n rounded to floor. */
			if(fmpz_cmp_ui(m,2000)>0) fmpz_sub_ui(m,m,rand()%100);
			else if(fmpz_cmp_ui(m,100)>0) fmpz_sub_ui(m,m,rand()%10);
		} while(fmpz_divisible(n,m));

		fmpz_set(N,n); /* Make a copy of n, then compute the coefficients of f by repeatedly
		deviding N by m and take the remainders. */
		for(int i = 0; i <= d; i++)
		{
			fmpz_mod(Nmm,N,m);
			fmpz_poly_set_coeff_fmpz(f,i,Nmm);
			fmpz_fdiv_q(N,N,m);
		}
		fmpz_poly_get_coeff_fmpz(lead,f,d);
		fmpz_poly_get_coeff_fmpz(tail,f,0);
		/* f is required to be monic, nonzero constant term, and the constant term is not 1. */
		if(fmpz_is_one(lead) && !fmpz_is_zero(tail) && !fmpz_is_one(tail)) break;
	}
	fmpz_clear(lead);
	fmpz_clear(tail);
	fmpz_clear(N);
	fmpz_clear(Nmm);
}

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


/*
 * sieve_array: The array to accumulate
 * RB: Factor base
 * nRB: Size of factor base
 * sieve_len: sieve length
 */
void rationalSieve(double *sieve_array, ulong sieve_len, const ulong *RB, const double *lRB, ulong nRB, slong A, const fmpz_t bm)
{
	fmpz_t Abm,Abmi,fA,Abmmp;
	fmpz_init(Abm);
	fmpz_init(Abmi);
	fmpz_init(fA);
	fmpz_init(Abmmp);
	fmpz_set_si(fA,A);
	fmpz_add(Abm,bm,fA);
	/*Zerolize the sieve array*/
	for(ulong i = 0; i < sieve_len; i++)
	{
		fmpz_add_ui(Abmi,Abm,i);
		fmpz_abs(Abmi,Abmi);

		if(fmpz_is_zero(Abmi)) continue;
		sieve_array[i] = -fmpz_dlog(Abmi);
	}
	
	for(ulong i = 0; i < nRB; i++)
	{
		ulong p = RB[i];
		ulong f = fmpz_mod_ui(Abmmp,Abm,p);
		ulong loc = (f == 0? 0: p - f);
		
		while(loc < sieve_len)
		{
			fmpz_add_ui(Abmi,Abm,loc);
			sieve_array[loc] += lRB[i];
			loc += p;
		}
	}
	fmpz_clear(Abm);
	fmpz_clear(Abmi);
	fmpz_clear(fA);
	fmpz_clear(Abmmp);
}

void algebraicSieve(double *nf_sieve_array, const fmpz_poly_t f, const MyPair *AB, const double *lAB, ulong nAB,
					ulong sieve_len, slong A, ulong b)
{
	fmpz_t normab,fA,fb,nm,Ai,r,Abr,Abrmp;
	fmpz_init(normab);
	fmpz_init(fA);
	fmpz_init(fb);
	fmpz_init(nm);
	fmpz_init(Ai);
	fmpz_init(r);
	fmpz_init(Abr);
	fmpz_init(Abrmp);
	fmpz_set_si(fA,A);
	fmpz_set_si(fb,b);
	/*Zerolize the sieve array*/
	for(ulong i = 0; i < sieve_len; i++)
	{
		fmpz_add_ui(Ai,fA,i);
		norm(nm,f,Ai,fb);
		fmpz_abs(nm,nm);
		if(fmpz_is_zero(nm)) continue;
		nf_sieve_array[i] = -fmpz_dlog(nm);
	}
	
	for(ulong i = 0; i < nAB; i++)
	{
		ulong p = AB[i].p;
		fmpz_set_si(r,AB[i].r);
		fmpz_set(Abr,fA);
		fmpz_addmul(Abr,fb,r);
		ulong f = fmpz_mod_ui(Abrmp,Abr,p);
		ulong loc = (f == 0? 0: p - f);
		
		while(loc < sieve_len)
		{
			nf_sieve_array[loc] += lAB[i];
			loc += p;
		}
	}
	fmpz_clear(normab);
	fmpz_clear(fA);
	fmpz_clear(fb);
	fmpz_clear(nm);
	fmpz_clear(Ai);
	fmpz_clear(r);
	fmpz_clear(Abr);
	fmpz_clear(Abrmp);
}

void sieve(fmpz_poly_t f, const ulong *RB, const double *lRB, ulong nRB,
		   const MyPair *AB, const double *lAB, ulong nAB, MyPair *abPairs, ulong num, slong N, fmpz_t m)
{
	fmpz_t bm,fa,fb,nm,abm,gcd;
	fmpz_init(bm);
	fmpz_init(fa);
	fmpz_init(fb);
	fmpz_init(nm);
	fmpz_init(abm);
	fmpz_init(gcd);

	ulong loc = 0;
	ulong I = 2*N+1;
	double *r_sieve_array = new double[I];
	double *a_sieve_array = new double[I];
	for(ulong b = 1; b < MaxB; b++)
	{
		fmpz_set_ui(fb,b);
		fmpz_mul_ui(bm,m,b);
		rationalSieve(r_sieve_array, 2*N+1, RB, lRB, nRB, -N, bm);
		algebraicSieve(a_sieve_array, f, AB, lAB, nAB, 2*N+1, -N, b);
		for(ulong i = 0; i < I; i++)
		{
			slong a = i - N;
			fmpz_set_si(fa,a);
			fmpz_add(abm,bm,fa);
			fmpz_gcd(gcd,fa,fb);
			norm(nm,f,fa,fb);
			fmpz_abs(abm,abm);
			fmpz_abs(nm,nm);
			if(r_sieve_array[i] >= -5.0 && a_sieve_array[i] >= -5.0 && fmpz_is_one(gcd))
			{
				if(isSmooth(abm,RB,nRB) && isSmooth(nm,AB,nAB))
				{
					abPairs[loc] = MyPair(a,b);
					loc++;
					if(loc >= num) break;
				}
#ifdef PRINT_PROCESS
				cerr << "\r" << loc << "/" << num;
				cerr << "        ";
				cerr << b << "/" << MaxB;
				cerr << "        ";
				cerr.flush();
#endif
			}
		}
		if(loc >= num) break;
	}
#ifdef PRINT_PROCESS
	cerr << endl;
#endif
	assert(loc == num);
	num = loc;
	delete []r_sieve_array;
	delete []a_sieve_array;

	fmpz_clear(bm);
	fmpz_clear(fa);
	fmpz_clear(fb);
	fmpz_clear(nm);
	fmpz_clear(abm);
	fmpz_clear(gcd);
}

/**
 *	Form the matrix from the sieved (a,b) pairs
 *
 *	Each column is divided into four parts (sign, rational, algebraic, quadratic)
 *
 *	1. sign. Set to 0 if a+bm >= 0, to 1 otherwise.
 *	2. rational. (e1 mod 2, e2 mod 2, ... , eN mod 2), N is size of rational base.
 *		ei is power of i'th prime in a+bm.
 *	3. algebraic. (f1 mod 2, f2 mod 2, ... , fN mod 2), N is size of algebraic base.
 *		fi is power of i'th pair (p,r) in N(a,b).
 *	4. quadratic. (h1, h2, ... , hN), N is size of quadratic base.
 *		hi is 1 if legendre symbol(a+bs,q) == 0, and 0 otherwise.
 */
void formMatrix(nmod_mat_t mat, ulong I, ulong J, const fmpz_t m, const fmpz_poly_t f, const MyPair *abPairs,
				const ulong *RB, ulong nRB, const MyPair* AB, ulong nAB, const MyPair* QB, ulong nQB)
{
	fmpz_t fa,fb,A,r,s,bm;
	fmpz_init(fa);
	fmpz_init(fb);
	fmpz_init(A);
	fmpz_init(r);
	fmpz_init(s);
	fmpz_init(bm);
	nmod_mat_init(mat,I,J,2);
	for(ulong j = 0; j < J; j++)
	{
		slong a = abPairs[j].r;
		slong b = abPairs[j].p;
		fmpz_set_si(fa,a);
		fmpz_set_ui(fb,b);
		fmpz_mul_ui(bm,m,b);
		fmpz_add(A,bm,fa);
		*nmod_mat_entry_ptr(mat,0,j) = fmpz_cmp_si(A,0) >= 0? 0: 1;
		for(ulong i = 0; i < nRB; i++)
		{
			ulong p = RB[i];
			ulong e = 0;
			while(!fmpz_is_zero(A) && fmpz_divisible_si(A,p))
			{
				fmpz_fdiv_q_ui(A,A,p);
				e++;
			}
			*nmod_mat_entry_ptr(mat,i+1,j) = e % 2;
		}
#ifdef DEBUG
		assert(fmpz_is_pm1(A));
#endif
		norm(A,f,fa,fb);
		for(ulong i = 0; i < nAB; i++)
		{
			slong p = AB[i].p;
			slong r = AB[i].r;
			ulong e = 0;
			while(!fmpz_is_zero(A) && fmpz_divisible_si(A,p) && (a+b*r)%p == 0)
			{
				fmpz_fdiv_q_ui(A,A,p);
				e++;
			}
			*nmod_mat_entry_ptr(mat,i+1+nRB,j) = e % 2;
		}
#ifdef DEBUG
		assert(fmpz_is_pm1(A));
#endif
		for(ulong i = 0; i < nQB; i++)
		{
			slong s = QB[i].r;
			slong q = QB[i].p;
			int l = Leg(a+b*s,q);
			if((a+b*s)%q==0)
			{
				cout << q << " divides " << a+b*s << endl;
				cout << a << ' ' << b << endl;
				norm(A,f,fa,fb);
				fmpz_print(A); printf("\n");
			}
#ifdef DEBUG
			assert((a+b*s)%q != 0);
#endif
			*nmod_mat_entry_ptr(mat,i+1+nRB+nAB,j) = (l == 1? 0: 1);
		}
	}
	fmpz_clear(fa);
	fmpz_clear(fb);
	fmpz_clear(A);
	fmpz_clear(r);
	fmpz_clear(s);
	fmpz_clear(bm);
}


/**
 *	Solving linear system procedure: get a nonzero 0,1 solution of AX=O
 */
void solveMatrix(nmod_mat_t mat, ulong I, ulong J, int *vec)
{
	nmod_mat_rref(mat);
	ulong minI = 0;
	for(minI = 0; minI < I; minI++)
	{
		mp_limb_t b = nmod_mat_entry(mat,minI,minI);
		if(b == 0) break;
	}
	vec[minI] = 1;
	for(ulong i = minI+1; i < J; i++)
		vec[i] = 0;
	for(ulong i = 0; i < minI; i++)
		vec[i] = nmod_mat_entry(mat,i,minI);
}

/**
 *	Compute the square root of Prod(a+bm) for all the pairs (a,b) in an array
 */
void sqrtProductOfPairs(fmpz_t s, const MyPair *pairs, ulong num, const fmpz_t m)
{
	fmpz_t abm,bm,fa,r,e;
	fmpz_init(abm);
	fmpz_init(bm);
	fmpz_init(fa);
	fmpz_init(r);
	fmpz_init(e);
	fmpz_one(s);
	for(long i = 0; i < num; i++)
	{
		slong a = pairs[i].r;
		ulong b = pairs[i].p;
		fmpz_set_si(fa,a);
		fmpz_mul_ui(bm,m,b);
		fmpz_add(abm,bm,fa);
		fmpz_mul(s,s,abm);
	}
	fmpz_sqrtrem(r,e,s);
	/*Check that s is indeed a square number.*/
	assert(fmpz_is_zero(e));
	fmpz_set(s,r);

	fmpz_clear(abm);
	fmpz_clear(bm);
	fmpz_clear(fa);
	fmpz_clear(r);
	fmpz_clear(e);
}

/**
 *	Compute product of polynomials (a+bx) modular f for all the pairs (a,b) in the array.
 */
void productOfPairs(fmpz_poly_t sx, const MyPair *abPairs, ulong num, const fmpz_poly_t f, fmpz_t Nm)
{
	fmpz_t r,e,nm,fa,fb;
	fmpz_poly_t gx,q;
	fmpz_init(r);
	fmpz_init(e);
	fmpz_init(nm);
	fmpz_init(fa);
	fmpz_init(fb);
	fmpz_poly_init(gx);
	fmpz_poly_init(q);
	fmpz_poly_one(sx);
	fmpz_one(Nm);

	for(ulong i = 0; i < num; i++)
	{
		slong a = abPairs[i].r;
		ulong b = abPairs[i].p;
		fmpz_set_si(fa,a);
		fmpz_set_ui(fb,b);
		fmpz_poly_set_coeff_si(gx,0,a);
		fmpz_poly_set_coeff_si(gx,1,b);
		fmpz_poly_mul(sx,sx,gx);
		fmpz_poly_divrem(q,sx,sx,f);
		norm(nm,f,fa,fb);
		fmpz_mul(Nm,Nm,nm);
	}
	fmpz_sqrtrem(r,e,Nm);
	/*Check that Nm is indeed a square number.*/
	assert(fmpz_is_zero(e));
	fmpz_set(Nm,r);

	fmpz_clear(r);
	fmpz_clear(e);
	fmpz_clear(nm);
	fmpz_clear(fa);
	fmpz_clear(fb);
	fmpz_poly_clear(gx);
	fmpz_poly_clear(q);
}

/**
 *	Compute the square root of a polynomial sx modular polynomial f in Fp[x].
 */
ulong computeSquareRoot(const fmpz_poly_t sx, const fmpz_poly_t f, slong p, const fmpz_t m, const fmpz_t Nm)
{
#ifdef DEBUG
	assert(p > 1);
#endif

	slong d = fmpz_poly_degree(f);
	fmpz_t c,fp,q,q1,q2,s,s1,Nmp;
	mpz_t eq1,eq2,es,es1;
	mpz_init(eq1);
	mpz_init(eq2);
	mpz_init(es);
	mpz_init(es1);
	nmod_poly_t Fx, Sx, Sxp,lambda,omega,zeta,eta,pzeta1,pzeta2,omega2,NN;
	fmpz_init(c);
	fmpz_init_set_ui(fp,p);
	fmpz_init(q);
	fmpz_init(q1);
	fmpz_init(q2);
	fmpz_init(s);
	fmpz_init(s1);
	fmpz_init(Nmp);
	nmod_poly_init(Fx,p);
	nmod_poly_init(Sx,p);
	nmod_poly_init(Sxp,p);
	nmod_poly_init(lambda,p);
	nmod_poly_init(omega,p);
	nmod_poly_init(zeta,p);
	nmod_poly_init(eta,p);
	nmod_poly_init(pzeta1,p);
	nmod_poly_init(pzeta2,p);
	nmod_poly_init(omega2,p);
	nmod_poly_init(NN,p);

	for(int i = 0; i <= d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,f,i);
		nmod_poly_set_coeff_ui(Fx,i,fmpz_mod_ui(c,c,p));
	}
	for(int i = 0; i < d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,sx,i);
		nmod_poly_set_coeff_ui(Sx,i,fmpz_mod_ui(c,c,p));
	}
	fmpz_pow_ui(q,fp,d);
	fmpz_sub_ui(q,q,1);
	fmpz_fdiv_q_ui(q1,q,2);
	fmpz_fdiv_q_ui(q2,q,p-1);
	fmpz_get_mpz(eq1,q1);
	fmpz_get_mpz(eq2,q2);
	/*Check Sx is quadratic residual*/
	nmod_poly_powmod_mpz_binexp(Sxp,Sx,eq1,Fx);
#ifdef DEBUG
	assert(nmod_poly_is_one(Sxp));
#endif
	/********************************/
	slong r = 1;
	fmpz_set(s,q1);
	while(true)
	{
		if(!fmpz_divisible_si(s,2)) break;
		fmpz_fdiv_q_ui(s,s,2);
		r++;
	}
	fmpz_get_mpz(es,s);
	fmpz_add_ui(s1,s,1);
	fmpz_fdiv_q_ui(s1,s1,2);
	fmpz_get_mpz(es1,s1);
	nmod_poly_powmod_mpz_binexp(lambda,Sx,es,Fx);
	nmod_poly_powmod_mpz_binexp(omega,Sx,es1,Fx);
	selectNonResidual(eta,Fx,p,q1,d);
	nmod_poly_powmod_mpz_binexp(zeta,eta,es,Fx);

	while(!nmod_poly_is_one(lambda))
	{
		slong m = computeOrder(lambda,Fx);
		assert(r > m);
		nmod_poly_set(pzeta1,zeta);
		for(slong i = 0; i < r-m-1; i++)
			nmod_poly_powmod_ui_binexp(pzeta1,pzeta1,2,Fx);
		nmod_poly_powmod_ui_binexp(pzeta2,pzeta1,2,Fx);
		nmod_poly_mulmod(lambda,lambda,pzeta2,Fx);
		nmod_poly_mulmod(omega,omega,pzeta1,Fx);

		nmod_poly_t lambdadelta;
		nmod_poly_init(lambdadelta,p);
		nmod_poly_powmod_ui_binexp(omega2,omega,2,Fx);
		nmod_poly_mulmod(lambdadelta,lambda,Sx,Fx);
		assert(nmod_poly_equal(omega2,lambdadelta));
		nmod_poly_clear(lambdadelta);
	}
	/*Check that omega_n^2 = delta*/
	nmod_poly_powmod_ui_binexp(omega2,omega,2,Fx);
	assert(nmod_poly_equal(omega2,Sx));
	/*If the norm of omega is not equivalent to norm, negate it*/
	nmod_poly_powmod_mpz_binexp(NN,omega,eq2,Fx);
	assert(nmod_poly_degree(NN)==0);
	slong pp = nmod_poly_get_coeff_ui(NN,0);
	slong Nmmp = fmpz_mod_ui(Nmp,Nm,p);
	if(pp != Nmmp)
	{
		nmod_poly_neg(omega,omega);
		nmod_poly_powmod_mpz_binexp(NN,omega,eq2,Fx);
		assert(nmod_poly_degree(NN)==0);
		pp = nmod_poly_get_coeff_ui(NN,0);
		assert(pp == Nmmp);
	}
	/*Substitute m in*/
	slong M = nmod_poly_evaluate_nmod(omega,fmpz_mod_ui(Nmp,m,p));

	mpz_init(eq1);
	mpz_init(eq2);
	mpz_init(es);
	mpz_init(es1);
	fmpz_clear(c);
	fmpz_clear(fp);
	fmpz_clear(q);
	fmpz_clear(q1);
	fmpz_clear(q2);
	fmpz_clear(s);
	fmpz_clear(s1);
	fmpz_clear(Nmp);
	nmod_poly_clear(Fx);
	nmod_poly_clear(Sx);
	nmod_poly_clear(Sxp);
	nmod_poly_clear(lambda);
	nmod_poly_clear(omega);
	nmod_poly_clear(zeta);
	nmod_poly_clear(eta);
	nmod_poly_clear(pzeta1);
	nmod_poly_clear(pzeta2);
	nmod_poly_clear(omega2);
	nmod_poly_clear(NN);
	return M;
}

/**
 *	We use Chinese remainder theorem to compute x = phi(beta), first we need to find
 *	a set of primes p, such that their products is larger than x.
 */

/**
 *	First estimate an upperbound for x.
 */
void estimateUpperBoundForX(fmpz_t res, const fmpz_poly_t delta, const fmpz_t m, ulong d)
{
	fmpz_zero(res);
	fmpz_t pom,c;
	fmpz_init_set_ui(pom,1);
	fmpz_init(c);
	for(ulong i = 0; i < d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,delta,i);
		fmpz_abs(c,c);
		fmpz_sqrt(c,c);
		fmpz_addmul(res,c,pom);
		fmpz_mul(pom,pom,m);
	}
	fmpz_mul_ui(res,res,100);
	fmpz_clear(pom);
	fmpz_clear(c);
}

/**
 *	Then select the primes by repeatedly dividing the upperbound by larger and
 *	larger primes.
 */
bool selectPrimesCoverX(ulong *primes, ulong &nprimes, fmpz_t upperBound, ulong d, const fmpz_poly_t f)
{
	n_primes_t iter;
	n_primes_init(iter);
	n_primes_jump_after(iter,1000);
	nprimes = 0;
	while(fmpz_cmp_ui(upperBound,1) > 0)
	{
		ulong p = n_primes_next(iter);
		assert(p <= MaxPrime);
		assert(nprimes < MaxSelectedPrimes);
		if(irreducibleMod(f,p))
		{
			primes[nprimes++] = p;
			fmpz_cdiv_q_ui(upperBound,upperBound,p);
		}
	}
	n_primes_clear(iter);
	return fmpz_cmp_ui(upperBound,1)<=0;
}

/**
 *	For each selected prime p, compute square root of polynomial delta mod f mod p,
 *	then substitute m in the result polynomial, we will get a square root of
 *	phi(beta^2) mod p, which are the x mod p[i]'s used in the CRT.
 */
void computeSquareRoots(ulong *XmodPi, ulong *primes, ulong nprimes, fmpz_poly_t delta,
						const fmpz_poly_t f, const fmpz_t m, const fmpz_t Nm)
{
	for(ulong i = 0; i < nprimes; i++)
	{
		ulong p = primes[i];
		XmodPi[i] = computeSquareRoot(delta,f,p,m,Nm);
	}
}

/**
 *	A step in the CRT: compute the inverse of P[i] modular p[i], where P[i] is P/p[i],
 *	and P is the product of all the p[i]'s.
 */
void computePinvs(ulong *Pinv, const ulong *primes, ulong nprimes)
{
	fmpz_t fp,fq,f;
	fmpz_init(fp);
	fmpz_init(fq);
	fmpz_init(f);
	for(ulong i = 0; i < nprimes; i++)
	{
		ulong p = primes[i];
		fmpz_set_ui(fp,p);
		Pinv[i] = 1;
		for(ulong j = 0; j < nprimes; j++)
		{
			if(j == i) continue;
			ulong q = primes[j];
			fmpz_set_ui(fq,q);
			fmpz_invmod(f,fq,fp);
			ulong qinv = fmpz_get_ui(f);
			Pinv[i] = (Pinv[i]*qinv)%p;
		}
	}
	fmpz_clear(fp);
	fmpz_clear(fq);
	fmpz_clear(f);
}

void computePinvsModn(fmpz_t *Pinvmodn, const ulong *primes, ulong nprimes, const fmpz_t n)
{
	fmpz_t fp,f;
	fmpz_init(fp);
	fmpz_init(f);
	for(ulong i = 0; i < nprimes; i++)
	{
		fmpz_init(Pinvmodn[i]);
		ulong p = primes[i];
		fmpz_set_ui(fp,p);
		fmpz_invmod(Pinvmodn[i],fp,n);
	}
	fmpz_clear(fp);
	fmpz_clear(f);
}

/**
 *	For each p[i], compute (double) a[i]x[i]/p[i], where a[i] is the inverse of P[i] mod p[i].
 */
void computeAXoP(double *AXoP,const ulong *Pinv,const ulong *XmodPi,const ulong *primes,ulong nprimes)
{
	for(ulong i = 0; i < nprimes; i++)
		AXoP[i] = Pinv[i] * XmodPi[i] / (double) primes[i];
}

/**
 *	Summation of a[i]x[i]p[i].
 */
void sumOfAXPmodN(fmpz_t res, const ulong *Pinv, const ulong *XmodPi, const fmpz_t *pinvmodn, const fmpz_t Pmodn,
				  const ulong *primes, ulong nprimes, const fmpz_t n)
{
	fmpz_t s;
	fmpz_init(s);
	fmpz_zero(res);
	for(ulong i = 0; i < nprimes; i++)
	{
		fmpz_set(s,pinvmodn[i]);
		fmpz_mul(s,s,Pmodn);
		fmpz_mul_ui(s,s,XmodPi[i]);
		fmpz_mul_ui(s,s,Pinv[i]);
		fmpz_add(res,res,s);
		fmpz_mod(res,res,n);
	}
	fmpz_clear(s);
}

bool NFS(const fmpz_t n)
{
	/*--------------------Select polynomial-----------------------------------*/
	fmpz_poly_t f,delta;
	fmpz_t m, y, Nm, upperBoundOfX,r,x,Pmodn,rPmodn,xxmodn,yymodn;
	ulong d;
	fmpz_poly_init(delta);
#ifdef PRINT_PROCESS
	cout << "Selecting polynomial..." << endl;
#endif
	selectPolynomial(n,f,m,d);

#ifdef PRINT_MDF
	printf("m = "); fmpz_print(m); printf("\n");
	printf("d = %lu",d); printf("\n");
	printf("f(x) = "); fmpz_poly_print(f); printf("\n");
#endif
	
    /*--choose the bound for smoothness---------------------------------------*/
#ifdef PRINT_PROCESS
	cout << "Choosing smoothness bound..." << endl;
#endif
	ulong smoothBound = boundForSmoothness(d,n);
#ifdef PRINT_SMOOTH_BOUND
	cout << "Smoothness bound is " << smoothBound << endl;
#endif
	
	/*-Prepare the rational base----------------------------------------------*/
	ulong RB[MaxPrimeBufSize], nRB = 0;
	double lRB[MaxPrimeBufSize];
#ifdef PRINT_PROCESS
	cout << "Preparing the rational base..." << endl;
#endif
	prepareRationalBase(RB,lRB,nRB,smoothBound);
#ifdef PRINT_RATIONAL_BASE
	cout << "Rational base: " << endl;
	printListOfNumbers(RB,nRB,10);
#endif

	/*-Prepare the algebraic base---------------------------------------------*/
	MyPair AB[MaxPrimeBufSize];
	double lAB[MaxPrimeBufSize];
	ulong nAB = 0;
#ifdef PRITN_PROCESS
	cout << "Preparing the algebraic base..." << endl;
#endif
	prepareAlgebraicBase(AB,lAB,nAB,nRB,f);
#ifdef PRINT_ALGEBRAIC_BASE
	cout << "Algebraic base: " << endl;
	printListOfPairs(AB,nAB,10);
#endif

	/*-Prepare the quadratic base---------------------------------------------*/
	MyPair QB[MaxPrimeBufSize];
	ulong nQB = 0;
#ifdef PRINT_PROCESS
	cout << "Preparing the quadratic base..." << endl;
#endif
	prepareQuadraticBase(QB,nQB,smoothBound,smoothBound+20*log(smoothBound),f);
#ifdef PRINT_QUADRATIC_BASE
	cout << "Quadratic base: " << endl;
	printListOfPairs(QB,nQB,10);
#endif

	/*----------Sieve---------------------------------------------------------*/
	MyPair abPairs[2*MaxPrimeBufSize+1];
	ulong num = 2+nRB+nAB+nQB; /*Number of (a,b) pairs to search*/
	ulong N = smoothBound*20;
	time_t start = clock();
#ifdef PRINT_PROCESS
	cout << "Sieving for " << num << " (a,b) pairs in region [-" << N << "," << N << "]-[1," << MaxB << "]"  << endl;
#endif
	sieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs, num, N, m);
	//latticeSieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs, num, A, B, m);
#ifdef PRINT_PROCESS
	cout << "Time used by the sieving process: " << (clock() - start)/1000 << "ms" << endl;
#endif
#ifdef PRINT_SELECTED_ABPAIRS
	cout << "Selected smooth (a,b) pairs: " << endl;
	printListOfPairs(abPairs,num,10);
#endif

	/*---------Form matrix----------------------------------------------------*/
	ulong I = num-1, J = num;
	nmod_mat_t matrix;
#ifdef PRINT_PROCESS
	cout << "Forming the matrix..." << endl;
#endif
	formMatrix(matrix,I,J,m,f,abPairs,RB,nRB,AB,nAB,QB,nQB);
#ifdef PRINT_MATRIX
	cout << "The matrix is: " << endl;
	nmod_mat_print_pretty(matrix);
#endif

	/*-------Solve the linear system------------------------------------------*/
	int *vec = new int[J];
#ifdef PRINT_PROCESS
	cout << "Solving the " << I << " by " << J << " linear system..." << endl;
#endif
	solveMatrix(matrix,I,J,vec);
	select(abPairs,vec,J,num); /*Select the pairs corresponding to 1 in vec*/
	nmod_mat_clear(matrix);
	delete[] vec;

#ifdef PRINT_SELECTED_SQUARE_ABPAIRS
	cout << "The selected (a,b) pairs whose product is square in both Z and Z[X]:" << endl;
	printListOfPairs(abPairs,num,10);
#endif

	/*---------Calculate prod(a+bm)-------------------------------------------*/
#ifdef PRINT_PROCESS
	cout << "Computing prod(a+bm)..." << endl;
#endif
	sqrtProductOfPairs(y,abPairs,num,m);
#ifdef PRINT_PROD_ABM
	fmpz_print(y);cout << endl;
#endif
	fmpz_mod(y,y,n);

	/*---------Calculate prod(a+b theta)--------------------------------------*/
	/*Compute the product of the norm of (a,b) pairs, used to select
	  beta or -beta when computing square root of delta mod p*/
#ifdef PRINT_PROCESS
	cout << "Computing prod(a+b theta)/f(theta)..." << endl;
#endif
	productOfPairs(delta,abPairs,num,f,Nm);









	/*=================Calculate x = phi(beta)================================*/
	/*This is a big project, divide it into several parts*/
#ifdef PRINT_PROCESS
	cout << "Computing phi(beta) mod n..." << endl;
#endif
	/*1. Estimate an upper bound for x~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	estimateUpperBoundForX(upperBoundOfX,delta,m,d);
#ifdef PRINT_UPPDER_BOUND
	cout << "Upper bound for X is ";
	fmpz_print(upperBoundOfX); cout << endl;
#endif
	/************************************************************************/



	/*2. Select p_i that Prod p_i > x~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	ulong primes[MaxSelectedPrimes]; ulong nprimes;
#ifdef PRINT_PROCESS
	cout << "----Selecting primes p_i such that prod p_i > x..." << endl;
#endif
	if(!selectPrimesCoverX(primes,nprimes,upperBoundOfX,d,f)) return false;
#ifdef PRINT_PRIMES
	cout << "--------Selected primes: " << endl;
	printListOfNumbers(primes,nprimes,10);
#endif
	/************************************************************************/



	/*3. Compute x_i = x mod p_i~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	ulong XmodPi[MaxSelectedPrimes];
#ifdef PRINT_PROCESS
	printf("----Computing x_i = x mod p_i...\n");
#endif
	computeSquareRoots(XmodPi,primes,nprimes,delta,f,m,Nm);
#ifdef PRINT_XI
	printListOfNumbers(XmodPi,nprimes,10);
#endif
	/************************************************************************/



	/**************************/
	/*4. Compute x mod n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	ulong Pinv[MaxSelectedPrimes]; /*Inverse of P_i mod p_i, where P_i = P/p_i, P = Prod p_i*/
	double AXoP[MaxSelectedPrimes];/*a_i*x_i/p_i, where a_i = Pinv[i]*/
	fmpz_t Pinvmodn[MaxSelectedPrimes];
#ifdef PRINT_PROCESS
	cout << "----Computing x from (x mod p_i) by Chinese remainder..." << endl;
#endif
	computePinvs(Pinv,primes,nprimes);
	computeAXoP(AXoP,Pinv,XmodPi,primes,nprimes);
	computePinvsModn(Pinvmodn,primes,nprimes,n);
	/*Let z be the result of Chinese remainder, then x = z mod P,
	  and z = x + rP, so r = (z-x)/P, since x < P, r = Floor(z/P)*/
	fmpz_set_d(r,doublesum(AXoP,0,nprimes));
	productMod(Pmodn,primes,nprimes,n);
	sumOfAXPmodN(x,Pinv,XmodPi,Pinvmodn,Pmodn,primes,nprimes,n);
	fmpz_mul(rPmodn,r,Pmodn);
	fmpz_mod(rPmodn,rPmodn,n);
	fmpz_sub(x,x,rPmodn);
	if(fmpz_cmp_ui(x,0)<0) fmpz_add(x,x,n);
	/*There might be cases where x < 0, then the x obtained above is not
	 * the real one, subtract P from it and mod n again.*/
	fmpz_mul(xxmodn,x,x);
	fmpz_mod(xxmodn,xxmodn,n);
	fmpz_mul(yymodn,y,y);
	fmpz_mod(yymodn,yymodn,n);
	if(!fmpz_equal(xxmodn,yymodn))
	{
		fmpz_sub(x,x,Pmodn);
		fmpz_mul(xxmodn,x,x);
		fmpz_mod(xxmodn,xxmodn,n);
	}
	/************************************************************************/
	fmpz_clear(m);
	fmpz_clear(Nm);
	fmpz_clear(upperBoundOfX);
	fmpz_clear(rPmodn);
	fmpz_poly_clear(f);
	fmpz_poly_clear(delta);
	freeArrayOfFmpz(Pinvmodn,nprimes);



	fmpz_t xpy, xmy, f1, f2;
	fmpz_init(xpy);
	fmpz_init(xmy);
	fmpz_init(f1);
	fmpz_init(f2);
	fmpz_add(xpy,x,y);
	fmpz_sub(xmy,x,y);
	fmpz_gcd(f1,xpy,n);
	fmpz_gcd(f2,xmy,n);
	/*Finally, we get our x mod n and y mod n. Time to sum up.*/
#if(!ONLY_RESULT)
	printf("x mod n = "); fmpz_print(x); printf("\n");
	printf("y mod n = "); fmpz_print(y); printf("\n");
	/*Check square of x and y*/
	printf("x^2 mod n = "); fmpz_print(xxmodn); printf("\n");
	printf("y^2 mod n = "); fmpz_print(yymodn); printf("\n");
	printf("x + y = "); fmpz_print(xpy); printf("\n");
	printf("x - y = "); fmpz_print(xmy); printf("\n");
	printf("GCD(x+y,n) = "); fmpz_print(f1); printf("\n");
	printf("GCD(x-y,n) = "); fmpz_print(f2); printf("\n");
#endif
	fmpz_clear(xpy);
	fmpz_clear(xmy);
	fmpz_clear(y);
	fmpz_clear(xxmodn);
	fmpz_clear(yymodn);
	/*Return true if any of f1 and f2 is a proper factor of n*/
	if(fmpz_cmp_ui(f1,1)>0 && fmpz_cmp(f1,n)<0)
	{
		fmpz_t nof1;
		fmpz_init(nof1);
		fmpz_fdiv_q(nof1,n,f1);
		fmpz_print(n);
		printf(" = ");
		fmpz_print(f1);
		printf(" * ");
		fmpz_print(nof1);
#if(!ONLY_RESULT)
		printf("\n");
#endif
		fmpz_clear(nof1);
		fmpz_clear(f1);
		fmpz_clear(f2);
		return true;
	}
	if(fmpz_cmp_ui(f2,1)>0 && fmpz_cmp(f2,n)<0)
	{
		fmpz_t nof2;
		fmpz_init(nof2);
		fmpz_fdiv_q(nof2,n,f2);
		fmpz_print(n);
		printf(" = ");
		fmpz_print(f2);
		printf(" * ");
		fmpz_print(nof2);
#if(!ONLY_RESULT)
		printf("\n");
#endif
		fmpz_clear(nof2);
		fmpz_clear(f1);
		fmpz_clear(f2);
		return true;
	}
	return false;
}

void testNFS(ulong e)
{
	time_t t = clock();
	fmpz_t n;
	fmpz_init_set_ui(n,2);
	fmpz_pow_ui(n,n,e);
	fmpz_add_ui(n,n,1);
	printf("2^%lu+1 = \t",e);
	for(int i = 0; i < 50; i++)
		if(NFS(n)) break;
	printf("\t");
	t = clock() - t;
	printf("%ld\n",t);
	fmpz_clear(n);
}

void init()
{
	srand((unsigned)time(0));
}

int main(int argc, char *argv[])
{
	init();
#if(!ONLY_RESULT)
	fmpz_t n; fmpz_init(n);
	char *ptr = NULL;
	char num[] = "268435457";
	if(argc > 1) ptr = argv[1];
	else ptr = num;
	fmpz_set_str(n,ptr,10);
	printf("n = "); fmpz_print(n); printf("\n");
	int Tries = 50;
	for(int i = 0; i < Tries; i++)
	{
		printf("--------------------------------------------------------------------------------\n");
		printf("Trying for the %d time...\n",i+1);
		if(NFS(n)) break;
		printf("--------------------------------------------------------------------------------\n");
	}
	fmpz_clear(n);
#else
	testNFS(10);
	testNFS(23);
	testNFS(30);
	testNFS(41);
	testNFS(50);
	testNFS(61);
	testNFS(71);
	testNFS(82);
	testNFS(90);
	testNFS(101);
#endif
	return 0;
}
