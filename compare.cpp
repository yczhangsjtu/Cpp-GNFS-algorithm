#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <list>
#include <set>
#include <map>
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

#define DEBUG 1
#define PRINT_PROCESS 1
#define PRINT_SIEVE_PROCESS 1
#define SLOW_PRINT_SIEVE_PROCESS 0

#if(PRINT_PROCESS)
#define PRINT_MDF
#define PRINT_SMOOTH_BOUND
//#define PRINT_RATIONAL_BASE
//#define PRINT_ALGEBRAIC_BASE
//#define PRINT_QUADRATIC_BASE
//#define PRINT_SELECTED_ABPAIRS
#endif

using namespace std;

/*******************************************************************************
 *	Global Constants
 ******************************************************************************/
const int MaxPrimeBufSize = 30000;
const int MaxPrime = 200000;
const int MaxSelectedPrimes = 10000;
const int MaxB = 10240;
const ulong MaxT = 400000;
const int smoothfactor = 5;
const double threshold = 5.0;

const double abratio = 6000.0;
const double Afactor = 42;

/*
const double labratio = 2600.0;
const double lAfactor = 19;
*/
const double boundaryFactor = 1.0;

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
typedef set<MyPair> PairSet;
typedef map<MyPair,double> PairMap;

/**
 *	Pair used in the HashTable mapping integer pairs to doubles
 */
typedef struct PairDouble
{
	MyPair p;
	double v;
	PairDouble(MyPair a, double b){p=a;v=b;}
	PairDouble():p(){v=0.0;}
} PairDouble;

bool operator==(const MyPair &a, const MyPair &b)
{
	return a.r == b.r && a.p == b.p;
}

bool operator<(const MyPair &a, const MyPair &b)
{
	if(a.r < b.r) return true;
	if(a.r > b.r) return false;
	return a.p < b.p;
}

slong length(MyPair u)
{
	return u.p*u.p+u.r*u.r;
}

/*******************************************************************************
 *	Definition of HashTable class
 ******************************************************************************/
class HashTable{
public:
	ulong _T;
	HashTable(ulong T){
		if(T > MaxT) T = MaxT;
		_T = T;
		_data = new list<PairDouble>[T];
	}
	~HashTable(){
		if(size()/_T >= 3)
			cerr << endl << "Warning: Load factor = "
				 << (double)size()/_T << endl;
		delete []_data;
	}
	list<PairDouble> * _data;
	list<MyPair> _list;
	inline ulong size(){return _list.size();}
	slong hash(const MyPair &p)
	{
		return (p.r+p.p+p.r*(_T/2))%_T;
	}
	void set(const MyPair &p, double v)
	{
		/* Map pair p to value v, if not exist, insert it into the table */
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
		{
			if(iter->p==p)
			{
				iter->v=v;
				return;
			}
		}
		_data[s].push_back(PairDouble(p,v));
		_list.push_back(p);
	}
	inline void insert(const MyPair &p)
	{
		/* Insert pair p into table, used when the table is used as a set and
		 * the value is unimportant. */
		set(p,0.0);
	}
	void add(const MyPair &p, double v)
	{
		/* Add v to the value corresponding to p, if p does not exist, insert it
		 * and set the value to v. */
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
			if(iter->p==p)
			{
				iter->v+=v;
				return;
			}
		_data[s].push_back(PairDouble(p,v));
		_list.push_back(p);
	}
	double get(const MyPair &p)
	{
		/* Get the value of pair p. If not exist, return 0. */
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
			if(iter->p==p) return iter->v;
		return 0.0;
	}
	double find(const MyPair &p)
	{
		/* Return if p is contained in the table. */
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
			if(iter->p==p) return true;
		return false;
	}
};

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

bool isSmooth(MyPair p, const fmpz_t &m, const ulong *ps, ulong np)
{
	fmpz_t fa,fb,bm,abm,gcd;
	fmpz_init(fa);
	fmpz_init(fb);
	fmpz_init(bm);
	fmpz_init(abm);
	fmpz_init(gcd);
	fmpz_set_si(fa,p.r);
	fmpz_set_si(fb,p.p);
	fmpz_mul(bm,fb,m);
	fmpz_add(abm,bm,fa);
	fmpz_gcd(gcd,fa,fb);
	if(!fmpz_is_one(gcd)) return false;
	return isSmooth(abm,ps,np);
}

bool isSmooth(MyPair p, const fmpz_poly_t &f, const MyPair *ps, ulong np)
{
	fmpz_t fa,fb,nm,gcd;
	fmpz_init(fa);
	fmpz_init(fb);
	fmpz_init(nm);
	fmpz_init(gcd);

	fmpz_set_si(fa,p.r);
	fmpz_set_si(fb,p.p);
	fmpz_gcd(gcd,fa,fb);
	if(!fmpz_is_one(gcd)) return false;

	norm(nm,f,fa,fb);
	return isSmooth(nm,ps,np);
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
 * Test if a polynomial is irreducible.
 */
bool testPolynomial(const fmpz_poly_t f)
{
	n_primes_t iter;
	n_primes_init(iter);
	n_primes_jump_after(iter,1000);
	while(true)
	{
		ulong p = n_primes_next(iter);
		if(irreducibleMod(f,p)) return true;
		if(p > MaxPrime) break;
	}
	n_primes_clear(iter);
	return false;
}

/* Tools about the lattice.****************************************************/

/*Gaussian lattice reduce: given two vectors definint a lattice,
 * find the smallest basis that generates the lattice. */
#if 0
void gaussianLatticeReduce(MyPair &u, MyPair &v)
{
	slong lu = length(u), lv = length(v);
	do{
		if(lu > lv)
		{
			MyPair tmp = u;
			u = v;
			v = tmp;
			lu = length(u);
			lv = length(v);
		}
		if(lu == 0) return;
		slong k = (u.p*v.p+u.r*v.r)/lu;
		slong k1 = k + 1;
		slong l = (v.p-k*u.p)*(v.p-k*u.p)+(v.r-k*u.r)*(v.r-k*u.r);
		slong l1 = (v.p-k1*u.p)*(v.p-k1*u.p)+(v.r-k1*u.r)*(v.r-k1*u.r);
		if(l > l1)
		{
			v.p -= k1*u.p;
			v.r -= k1*u.r;
		}
		else
		{
			v.p -= k*u.p;
			v.r -= k*u.r;
		}
		lv = length(v);
	} while(lu > lv);
}
#else
void gaussianLatticeReduce(MyPair &u, MyPair &v)
{
	v.r %= u.r;
}
#endif

/**
 *	Given C and D, and two vectors s, t, find the bound E1,E2 for s, F1,F2 for t, such that
 *	the region {es+ft | E1<=e<=E2, F1<=f<=F2} covers the region {(c,d) | -C<=c<=C, -D<=d<=D}
 *	as 'good' as possible.
 *
 *	Here the algorithm is simply taken as: stretch the vectors s,t until reaching the bound
 */
void getBound(slong &E1, slong &E2, slong &F1, slong &F2, slong C, slong D, MyPair s, MyPair t)
{
	if(s.r == 0) E1 = D/abs(s.p);
	else if(s.p == 0) E1 = C/abs(s.r);
	else
	{
		E1 = C/abs(s.r);
		if(D/abs(s.p) < E1) E1 = D/abs(s.p);
	}
	E1 *= boundaryFactor;
	E2 = E1;
	E1 = -E1;
	if(t.r == 0) F2 = D/abs(t.p);
	else if(t.p == 0) F2 = C/abs(t.r);
	else
	{
		F2 = C/abs(t.r);
		if(D/abs(t.p) < F2) F2 = D/abs(t.p);
	}
	F2 *= boundaryFactor;
	F1 = 1;
#ifdef DEBUG
	assert(E2 > 0 && F2 > 0);
#endif
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
	fmpz_root(m,n,d); /* Set m to the d'th root of n rounded to floor. */
	/* Loop until proper f,m are found. */
	while(true)
	{
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
		fmpz_poly_print(f);
		cout << endl;
		/* f is required to be monic, nonzero constant term, and the constant term is not 1. */
		if(fmpz_is_one(lead) && !fmpz_is_zero(tail) && !fmpz_is_one(tail) && testPolynomial(f)) break;
		do
		{
			fmpz_sub_ui(m,m,1);
			fmpz_print(m);
			cout << endl;
		} while(fmpz_divisible(n,m));
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

/**
 *	Sieving: Lattice sieve in the rational part.
 */
void latticeRationalSieve(HashTable &pairs, const ulong *RB, const double *lRB,
	ulong iRB, ulong nRB, MyPair u, MyPair v, slong C, slong D, double logq, const fmpz_t m,
	int &count)
{
	time_t start = clock();
	HashTable result(300*C*log(D));
	count += clock()-start;
	slong im = fmpz_get_si(m);
	/* Accumulate log(p) for each p < q in the base. */
	for(ulong j = 0; j < iRB; j++)
	{
		slong p = RB[j];
		double lp = lRB[j];
		slong h = (u.r+u.p*im)%p;
		slong k = (v.r+v.p*im)%p;
		/* For each p, find all the (c,d) pairs such that
		 * (c*h+d*k) mod p == 0 */
		if(k)
		{
			fmpz_t ki,fk,fp;
			fmpz_init(ki);
			fmpz_init_set_si(fk,k);
			fmpz_init_set_si(fp,p);
			fmpz_invmod(ki,fk,fp);
			for(slong c = -C; c <= C; c++)
			{
				slong ch = c*h;
				slong kinv = fmpz_get_si(ki);
				slong d = -ch * kinv;
#ifdef DEBUG
				assert((c*h+d*k)%p==0);
#endif
				if(d > 0) d -= (d/p)*p;
				if(d < 0) d += ((-d-1)/p+1)*p;
#ifdef DEBUG
				assert(d >= 0);
#endif
				slong pvr = p*v.r;
				slong pvp = p*v.p;
				for(;d <= D; d+=p)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#ifdef DEBUG
					assert((a+b*im)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					result.add(MyPair(a,b),lp);
				}
			}
			fmpz_clear(ki);
			fmpz_clear(fk);
			fmpz_clear(fp);
		}
		else if(h)
		{
			slong c = -(C/p)*p;
			for(;c <= C; c+=p)
			{
				slong d = 0;
				for(; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#ifdef DEBUG
					assert((a+b*im)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					result.add(MyPair(a,b),lp);
				}
			}
		}
		else
		{
			for(slong c = -C; c <= C; c++)
			{
				slong d = 0;
				for(; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#ifdef DEBUG
					assert((a+b*im)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					result.add(MyPair(a,b),lp);
				}
			}
		}
	}
	for(list<MyPair>::iterator iter = result._list.begin(); iter != result._list.end(); iter++)
	{
		fmpz_t bm,fa,fb,abm,gcd;
		fmpz_init(fa);
		fmpz_init(fb);
		fmpz_init(bm);
		fmpz_init(abm);
		fmpz_init(gcd);

		slong a = iter->r;
		slong b = iter->p;
		fmpz_set_si(fa,a);
		fmpz_set_si(fb,b);
		fmpz_mul_si(bm,m,b);
		fmpz_add(abm,bm,fa);
		fmpz_gcd(gcd,fa,fb);
		if(!fmpz_is_one(gcd)) continue;
		fmpz_abs(abm,abm);
		double labm = fmpz_dlog(abm);
		if(result.get(*iter)+logq - labm < -threshold)
			continue;
		pairs.insert(*iter);

		fmpz_clear(bm);
		fmpz_clear(fa);
		fmpz_clear(fb);
		fmpz_clear(abm);
		fmpz_clear(gcd);
	}
}

/**
 *	Sieving: Lattice sieve in the algebraic part.
 */
void latticeAlgebraicSieve(HashTable &abPairs, ulong &loc, HashTable &pairs, const fmpz_poly_t f,
	const MyPair *AB, const double *lAB, ulong iAB, ulong nAB, MyPair u, MyPair v, slong C, slong D,
	const ulong *RB, ulong nRB, const fmpz_t m, int &count)
{
	time_t start = clock();
	HashTable result(300*C*log(D));
	count += clock()-start;
	for(ulong j = 0; j < nAB; j++)
	{
		slong p = AB[j].p;
		double lp = lAB[j];
		slong r = AB[j].r;
		slong h = (u.r+u.p*r)%p;
		slong k = (v.r+v.p*r)%p;
		if(k)
		{
			fmpz_t ki,fk,fp;
			fmpz_init(ki);
			fmpz_init_set_si(fk,k);
			fmpz_init_set_si(fp,p);
			fmpz_invmod(ki,fk,fp);
			for(slong c = -C; c <= C; c++)
			{
				slong ch = c*h;
				slong kinv = fmpz_get_si(ki);
				slong d = -ch * kinv;
#ifdef DEBUG
				assert((c*h+d*k)%p==0);
#endif
				if(d > 0) d -= (d/p)*p;
				if(d < 0) d += ((-d-1)/p+1)*p;
#ifdef DEBUG
				assert(d >= 0);
#endif
				for(;d <= D; d+=p)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#ifdef DEBUG
					assert((a+b*r)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					if(!pairs.find(MyPair(a,b))) continue;
					result.add(MyPair(a,b),lp);
				}
			}
			fmpz_clear(ki);
			fmpz_clear(fk);
			fmpz_clear(fp);
		}
		else if(h)
		{
			slong c = -(C/p)*p;
			for(;c <= C; c+=p)
			{
				slong d = 0;
				for(; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#ifdef DEBUG
					assert((a+b*r)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					if(!pairs.find(MyPair(a,b))) continue;
					result.add(MyPair(a,b),lp);
				}
			}
		}
		else
		{
			for(slong c = -C; c <= C; c++)
			{
				slong d = 0;
				for(; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#ifdef DEBUG
					assert((a+b*r)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					if(!pairs.find(MyPair(a,b))) continue;
					result.add(MyPair(a,b),lp);
				}
			}
		}
	}
	for(list<MyPair>::iterator iter = result._list.begin(); iter != result._list.end(); iter++)
	{
		count++;
		fmpz_t fa,fb,nm,gcd,bm,abm;
		fmpz_init(fa);
		fmpz_init(fb);
		fmpz_init(nm);
		fmpz_init(bm);
		fmpz_init(abm);
		fmpz_init(gcd);

		slong a = iter->r;
		slong b = iter->p;
		fmpz_set_si(fa,a);
		fmpz_set_si(fb,b);
		fmpz_mul_si(bm,m,b);
		fmpz_add(abm,bm,fa);
		norm(nm,f,fa,fb);
		fmpz_gcd(gcd,fa,fb);
		if(!fmpz_is_one(gcd)) continue;
		fmpz_abs(nm,nm);

		double lnm = fmpz_dlog(nm);
		if(result.get(*iter) - lnm < -threshold) continue;
		if(!isSmooth(nm,AB,nAB)) continue;
		if(!isSmooth(abm,RB,nRB)) continue;
		if(pairs.find(*iter))
		{
			abPairs.insert(*iter);
			loc = abPairs.size();
		}

		fmpz_clear(fa);
		fmpz_clear(fb);
		fmpz_clear(nm);
		fmpz_clear(gcd);
	}
}

/* Traditional method: used if lattice sieve failed to find enough number of pairs. */
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

/**
 *	Sieving: The main procedure.
 */
void latticeSieve(const fmpz_poly_t f, const ulong *RB, const double *lRB, ulong nRB,
		   const MyPair *AB, const double *lAB, ulong nAB, PairSet &abPairs, slong A, slong B, fmpz_t m)
{
	ulong loc = 0;
	/* Use the HashTable to store the found (a,b) pairs to prevent repeat. */
	HashTable abpairs(A);
	/* Loop for each special-q. */
	int count1 = 0, count2 = 0, count3 = 0;
	for(ulong i = nRB/smoothfactor; i < nRB; i++)
	{
		HashTable pairs(A);

		slong q = RB[i];
		double logq = lRB[i];
		slong im = fmpz_get_si(m);
		MyPair u(q,0), v(im,-1);
		gaussianLatticeReduce(u,v);
		slong C1,C2,D1,D2;
		getBound(C1,C2,D1,D2,A,B,u,v);
		slong C = C2;
		slong D = D2;
		//cout << u.r*C << ' ' << u.p*C << ' ' << v.r*D << ' ' << v.p*D << endl;
		latticeRationalSieve(pairs, RB, lRB, i, nRB, u, v, C, D, logq, m, count1);
		latticeAlgebraicSieve(abpairs, loc, pairs, f, AB, lAB, i, nAB, u, v, C, D, RB, nRB, m, count2);

#if(PRINT_PROCESS && PRINT_SIEVE_PROCESS)
		loc = abpairs.size();
		cerr << "\r" << loc << "       ";
		cerr << i << "/" << nRB;
		cerr << "                       "; cerr.flush();
#endif
#if(PRINT_PROCESS && SLOW_PRINT_SIEVE_PROCESS)
		cerr << loc << "       ";
		cerr << i << "/" << nRB;
		cerr << endl;
#endif
	}
#if(PRINT_PROCESS && PRINT_SIEVE_PROCESS)
	cerr << endl;
#endif
	slong k = 0;
	for(list<MyPair>::iterator iter = abpairs._list.begin(); iter!=abpairs._list.end(); iter++)
	{
		abPairs.insert(*iter);
	}
	cout << count1 << endl;
	cout << count2 << endl;
	cout << count3 << endl;
}

void sieve(fmpz_poly_t f, const ulong *RB, const double *lRB, ulong nRB,
		   const MyPair *AB, const double *lAB, ulong nAB, PairSet &abPairs, slong A, slong B, fmpz_t m)
		   
{
	fmpz_t bm,fa,fb,nm,abm,gcd;
	fmpz_init(bm);
	fmpz_init(fa);
	fmpz_init(fb);
	fmpz_init(nm);
	fmpz_init(abm);
	fmpz_init(gcd);

	ulong loc = 0;
	ulong I = 2*A+1;
	double *r_sieve_array = new double[I];
	double *a_sieve_array = new double[I];
	int count = 0;
	for(ulong b = 1; b < B; b++)
	{
		fmpz_set_ui(fb,b);
		fmpz_mul_ui(bm,m,b);
		rationalSieve(r_sieve_array, 2*A+1, RB, lRB, nRB, -A, bm);
		algebraicSieve(a_sieve_array, f, AB, lAB, nAB, 2*A+1, -A, b);
		for(ulong i = 0; i < I; i++)
		{
			slong a = i - A;
			fmpz_set_si(fa,a);
			fmpz_add(abm,bm,fa);
			fmpz_gcd(gcd,fa,fb);
			norm(nm,f,fa,fb);
			fmpz_abs(abm,abm);
			fmpz_abs(nm,nm);
			count++;
			if(r_sieve_array[i] >= -5.0 && a_sieve_array[i] >= -5.0 && fmpz_is_one(gcd))
			{
				if(isSmooth(abm,RB,nRB))
				{
					if(isSmooth(nm,AB,nAB))
					{
						abPairs.insert(MyPair(a,b));
						loc++;
					}
				}
#if(PRINT_PROCESS && PRINT_SIEVE_PROCESS)
				cerr << "\r" << loc << "       ";
				cerr << b << "/" << B;
				cerr << "        ";
				cerr.flush();
#endif
			}
		}
#if(PRINT_PROCESS && SLOW_PRINT_SIEVE_PROCESS)
		cerr << loc;
		cerr << "        ";
		cerr << b << "/" << B;
		cerr << "        ";
		cerr << endl;
#endif
	}
#ifdef PRINT_PROCESS
	cerr << endl;
#endif
	delete []r_sieve_array;
	delete []a_sieve_array;

	fmpz_clear(bm);
	fmpz_clear(fa);
	fmpz_clear(fb);
	fmpz_clear(nm);
	fmpz_clear(abm);
	fmpz_clear(gcd);
	cout << count << endl;
}

void NFS(const fmpz_t n)
{
	/*--------------------Select polynomial-----------------------------------*/
	fmpz_poly_t f;
	fmpz_t m;
	ulong d;
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
	//MyPair abPairs[2*MaxPrimeBufSize+1];
	PairSet abPairs1, abPairs2, abPairs3;
	ulong num = 2+nRB+nAB+nQB; /*Number of (a,b) pairs to search*/
	slong A = smoothBound*Afactor, B = A/abratio;
	cout << "Sieving for " << num << " (a,b) pairs in region [" << -A << "," << A << "]-[1," << B << "]"  << endl;
	time_t t1, t2, t3;
	time_t start = clock();
	sieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs1, A, B, m);
	cout << "Find " << abPairs1.size() << " pairs." << endl;
	cout << "Time used by the sieving process: " << (clock() - start)/1000 << "ms" << endl;
	t1 = clock() - start;

/*
	start = clock();
	sieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs2, A, B, m);
	cout << "Find " << abPairs2.size() << " pairs." << endl;
	cout << "Time used by the sieving process: " << (clock() - start)/1000 << "ms" << endl;
	t2 = clock() - start;
*/

	start = clock();
	A = smoothBound*Afactor, B = A/abratio;
	cout << "Sieving for " << num << " (a,b) pairs in region [" << -A << "," << A << "]-[1," << B << "]"  << endl;
	latticeSieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs3, A, B, m);
	cout << "Find " << abPairs3.size() << " pairs." << endl;
	cout << "Time used by the sieving process: " << (clock() - start)/1000 << "ms" << endl;
	t3 = clock() - start;

	fmpz_clear(m);
	fmpz_poly_clear(f);
}

void init()
{
	srand((unsigned)time(0));
}

int main(int argc, char *argv[])
{
	init();
	fmpz_t n; fmpz_init(n);
	char *ptr = NULL;
	char num[] = "268435457";
	if(argc > 1) ptr = argv[1];
	else ptr = num;
	fmpz_set_str(n,ptr,10);
	printf("n = "); fmpz_print(n); printf("\n");
	NFS(n);
	fmpz_clear(n);
	return 0;
}
