#include "GNFS.h"
#include "poly.h"

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

void getMaxCoeff(fmpz_t m, const fmpz_poly_t f)
{
	ulong d = fmpz_poly_degree(f);
	fmpz_t c;
	fmpz_init(c);
	fmpz_poly_get_coeff_fmpz(m,f,0);
	for(slong i = 1; i <= d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,f,i);
		if(fmpz_cmp(c,m) > 0) fmpz_set(m,c);
	}
	fmpz_clear(c);
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
