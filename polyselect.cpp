#include "GNFS.h"

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
		/* f is required to be monic, nonzero constant term, and the constant term is not 1. */
		if(fmpz_is_one(lead) && !fmpz_is_zero(tail) && !fmpz_is_one(tail) && testPolynomial(f)) break;
		fmpz_sub_ui(m,m,1);
		while(fmpz_divisible(n,m))
		{
			fmpz_sub_ui(m,m,1);
		}
	}
	fmpz_clear(lead);
	fmpz_clear(tail);
	fmpz_clear(N);
	fmpz_clear(Nmm);
}

