#include "GNFS.h"
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
	slong t = minI;
	for(ulong i = minI; i < J; i++)
		vec[i] = 0;
	vec[t] = 1;
	for(slong i = 0; i < minI; i++)
		vec[i] = nmod_mat_entry(mat,i,t);
}

