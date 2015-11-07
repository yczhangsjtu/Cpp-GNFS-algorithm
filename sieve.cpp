#include "GNFS.h"

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
			if(r_sieve_array[i] < -5.0 || a_sieve_array[i] < -5.0) continue;
			slong a = i - N;
			fmpz_set_si(fa,a);
			fmpz_add(abm,bm,fa);
			fmpz_gcd(gcd,fa,fb);
			norm(nm,f,fa,fb);
			fmpz_abs(abm,abm);
			fmpz_abs(nm,nm);
			if(!fmpz_is_one(gcd)) continue;
			if(isSmooth(abm,RB,nRB) && isSmooth(nm,AB,nAB))
			{
				abPairs[loc] = MyPair(a,b);
				loc++;
				if(loc >= num) break;
			}
#if(PRINT_PROCESS && PRINT_SIEVE_PROCESS)
			cerr << "\r" << loc << "/" << num;
			cerr << "        ";
			cerr << b << "/" << MaxB;
			cerr << "        ";
			cerr.flush();
#endif
		}
#if(PRINT_PROCESS && SLOW_PRINT_SIEVE_PROCESS)
		cerr << loc << "/" << num;
		cerr << "        ";
		cerr << b << "/" << MaxB;
		cerr << "        ";
		cerr << endl;
#endif
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

