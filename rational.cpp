#include "rational.h"

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

