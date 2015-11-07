#include "GNFS.h"

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

