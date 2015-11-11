#include <algorithm>
#include "GNFS.h"
#include "checkredundancy.h"

using namespace std;
int main(int argc, char *argv[])
{
	if(argc < 2)
	{
		cerr << "Usage: checkredundancy inputfile" << endl;
		exit(-1);
	}
	FILE *input = fopen(argv[1],"r");
	if(!input) perror(argv[1]);

	fmpz_t n, m;
	fmpz_poly_t f;
	ulong RB[MaxPrimeBufSize], nRB = 0, nAB = 0, nQB = 0, num = 0;
	MyPair AB[MaxPrimeBufSize], QB[MaxPrimeBufSize], abPairs[2*MaxPrimeBufSize+1];
	int *vec = NULL;

	fmpz_init(n);
	fmpz_init(m);
	fmpz_poly_init(f);
	fmpz_fread(input,n);
	fmpz_fread(input,m);
	fmpz_poly_fread(input,f);
	fscanf(input,"%lu",&nRB);
	for(slong i = 0; i < nRB; i++)
		fscanf(input,"%lu",&RB[i]);
	fscanf(input,"%lu",&nAB);
	for(slong i = 0; i < nAB; i++)
		fscanf(input,"%ld%ld",&AB[i].r,&AB[i].p);
	fscanf(input,"%lu",&nQB);
	for(slong i = 0; i < nQB; i++)
		fscanf(input,"%ld%ld",&QB[i].r,&QB[i].p);
	fscanf(input,"%lu",&num);
	for(slong i = 0; i < num; i++)
		fscanf(input,"%ld%ld",&abPairs[i].r,&abPairs[i].p);
	
	sort(&abPairs[0],&abPairs[num]);

	fmpz_clear(n);
	fmpz_clear(m);
	fmpz_poly_clear(f);

	for(slong i = 0; i < num-1; i++)
	{
		if(abPairs[i].r == abPairs[i+1].r && abPairs[i].p == abPairs[i+1].p)
		{
			fprintf(stderr,"(%ld,%ld) is repeated.\n",abPairs[i].r,abPairs[i].p);
			exit(-1);
		}
	}
	for(slong i = 0; i < num; i++)
	{
		fmpz_t fa, fb, gcd;
		fmpz_init_set_si(fa,abPairs[i].r);
		fmpz_init_set_si(fb,abPairs[i].p);
		fmpz_init(gcd);
		fmpz_gcd(gcd,fa,fb);
		fmpz_abs(gcd,gcd);
		if(!fmpz_is_one(gcd))
		{
			fprintf(stderr,"(%ld,%ld) is not coprime.\n",abPairs[i].r,abPairs[i].p);
			exit(-1);
		}
		fmpz_clear(fa);
		fmpz_clear(fb);
		fmpz_clear(gcd);
	}

	fclose(input);
	return 0;
}
