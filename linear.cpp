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
	ulong *fs = new ulong[I];
	ulong *bufi = new ulong[J];
	ulong *bufj = new ulong[J];
	ulong i0 = 0, j0 = 0, pbuf = 0;
	while(i0 < I && j0 < J)
	{
		mp_limb_t b = nmod_mat_entry(mat,i0,j0);
		if(b == 1)
		{
			fs[i0] = j0;
			i0++;
			j0++;
		}
		else
		{
			bufi[pbuf] = i0;
			bufj[pbuf] = j0++;
			pbuf++;
		}
	}
	while(j0 < J)
	{
		bufi[pbuf] = i0;
		bufj[pbuf] = j0++;
	}
	slong k = rand()%pbuf;
	slong t = bufj[k];
	slong s = bufi[k];
	for(ulong i = 0; i < J; i++)
		vec[i] = 0;
	vec[t] = 1;
	for(slong i = 0; i < s; i++)
	{
		if(nmod_mat_entry(mat,i,t) == 1)
			vec[fs[i]] = 1;
	}
	
	delete []fs;
	delete []bufi;
	delete []bufj;
}

int main(int argc, char *argv[])
{
	if(argc < 3)
	{
		cerr << "Usage: linear inputfile outputfile" << endl;
		exit(-1);
	}
	FILE *input = fopen(argv[1],"r");
	if(!input) perror(argv[1]);
	FILE *output = fopen(argv[2],"w");
	if(!output) perror(argv[2]);

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

	ulong I = num-1, J = num;
	nmod_mat_t matrix;
	formMatrix(matrix,I,J,m,f,abPairs,RB,nRB,AB,nAB,QB,nQB);
	vec = new int[J];
	solveMatrix(matrix,I,J,vec);
	select(abPairs,vec,J,num); /*Select the pairs corresponding to 1 in vec*/
	nmod_mat_clear(matrix);
	delete[] vec;

	fmpz_fprint(output,n); fprintf(output,"\n");
	fmpz_fprint(output,m); fprintf(output,"\n");
	fmpz_poly_fprint(output,f); fprintf(output,"\n");
	fprintf(output,"%lu\n",num); fprintf(output,"\n");
	printListOfPairs(output,abPairs,num,5);

	fmpz_clear(n);
	fmpz_clear(m);
	fmpz_poly_clear(f);

	fclose(input);
	fclose(output);
	return 0;
}
