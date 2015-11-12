#include <mpi.h>
#include <unistd.h>
#include "GNFS.h"
#include "algebraic.h"
#include "rational.h"
#include "sieve.h"
#include "poly.h"

int MaxPrime = DefaultMaxPrime;
int Afactor = DefaultAfactor;
long int MaxB = 20240;
double threshold = DefaultThreshold;

void sieve(fmpz_poly_t f, const ulong *RB, const double *lRB, ulong nRB,
		   const MyPair *AB, const double *lAB, ulong nAB, MyPair *abPairs,
		   ulong num, slong N, fmpz_t m,
		   int mynode, int totalnodes, MPI_Status *status)
{
	fmpz_t bm,fa,fb,nm,abm,gcd;
	fmpz_init(bm);
	fmpz_init(fa);
	fmpz_init(fb);
	fmpz_init(nm);
	fmpz_init(abm);
	fmpz_init(gcd);

	int end = 0;
	ulong loc = 0;
	ulong I = 2*N+1;
	double *r_sieve_array = new double[I];
	double *a_sieve_array = new double[I];
	for(ulong b = mynode+1; b < MaxB; b+=totalnodes)
	{
		fmpz_set_ui(fb,b);
		fmpz_mul_ui(bm,m,b);
		rationalSieve(r_sieve_array, 2*N+1, RB, lRB, nRB, -N, bm);
		algebraicSieve(a_sieve_array, f, AB, lAB, nAB, 2*N+1, -N, b);
		if(mynode) loc = 0;
		for(ulong i = 0; i < I; i++)
		{
			if(r_sieve_array[i] < -threshold || a_sieve_array[i] < -threshold) continue;
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
			if(mynode == 0)
			{
				std::cerr << "\r" << loc << "/" << num;
				std::cerr << "        ";
				std::cerr << b << "/" << MaxB;
				std::cerr << "        ";
				std::cerr.flush();
			}
#endif
		}
		if(mynode == 0)
		{
#if(PRINT_PROCESS && SLOW_PRINT_SIEVE_PROCESS)
			std::cerr << loc << "/" << num;
			std::cerr << "        ";
			std::cerr << b << "/" << MaxB;
			std::cerr << "        ";
			std::cerr << std::endl;
#endif
			MyPair *recvbuf = new MyPair[I];
			for(int k = 1; k < totalnodes; k++)
			{
				int recvsize;
				MPI_Recv(&recvsize,1,MPI_UNSIGNED_LONG,k,0,MPI_COMM_WORLD,status);
				MPI_Recv(recvbuf,I*2,MPI_LONG,k,0,MPI_COMM_WORLD,status)/2;
#if(DEBUG)
				std::cerr << "Received " << recvsize << " pairs from node " << k << std::endl;
#endif
				for(int i = 0; i < recvsize; i++)
				{
					abPairs[loc] = recvbuf[i];
					loc++;
					if(loc >= num) break;
				}
				if(loc >= num) break;
				end = 0;
				MPI_Send(&end,1,MPI_INT,k,0,MPI_COMM_WORLD);
			}
			delete []recvbuf;
		}
		else
		{
#if(DEBUG)
			std::cerr << "Found " << loc << " pairs in b = " << b << "." << std::endl;
#endif
			MPI_Send(&loc,1,MPI_UNSIGNED_LONG,0,0,MPI_COMM_WORLD);
			MPI_Send(abPairs,loc*2,MPI_LONG,0,0,MPI_COMM_WORLD);
			MPI_Recv(&end,1,MPI_INT,0,0,MPI_COMM_WORLD,status);
			if(end) break;
		}
		if(mynode == 0 && loc >= num) break;
		if(mynode && end) break;
	}
	delete []r_sieve_array;
	delete []a_sieve_array;

	if(mynode == 0)
	{
#if(PRINT_PROCESS && PRINT_SIEVE_PROCESS)
		std::cerr << std::endl;
#endif
		assert(loc >= num);
		num = loc;
		end = 1;
		for(int k = 1; k < totalnodes; k++)
			MPI_Send(&end,1,MPI_INT,k,0,MPI_COMM_WORLD);
	}

	fmpz_clear(bm);
	fmpz_clear(fa);
	fmpz_clear(fb);
	fmpz_clear(nm);
	fmpz_clear(abm);
	fmpz_clear(gcd);
}

using namespace std;
int main(int argc, char *argv[])
{
	int mynode, totalnodes;
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&totalnodes);
	MPI_Comm_rank(MPI_COMM_WORLD,&mynode);


	if(argc < 3)
	{
		cerr << "Usage: sieve inputfile outputfile" << endl;
		exit(-1);
	}
	FILE *input = fopen(argv[1],"r");
	if(!input) perror(argv[1]);
	FILE *output = NULL;
	if(mynode == 0)
	{
		output = fopen(argv[2],"w");
		if(!output) perror(argv[2]);
	}

	int ch;
	while((ch = getopt(argc-2,argv+2,"a:")) != -1)
	{
		switch(ch)
		{
		case 'a': Afactor = atoi(optarg); break;
		}
	}

	fmpz_t n, m;
	fmpz_poly_t f;
	ulong RB[MaxPrimeBufSize], nRB = 0, nAB = 0, nQB = 0, smoothBound;
	MyPair AB[MaxPrimeBufSize], QB[MaxPrimeBufSize];
	double lRB[MaxPrimeBufSize], lAB[MaxPrimeBufSize];

	fmpz_init(n);
	fmpz_init(m);
	fmpz_poly_init(f);
	fmpz_fread(input,n);
	fmpz_fread(input,m);
	fscanf(input,"%lu",&smoothBound);
	fmpz_poly_fread(input,f);
	fscanf(input,"%lu",&nRB);
	for(slong i = 0; i < nRB; i++)
	{
		fscanf(input,"%lu",&RB[i]);
		lRB[i] = log(RB[i]);
	}
	fscanf(input,"%lu",&nAB);
	for(slong i = 0; i < nAB; i++)
	{
		fscanf(input,"%ld%ld",&AB[i].r,&AB[i].p);
		lAB[i] = log(AB[i].p);
	}
	fscanf(input,"%lu",&nQB);
	for(slong i = 0; i < nQB; i++)
	{
		fscanf(input,"%ld%ld",&QB[i].r,&QB[i].p);
	}

	MyPair abPairs[2*MaxPrimeBufSize+1];
	ulong num = 2+nRB+nAB+nQB; /*Number of (a,b) pairs to search*/
	ulong N = smoothBound*Afactor;

	if(mynode == 0) fprintf(stdout,"Sieving in region [-%lu,%lu]x[%d,%ld].\n",N,N,1,MaxB);

	sieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs, num, N, m, mynode, totalnodes, &status);

	if(mynode == 0)
	{
		fmpz_fprint(output,n); fprintf(output,"\n");
		fmpz_fprint(output,m); fprintf(output,"\n");
		fmpz_poly_fprint(output,f); fprintf(output,"\n");
		fprintf(output,"%lu\n",nRB);
		printListOfNumbers(output,RB,nRB,10);
		fprintf(output,"%lu\n",nAB);
		printListOfPairs(output,AB,nAB,5);
		fprintf(output,"%lu\n",nQB);
		printListOfPairs(output,QB,nQB,5);
		fprintf(output,"%lu\n",num);
		printListOfPairs(output,abPairs,num,5);

		fmpz_clear(n);
		fmpz_clear(m);
		fmpz_poly_clear(f);

		fclose(output);
	}
	fclose(input);

	MPI_Finalize();
	return 0;
}
