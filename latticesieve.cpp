#include <mpi.h>
#include <algorithm>
#include <unistd.h>
#include "GNFS.h"
#include "util.h"
#include "poly.h"
#include "latticesieve.h"
#include "latticeutil.h"

double abratio = 4.0;
int Afactor = DefaultAfactor;
int MaxPrime = DefaultMaxPrime;
double smoothfactor = 5;
double threshold = DefaultThreshold;

double **createCDTable(slong C, slong D)
{
	double **table = new double*[2*C+1];
	for(slong c = 0; c <= 2*C; c++)
		table[c] = new double[D+1];
	return table;
}

int **createMarkTable(slong C, slong D)
{
	int **table = new int*[2*C+1];
	for(slong c = 0; c <= 2*C; c++)
		table[c] = new int[D+1];
	return table;
}

void freeCDTable(double **table, slong C)
{
	for(slong c = 0; c <= 2*C; c++)
		delete []table[c];
	delete []table;
}

void freeMarkTable(int **table, slong C)
{
	for(slong c = 0; c <= 2*C; c++)
		delete []table[c];
	delete []table;
}

/**
 *	Sieving: Lattice sieve in the rational part.
 */
void latticeRationalSieve(double **cdTable, int **marks, const ulong *RB, const double *lRB,
	ulong iRB, ulong nRB, MyPair u, MyPair v, slong C, slong D, double logq, const fmpz_t m)
{
	for(slong c = -C; c <= C; c++)
	{
		for(slong d = 0; d <= D; d++)
		{
			slong a = c*u.r+d*v.r;
			slong b = c*u.p+d*v.p;
			fmpz_t bm,fa,fb,abm;
			fmpz_init(bm);
			fmpz_init_set_si(fa,a);
			fmpz_init_set_si(fb,b);
			fmpz_init(abm);
			fmpz_mul_si(bm,m,b);
			fmpz_add(abm,bm,fa);
			fmpz_abs(abm,abm);
			double labm = fmpz_dlog(abm);
			cdTable[c+C][d] = -labm;
			marks[c+C][d] = 0;
			fmpz_clear(bm);
			fmpz_clear(fa);
			fmpz_clear(fb);
			fmpz_clear(abm);
		}
	}
	slong im = fmpz_get_si(m);
	time_t start = clock();
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
#if(DEBUG)
				assert((c*h+d*k)%p==0);
#endif
				if(d > 1) d -= ((d-1)/p)*p;
				if(d < 1) d += ((-d)/p+1)*p;
#if(DEBUG)
				assert(d >= 1);
#endif
				for(;d <= D; d+=p)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#if(DEBUG)
					assert((a+b*im)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					cdTable[c+C][d] += lp;
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
				slong d = 1;
				for(; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#if(DEBUG)
					assert((a+b*im)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					cdTable[c+C][d] += lp;
				}
			}
		}
		else
		{
			for(slong c = -C; c <= C; c++)
			{
				slong d = 1;
				for(; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#if(DEBUG)
					assert((a+b*im)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					cdTable[c+C][d] += lp;
				}
			}
		}
	}
	start = clock();
	for(slong c = -C; c <= C; c++)
	{
		for(slong d = 0; d <= D; d++)
		{
			if(cdTable[c+C][d] >= -threshold - logq)
				marks[c+C][d] = 1;
		}
	}
}

/**
 *	Sieving: Lattice sieve in the algebraic part.
 */
void latticeAlgebraicSieve(double **cdTable, int **marks, ulong &loc, slong num, const fmpz_poly_t f,
	const MyPair *AB, const double *lAB, ulong iAB, ulong nAB, MyPair u, MyPair v, slong C, slong D)
{
	for(slong c = -C; c <= C; c++)
	{
		for(slong d = 0; d <= D; d++)
		{
			slong a = c*u.r+d*v.r;
			slong b = c*u.p+d*v.p;
			fmpz_t fa,fb,nm;
			fmpz_init_set_si(fa,a);
			fmpz_init_set_si(fb,b);
			fmpz_init(nm);
			norm(nm,f,fa,fb);
			fmpz_abs(nm,nm);
			double lnm = fmpz_dlog(nm);
			cdTable[c+C][d] = -lnm;
			fmpz_clear(fa);
			fmpz_clear(fb);
			fmpz_clear(nm);
		}
	}
	time_t start = clock();
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
#if(DEBUG)
				assert((c*h+d*k)%p==0);
#endif
				if(d > 1) d -= ((d-1)/p)*p;
				if(d < 1) d += ((-d)/p+1)*p;
#if(DEBUG)
				assert(d >= 1);
#endif
				for(;d <= D; d+=p)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#if(DEBUG)
					assert((a+b*r)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					if(!marks[c+C][d]) continue;
					cdTable[c+C][d] += lp;
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
				slong d = 1;
				for(; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#if(DEBUG)
					assert((a+b*r)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					if(!marks[c+C][d]) continue;
					cdTable[c+C][d] += lp;
				}
			}
		}
		else
		{
			for(slong c = -C; c <= C; c++)
			{
				slong d = 1;
				for(; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#if(DEBUG)
					assert((a+b*r)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					if(!marks[c+C][d]) continue;
					cdTable[c+C][d] += lp;
				}
			}
		}
	}
	start = clock();
	for(slong c = -C; c <= C; c++)
	{
		for(slong d = 0; d <= D; d++)
		{
			if(!marks[c+C][d]) continue;
			if(cdTable[c+C][d] < -threshold) continue;
			marks[c+C][d] = 2;
		}
	}
}

/**
 *	Sieving: The main procedure.
 */
void latticeSieve(const fmpz_poly_t f, const ulong *RB, const double *lRB, ulong nRB,
		   const MyPair *AB, const double *lAB, ulong nAB, MyPair *abPairs, const ulong num,
		   slong A, slong B, fmpz_t m, ulong &start, ulong &found,
		   int mynode, int totalnodes, MPI_Status *status)
{
	int end = 0, flag = 0;
	ulong loc = found;
	if(mynode) loc = 0;
	/* Loop for each special-q. */
	for(ulong i = start+mynode; i < nRB; i+=totalnodes)
	{
		slong q = RB[i];
		double logq = log(q);
		slong im = fmpz_get_si(m);
		MyPair u(q,0), v(im,-1);
		gaussianLatticeReduce(u,v);
		slong C1,C2,D1,D2;
		getBound(C1,C2,D1,D2,A,B,u,v);
		slong C = C2;
		slong D = D2;
		double **cdTable = createCDTable(C,D);
		int **marks = createMarkTable(C,D);

		latticeRationalSieve(cdTable, marks, RB, lRB, i, nRB, u, v, C, D, logq, m);
		latticeAlgebraicSieve(cdTable, marks, loc, num, f, AB, lAB, i, nAB, u, v, C, D);
		if(mynode) loc = 0;
		for(slong c = -C; c <= C; c++)
		{
			for(slong d = 0; d <= D; d++)
			{
				if(marks[c+C][d] != 2) continue;
				slong a = c*u.r+d*v.r;
				slong b = c*u.p+d*v.p;
				if(b < 0){a=-a;b=-b;}
				fmpz_t bm,fa,fb,nm,abm,gcd;
				fmpz_init(bm);
				fmpz_init(nm);
				fmpz_init(abm);
				fmpz_init(gcd);
				fmpz_init_set_si(fa,a);
				fmpz_init_set_si(fb,b);
				fmpz_mul_ui(bm,m,b);
				fmpz_add(abm,bm,fa);
				fmpz_gcd(gcd,fa,fb);
				fmpz_abs(gcd,gcd);
				if(!fmpz_is_one(gcd)) continue;
				norm(nm,f,fa,fb);
				fmpz_abs(abm,abm);
				fmpz_abs(nm,nm);
				if(isSmooth(abm,RB,nRB) && isSmooth(nm,AB,nAB))
				{
					abPairs[loc] = MyPair(a,b);
					loc++;
					if(loc >= num) break;
				}

#if(PRINT_PROCESS && PRINT_SIEVE_PROCESS )
				if(mynode == 0)
				{
					std::cerr << "\r" << loc << "/" << num << "       ";
					std::cerr << i << "/" << nRB;
					std::cerr << "                       "; std::cerr.flush();
				}
#endif
				if(mynode)
				{
					MPI_Iprobe(0,1,MPI_COMM_WORLD,&flag,status);
					if(flag) break;
				}
			}
			if(loc >= num) break;
		}
		if(mynode == 0) start = i+totalnodes;

		freeCDTable(cdTable,C);
		freeMarkTable(marks,C);

		if(mynode == 0)
		{
#if(PRINT_PROCESS && SLOW_PRINT_SIEVE_PROCESS)
			std::cerr << loc << "/" << num << "       ";
			std::cerr << i << "/" << nRB;
			std::cerr << std::endl;
#endif
			MyPair *recvbuf = new MyPair[num];
			for(int k = 1; k < totalnodes; k++)
			{
				int recvsize;
				MPI_Recv(&recvsize,1,MPI_UNSIGNED_LONG,k,0,MPI_COMM_WORLD,status);
				MPI_Recv(recvbuf,num*2,MPI_LONG,k,0,MPI_COMM_WORLD,status)/2;
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
			}
			delete []recvbuf;
			if(loc >= num) break;
		}
		else
		{
			MPI_Iprobe(0,1,MPI_COMM_WORLD,&flag,status);
			if(flag) break;
#if(DEBUG)
			std::cerr << "Found " << loc << " pairs in q = " << q << "." << std::endl;
#endif
			MPI_Send(&loc,1,MPI_UNSIGNED_LONG,0,0,MPI_COMM_WORLD);
			MPI_Send(abPairs,loc*2,MPI_LONG,0,0,MPI_COMM_WORLD);

			MPI_Iprobe(0,1,MPI_COMM_WORLD,&flag,status);
			if(flag) break;
		}
	}
	if(mynode == 0)
	{
#if(PRINT_PROCESS && PRINT_SIEVE_PROCESS )
		std::cerr << std::endl;
#endif
		found = loc;
		assert(found >= num);
		end = 1;
		for(int k = 1; k < totalnodes; k++)
			MPI_Send(&end,1,MPI_INT,k,1,MPI_COMM_WORLD);
	}
	else
	{
#if(DEBUG)
		std::cerr << "Node " << mynode << ": Received end message" << std::endl;
#endif
		MPI_Recv(&end,1,MPI_INT,0,1,MPI_COMM_WORLD,status);
	}
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
	while((ch = getopt(argc-2,argv+2,"a:b:s:")) != -1)
	{
		switch(ch)
		{
		case 'a': Afactor = atoi(optarg); break;
		case 'b': abratio = atof(optarg); break;
		case 's': smoothfactor = atof(optarg); break;
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
	slong A = smoothBound*Afactor, B = A/abratio;
	fprintf(stdout,"Sieving in region [%ld,%ld]x[%d,%ld].\n",-A,A,1,B);
	ulong start = nRB/smoothfactor, found = 0;
	while(true)
	{
		latticeSieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs, num, A, B, m, start,
			found, mynode, totalnodes, &status);
		if(mynode == 0)
		{
			sort(&abPairs[0],&abPairs[num]);
			found = unique(&abPairs[0],&abPairs[found])-&abPairs[0];
			int enough = 1;
			if(found >= num)
			{
				for(int k = 1; k < totalnodes; k++)
					MPI_Send(&enough,1,MPI_INT,k,0,MPI_COMM_WORLD);
				break;
			}
			else
			{
				enough = 0;
				for(int k = 1; k < totalnodes; k++)
				{
					MPI_Send(&enough,1,MPI_INT,k,0,MPI_COMM_WORLD);
					MPI_Send(&start,1,MPI_UNSIGNED_LONG,k,0,MPI_COMM_WORLD);
				}
			}
		}
		else
		{
			int enough = 1;
			MPI_Recv(&enough,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
			if(enough) break;
			MPI_Recv(&start,1,MPI_UNSIGNED_LONG,0,0,MPI_COMM_WORLD,&status);
		}
	}

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
