#include "GNFS.h"
#include "GNFS-lattice.h"

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
#ifdef DEBUG
				assert((c*h+d*k)%p==0);
#endif
				if(d > 1) d -= ((d-1)/p)*p;
				if(d < 1) d += ((-d)/p+1)*p;
#ifdef DEBUG
				assert(d >= 1);
#endif
				for(;d <= D; d+=p)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#ifdef DEBUG
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
#ifdef DEBUG
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
#ifdef DEBUG
					assert((a+b*im)%p == 0);
#endif
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					cdTable[c+C][d] += lp;
				}
			}
		}
	}
	//cout << "Rat accum: " << clock()-start << endl;
	start = clock();
	for(slong c = -C; c <= C; c++)
	{
		for(slong d = 0; d <= D; d++)
		{
			if(cdTable[c+C][d] >= -threshold - logq)
				marks[c+C][d] = 1;
		}
	}
	//cout << "Rat check: " << clock()-start << endl;
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
#ifdef DEBUG
				assert((c*h+d*k)%p==0);
#endif
				if(d > 1) d -= ((d-1)/p)*p;
				if(d < 1) d += ((-d)/p+1)*p;
#ifdef DEBUG
				assert(d >= 1);
#endif
				for(;d <= D; d+=p)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
#ifdef DEBUG
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
#ifdef DEBUG
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
#ifdef DEBUG
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
	// cout << "Alg accum: " << clock() - start << endl;
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
		   const MyPair *AB, const double *lAB, ulong nAB, MyPair *abPairs, ulong num, slong A, slong B, fmpz_t m)
{
	ulong loc = 0;
	/* Loop for each special-q. */
	for(ulong i = nRB/smoothfactor; i < nRB; i++)
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
				norm(nm,f,fa,fb);
				fmpz_abs(abm,abm);
				fmpz_abs(nm,nm);
				if(isSmooth(abm,RB,nRB) && isSmooth(nm,AB,nAB))
				{
					abPairs[loc] = MyPair(a,b);
					loc++;
					if(loc >= num) break;
				}

#if(PRINT_PROCESS && PRINT_SIEVE_PROCESS)
				cerr << "\r" << loc << "/" << num << "       ";
				cerr << i << "/" << nRB;
				cerr << "                       "; cerr.flush();
#endif
			}
			if(loc >= num) break;
		}

		freeCDTable(cdTable,C);
		freeMarkTable(marks,C);
#if(PRINT_PROCESS && SLOW_PRINT_SIEVE_PROCESS)
		cerr << loc << "/" << num << "       ";
		cerr << i << "/" << nRB;
		cerr << endl;
#endif
		if(loc >= num) break;
	}
	assert(loc == num);
	num = loc;
}

