#include "GNFS.h"
#include "GNFS-lattice.h"
#include "HashTable.h"
/**
 *	Sieving: Lattice sieve in the rational part.
 */
void latticeRationalSieve(HashTable &pairs, const ulong *RB, const double *lRB,
	ulong iRB, ulong nRB, MyPair u, MyPair v, slong C, slong D, double logq, const fmpz_t m)
{
	HashTable result(300*C*log(D));
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
					result.add(MyPair(a,b),lp);
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
					result.add(MyPair(a,b),lp);
					a += v.r;
					b += v.p;
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
					result.add(MyPair(a,b),lp);
					a += v.r;
					b += v.p;
				}
			}
		}
	}
	//cout << "Rat accum: " << clock()-start << endl;
	start = clock();
	for(list<MyPair>::iterator iter = result._list.begin(); iter != result._list.end(); iter++)
	{
		fmpz_t bm,fa,fb,abm,gcd;
		fmpz_init(bm);
		fmpz_init(fa);
		fmpz_init(fb);
		fmpz_init(abm);
		fmpz_init(gcd);

		slong a = iter->r;
		slong b = iter->p;
		fmpz_set_si(fa,a);
		fmpz_set_si(fb,b);
		fmpz_mul_si(bm,m,b);
		fmpz_add(abm,bm,fa);
		fmpz_gcd(gcd,fa,fb);
		if(!fmpz_is_one(gcd)) continue;
		fmpz_abs(abm,abm);
		double labm = fmpz_dlog(abm);
		if(result.get(*iter) - labm < -threshold - logq) continue;
		if(!isSmooth(abm,RB,nRB)) continue;
		pairs.insert(*iter);

		fmpz_clear(bm);
		fmpz_clear(fa);
		fmpz_clear(fb);
		fmpz_clear(abm);
		fmpz_clear(gcd);
	}
	//cout << "Rat check: " << clock()-start << endl;
}

/**
 *	Sieving: Lattice sieve in the algebraic part.
 */
void latticeAlgebraicSieve(HashTable &abPairs, ulong &loc, slong num, HashTable &pairs, const fmpz_poly_t f,
	const MyPair *AB, const double *lAB, ulong iAB, ulong nAB, MyPair u, MyPair v, slong C, slong D)
{
	HashTable result(300*C*log(D));
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
					if(!pairs.find(MyPair(a,b))) continue;
					result.add(MyPair(a,b),lp);
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
					if(!pairs.find(MyPair(a,b))) continue;
					result.add(MyPair(a,b),lp);
					a += v.r;
					b += v.p;
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
					if(!pairs.find(MyPair(a,b))) continue;
					result.add(MyPair(a,b),lp);
					a += v.r;
					b += v.p;
				}
			}
		}
	}
	// cout << "Alg accum: " << clock() - start << endl;
	start = clock();
	for(list<MyPair>::iterator iter = result._list.begin(); iter != result._list.end(); iter++)
	{
		fmpz_t fa,fb,nm,gcd;
		fmpz_init(fa);
		fmpz_init(fb);
		fmpz_init(nm);
		fmpz_init(gcd);

		slong a = iter->r;
		slong b = iter->p;
		fmpz_set_si(fa,a);
		fmpz_set_si(fb,b);
		norm(nm,f,fa,fb);
		fmpz_gcd(gcd,fa,fb);
		if(!fmpz_is_one(gcd)) continue;
		fmpz_abs(nm,nm);

		double lnm = fmpz_dlog(nm);
		if(result.get(*iter) - lnm < -threshold) continue;
		if(!isSmooth(nm,AB,nAB)) continue;
		if(pairs.find(*iter))
		{
			abPairs.insert(*iter);
			loc = abPairs.size();
			if(loc >= num) break;
		}

		fmpz_clear(fa);
		fmpz_clear(fb);
		fmpz_clear(nm);
		fmpz_clear(gcd);
	}
	// cout << "Alg check: " << clock() - start << endl;
}

/**
 *	Sieving: The main procedure.
 */
void latticeSieve(const fmpz_poly_t f, const ulong *RB, const double *lRB, ulong nRB,
		   const MyPair *AB, const double *lAB, ulong nAB, MyPair *abPairs, ulong num, slong A, slong B, fmpz_t m)
{
	ulong loc = 0;
	/* Use the HashTable to store the found (a,b) pairs to prevent repeat. */
	HashTable abpairs(num);
	/* Loop for each special-q. */
	for(ulong i = nRB/smoothfactor; i < nRB; i++)
	{
		HashTable pairs(A);

		slong q = RB[i];
		double logq = log(q);
		slong im = fmpz_get_si(m);
		MyPair u(q,0), v(im,-1);
		gaussianLatticeReduce(u,v);
		slong C1,C2,D1,D2;
		getBound(C1,C2,D1,D2,A,B,u,v);
		slong C = C2;
		slong D = D2;
		latticeRationalSieve(pairs, RB, lRB, i, nRB, u, v, C, D, logq, m);
		latticeAlgebraicSieve(abpairs, loc, num, pairs, f, AB, lAB, i, nAB, u, v, C, D);

#if(PRINT_PROCESS && PRINT_SIEVE_PROCESS)
		cerr << "\r" << loc << "/" << num << "       ";
		cerr << i << "/" << nRB;
		cerr << "                       "; cerr.flush();
#endif
#if(PRINT_PROCESS && SLOW_PRINT_SIEVE_PROCESS)
		cerr << loc << "/" << num << "       ";
		cerr << i << "/" << nRB;
		cerr << endl;
#endif
		if(loc >= num) break;
	}
	if(loc < num) /*If loc < num, then continue sieving with traditional method*/
	{
		A *= 2; /*Double the size of the sieving region.*/
		fmpz_t bm,fa,fb,nm,abm,gcd;
		fmpz_init(bm);
		fmpz_init(fa);
		fmpz_init(fb);
		fmpz_init(nm);
		fmpz_init(abm);
		fmpz_init(gcd);
		ulong I = 2*A+1;
		double *r_sieve_array = new double[I];
		double *a_sieve_array = new double[I];

		for(ulong b = 1; b <= 4*B; b++)
		{
			fmpz_set_ui(fb,b);
			fmpz_mul_ui(bm,m,b);
			rationalSieve(r_sieve_array,I,RB,lRB,nRB,-A,bm);
			algebraicSieve(a_sieve_array,f,AB,lAB,nAB,I,-A,b);
			for(slong i = 0; i < I; i++)
			{
				slong a = i - A;
				fmpz_set_si(fa,a);
				fmpz_gcd(gcd,fa,fb);
				norm(nm,f,fa,fb);
				fmpz_mul(abm,fb,m);
				fmpz_add(abm,abm,fa);
				fmpz_abs(abm,abm);
				fmpz_abs(nm,nm);
				if(r_sieve_array[i] >= -5.0 && a_sieve_array[i] >= -5.0 && fmpz_is_one(gcd))
				{
					if(abpairs.find(MyPair(a,b))) continue;
					if(isSmooth(abm,RB,nRB) && isSmooth(nm,AB,nAB))
					{
						abpairs.insert(MyPair(a,b));
						loc = abpairs.size();
						if(loc >= num) break;
					}
				}
#if(PRINT_PROCESS && PRINT_SIEVE_PROCESS)
				cerr << "\r" << loc << "/" << num; cerr.flush();
#endif
			}
#if(PRINT_PROCESS && SLOW_PRINT_SIEVE_PROCESS)
			cerr << loc << "/" << num << endl;
#endif
			if(loc >= num) break;
		}
		delete []r_sieve_array;
		delete []a_sieve_array;
		fmpz_clear(bm);
		fmpz_clear(fa);
		fmpz_clear(fb);
		fmpz_clear(nm);
		fmpz_clear(abm);
		fmpz_clear(gcd);
	}
	slong k = 0;
	for(list<MyPair>::iterator iter = abpairs._list.begin(); iter!=abpairs._list.end(); iter++)
	{
		abPairs[k++] = *iter;
	}
	assert(k==num);
	assert(loc == num);
	num = loc;
}

