#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <list>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <limits.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>

#define DEBUG 1
#define ONLY_RESULT 0

#if(DEBUG)
#define PRINT_MDF
//#define PRINT_SMOOTH_BOUND
#define PRINT_PROCESS
//#define PRINT_RATIONAL_BASE
//#define PRINT_ALGEBRAIC_BASE
//#define PRINT_QUADRATIC_BASE
//#define PRINT_SELECTED_ABPAIRS
//#define PRINT_MATRIX
//#define PRINT_SELECTED_SQUARE_ABPAIRS
//#define PRINT_PROD_ABM
//#define PRINT_UPPDER_BOUND
//#define PRINT_PRIMES
//#define PRINT_XI
#endif

using namespace std;

const int MaxPrimeBufSize = 30000;
const int MaxPrime = 20000000;
const int MaxSelectedPrimes = 10000;
const int MaxB = 10240;

typedef struct MyPair
{
	slong r;
	slong p;
	MyPair(slong a, slong b){r=a; p=b;}
	MyPair(){r=0;p=0;}
} MyPair;

typedef struct PairDouble
{
	MyPair p;
	double v;
	PairDouble(MyPair a, double b){p=a;v=b;}
	PairDouble():p(){v=0.0;}
} PairDouble;

bool operator==(const MyPair &a, const MyPair &b)
{
	return a.r == b.r && a.p == b.p;
}

class HashTable{
public:
	ulong _T;
	HashTable(ulong T):_T(T){_data = new list<PairDouble>[T];}
	~HashTable(){
		if(_list.size()/_T >= 3)cerr << endl << "Warning: Load factor = " << (double)_list.size()/_T << endl;
		delete []_data;
	}
	list<PairDouble> * _data;
	list<MyPair> _list;
	slong hash(const MyPair &p)
	{
		return (p.r*p.r+p.p+(slong)sqrt(p.r*p.r+p.p*p.p)+p.r*(_T/2))%_T;
	}
	void set(const MyPair &p, double v)
	{
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
			if(iter->p==p) break;
		if(iter!=_data[s].end()) iter->v=v;
		else
		{
			_data[s].push_back(PairDouble(p,v));
			_list.push_back(p);
		}
	}
	void add(const MyPair &p, double v)
	{
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
			if(iter->p==p) break;
		if(iter!=_data[s].end()) iter->v+=v;
		else
		{
			_data[s].push_back(PairDouble(p,v));
			_list.push_back(p);
		}
	}
	double get(const MyPair &p)
	{
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
			if(iter->p==p) break;
		if(iter!=_data[s].end()) return iter->v;
		else return 0.0;
	}
	double find(const MyPair &p)
	{
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
			if(iter->p==p) break;
		return iter!=_data[s].end();
	}
};

/*
 * Given integer n, return f, m, d
 * f: Integer array representing the coefficients of polynomial
 * m: the number selected
 * d: degree of the polynomial
 */
void selectPolynomial(const fmpz_t n, fmpz_poly_t f, fmpz_t m, ulong &d)
{
	fmpz_t lead,tail,N,Nmm;
	assert(fmpz_cmp_ui(n,0) > 0);
	assert(fmpz_dlog(n) > 0);
	fmpz_poly_init(f);
	fmpz_init(m);
	fmpz_init(N);
	fmpz_init(lead);
	fmpz_init(tail);
	fmpz_init(Nmm);
	d = pow(3*fmpz_dlog(n)/log(fmpz_dlog(n)),0.333);
	if(d % 2 == 0) d++;
	while(true)
	{
		fmpz_root(m,n,d);
		if(fmpz_cmp_ui(m,2000)>0)
		{
			do
			{
				ulong l;
				fmpz_t r,k;
				if(fmpz_cmp_ui(m,2000)>0)
					fmpz_fdiv_q_ui(r,m,1000);
				else
					fmpz_fdiv_q_ui(r,m,10);
				if((l=fmpz_get_ui(r)) > RAND_MAX)
					fmpz_sub_ui(m,m,rand());
				else
				{
					fmpz_sub_ui(m,m,rand()%l);
				}
			} while(fmpz_divisible(n,m));
		}

		fmpz_t mtod;
		fmpz_init(mtod);
		fmpz_pow_ui(mtod,m,d);
		assert(fmpz_cmp(mtod,n) < 0);
		fmpz_clear(mtod);

		fmpz_set(N,n);
		for(int i = 0; i <= d; i++)
		{
			fmpz_mod(Nmm,N,m);
			fmpz_poly_set_coeff_fmpz(f,i,Nmm);
			fmpz_fdiv_q(N,N,m);
		}
		fmpz_poly_get_coeff_fmpz(lead,f,d);
		fmpz_poly_get_coeff_fmpz(tail,f,0);
		if(fmpz_is_one(lead) && !fmpz_is_zero(tail) && !fmpz_is_one(tail)) break;
	}
	fmpz_clear(lead);
	fmpz_clear(tail);
	fmpz_clear(N);
	fmpz_clear(Nmm);
}

slong length(MyPair u)
{
	return u.p*u.p+u.r*u.r;
}

void gaussianLatticeReduce(MyPair &u, MyPair &v)
{
	slong lu = length(u), lv = length(v);
	do{
		if(lu > lv)
		{
			MyPair tmp = u;
			u = v;
			v = tmp;
			lu = length(u);
			lv = length(v);
		}
		if(lu == 0) return;
		slong k = (u.p*v.p+u.r*v.r)/lu;
		slong k1 = k + 1;
		slong l = (v.p-k*u.p)*(v.p-k*u.p)+(v.r-k*u.r)*(v.r-k*u.r);
		slong l1 = (v.p-k1*u.p)*(v.p-k1*u.p)+(v.r-k1*u.r)*(v.r-k1*u.r);
		if(l > l1)
		{
			v.p -= k1*u.p;
			v.r -= k1*u.r;
		}
		else
		{
			v.p -= k*u.p;
			v.r -= k*u.r;
		}
		lv = length(v);
	} while(lu > lv);
}

bool randombool(double p)
{
	return rand() < RAND_MAX * p;
}

void norm(fmpz_t nm, const fmpz_poly_t f, const fmpz_t a, const fmpz_t b)
{
	fmpz_t poa,pob,mb,c,ab;
	fmpz_init(poa);
	fmpz_init(pob);
	fmpz_init(mb);
	fmpz_init(c);
	fmpz_init(ab);
	fmpz_one(poa);
	fmpz_one(pob);
	fmpz_zero(nm);
	fmpz_neg(mb,b);
	ulong d = fmpz_poly_degree(f);
	if(fmpz_is_zero(a))
	{
		for(ulong i = 0; i < d; i++)
			fmpz_mul(pob,pob,mb);
		fmpz_poly_get_coeff_fmpz(c,f,0);
		fmpz_mul(nm,pob,c);
	}
	else
	{
		for(ulong i = 0; i < d; i++)
			fmpz_mul(poa,poa,a);
		for(slong i = d; i >= 0; i--)
		{
			fmpz_poly_get_coeff_fmpz(c,f,i);
			fmpz_mul(ab,poa,pob);
			fmpz_mul(ab,ab,c);
			fmpz_add(nm,nm,ab);
			fmpz_mul(pob,pob,mb);
			fmpz_fdiv_q(poa,poa,a);
		}
	}

	fmpz_clear(poa);
	fmpz_clear(pob);
	fmpz_clear(mb);
	fmpz_clear(c);
	fmpz_clear(ab);
}

bool isSmooth(const fmpz_t num, const ulong *ps, ulong np)
{
	fmpz_t k;
	fmpz_init_set(k,num);
	for(ulong i = 0; i < np; i++)
	{
		ulong p = ps[i];
		while(fmpz_cmp_ui(k,0)!=0 && fmpz_divisible_si(k,p))
		{
			fmpz_fdiv_q_ui(k,k,p);
		}
	}
	fmpz_clear(k);
	return fmpz_is_pm1(k);
}

bool isSmooth(const fmpz_t num, const MyPair *ps, ulong np)
{
	fmpz_t k;
	fmpz_init_set(k,num);
	for(ulong i = 0; i < np; i++)
	{
		ulong p = ps[i].p;
		while(fmpz_cmp_ui(k,0)!=0 && fmpz_divisible_si(k,p))
			fmpz_fdiv_q_ui(k,k,p);
	}
	fmpz_clear(k);
	return fmpz_is_pm1(k);
}

void getBound(slong &E1, slong &E2, slong &F1, slong &F2, slong C, slong D, MyPair s, MyPair t)
{
	if(s.r == 0) E1 = D/abs(s.p);
	else if(s.p == 0) E1 = C/abs(s.r);
	else
	{
		E1 = C/abs(s.r);
		if(D/abs(s.p) < E1) E1 = D/abs(s.p);
	}
	E2 = E1 * 1.2;
	E1 = -E1;
	if(t.r == 0) F2 = D/abs(t.p);
	else if(t.p == 0) F2 = C/abs(t.r);
	else
	{
		F2 = C/abs(t.r);
		if(D/abs(t.p) < F2) F2 = D/abs(t.p);
	}
	F2 = F2 * 1.2;
	F1 = 1;
	assert(E2 > 0 && F2 > 0);
}

void latticeRationalSieve(HashTable &pairs, const ulong *RB, const double *lRB,
	ulong iRB, ulong nRB, MyPair u, MyPair v, slong C, slong D, const fmpz_t m)
{
	HashTable result(800*C);
	slong im = fmpz_get_si(m);
	MyPair s,t;
	for(ulong j = 0; j < iRB; j++)
	{
		slong p = RB[j];
		slong h = (u.r+u.p*im)%p;
		slong k = (v.r+v.p*im)%p;
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
				assert((c*h+d*k)%p==0);
				if(d > 1) d -= ((d-1)/p)*p;
				if(d < 1) d += ((-d)/p+1)*p;
				assert(d >= 1);
				for(;d <= D; d+=p)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
					assert((a+b*im)%p == 0);
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					result.add(MyPair(a,b),lRB[j]);
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
				for(slong d = 1; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
					assert((a+b*im)%p == 0);
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					result.add(MyPair(a,b),lRB[j]);
				}
			}
		}
		else
		{
			for(slong c = -C; c <= C; c++)
				for(slong d = 1; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
					assert((a+b*im)%p == 0);
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					result.add(MyPair(a,b),lRB[j]);
				}
		}
	}
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
		if(result.get(*iter) - labm < -10.0) continue;
		if(!isSmooth(abm,RB,nRB)) continue;
		pairs.set(*iter,0.0);

		fmpz_clear(bm);
		fmpz_clear(fa);
		fmpz_clear(fb);
		fmpz_clear(abm);
		fmpz_clear(gcd);
	}
}

void latticeAlgebraicSieve(HashTable &abPairs, ulong &loc, slong num, HashTable &pairs, const fmpz_poly_t f,
	const MyPair *AB, const double *lAB, ulong iAB, ulong nAB, MyPair u, MyPair v, slong C, slong D)
{
	HashTable result(800*C);
	MyPair s,t;
	for(ulong j = 0; j < nAB; j++)
	{
		slong p = AB[j].p;
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
				assert((c*h+d*k)%p==0);
				if(d > 1) d -= ((d-1)/p)*p;
				if(d < 1) d += ((-d)/p+1)*p;
				assert(d >= 1);
				for(;d <= D; d+=p)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
					assert((a+b*r)%p == 0);
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					result.add(MyPair(a,b),lAB[j]);
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
				for(slong d = 1; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
					assert((a+b*r)%p == 0);
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					result.add(MyPair(a,b),lAB[j]);
				}
			}
		}
		else
		{
			for(slong c = -C; c <= C; c++)
				for(slong d = 1; d <= D; d++)
				{
					slong a = c*u.r+d*v.r;
					slong b = c*u.p+d*v.p;
					assert((a+b*r)%p == 0);
					if(b==0) continue;
					if(b < 0){a=-a;b=-b;}
					if(!pairs.find(MyPair(a,b))) continue;
					result.add(MyPair(a,b),lAB[j]);
				}
		}
	}
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
		if(result.get(*iter) - lnm < -10.0) continue;
		if(!isSmooth(nm,AB,nAB)) continue;
		if(pairs.find(*iter))
		{
			abPairs.set(*iter,0.0);
			loc = abPairs._list.size();
			if(loc >= num) break;
		}

		fmpz_clear(fa);
		fmpz_clear(fb);
		fmpz_clear(nm);
		fmpz_clear(gcd);
	}
}
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
			slong a = i - N;
			fmpz_set_si(fa,a);
			fmpz_add(abm,bm,fa);
			fmpz_gcd(gcd,fa,fb);
			norm(nm,f,fa,fb);
			fmpz_abs(abm,abm);
			fmpz_abs(nm,nm);
			if(r_sieve_array[i] >= -5.0 && a_sieve_array[i] >= -5.0 && fmpz_is_one(gcd))
			{
				if(randombool(0.9))
				{
					if(isSmooth(abm,RB,nRB) && isSmooth(nm,AB,nAB))
					{
						abPairs[loc] = MyPair(a,b);
						loc++;
						if(loc >= num) break;
#ifdef PRINT_PROCESS
						cerr << "\r" << loc*100/num << "%"; cerr.flush();
#endif
					}
				}
			}
		}
		if(loc >= num) break;
	}
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

void latticeSieve(const fmpz_poly_t f, const ulong *RB, const double *lRB, ulong nRB,
		   const MyPair *AB, const double *lAB, ulong nAB, MyPair *abPairs, ulong num, slong A, slong B, fmpz_t m)
{
	ulong loc = 0;
	HashTable abpairs(num);

	for(ulong i = nRB/2; i < nRB; i++)
	{
		HashTable pairs(A*20);

		slong q = RB[i];
		double lq = log(q);
		slong im = fmpz_get_si(m);
		MyPair u(q,0), v(im,-1);
		gaussianLatticeReduce(u,v);
		slong C1,C2,D1,D2;
		getBound(C1,C2,D1,D2,A,B,u,v);
		slong C = C2/2;
		slong D = D2/2;
		latticeRationalSieve(pairs, RB, lRB, i, nRB, u, v, C, D, m);
		latticeAlgebraicSieve(abpairs, loc, num, pairs, f, AB, lAB, i, nAB, u, v, C, D);
		
		if(loc >= num) break;
#ifdef PRINT_PROCESS
		cerr << "\r" << loc << "/" << num << "       ";
		cerr << i << "/" << nRB;
		cerr << "                       "; cerr.flush();
#endif
	}
	if(loc < num) /*If loc < num, then continue sieving with traditional method*/
	{
		A *= 2;
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
						abpairs.set(MyPair(a,b),0.0);
						loc = abpairs._list.size();
						if(loc >= num) break;
					}
				}
#ifdef PRINT_PROCESS
				cerr << "\r" << loc << "/" << num; cerr.flush();
#endif
			}
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

void rootsMod(const fmpz_poly_t f, ulong p, ulong *roots, ulong &nroot)
{
	assert(p > 1);
	ulong d = fmpz_poly_degree(f);
	fmpz_t c;
	nmod_poly_factor_t fac;
	nmod_poly_t g;
	fmpz_init(c);
	nmod_poly_factor_init(fac);
	nmod_poly_factor_fit_length(fac,d);
	nmod_poly_init(g,p);
	for(int i = 0; i <= d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,f,i);
		nmod_poly_set_coeff_ui(g,i,fmpz_mod_ui(c,c,p));
	}
	nmod_poly_factor(fac,g);

	ulong n = fac->num;
	nroot = 0;
	for(long i = 0; i < n; i++)
	{
		nmod_poly_struct* h = fac->p+i;
		if(nmod_poly_degree(h) == 1)
		{
			ulong a = nmod_poly_get_coeff_ui(h,0);
			roots[nroot] = (p-a)%p;
			nroot++;
		}
	}
	fmpz_clear(c);
	nmod_poly_clear(g);
	nmod_poly_factor_clear(fac);
}

int Leg(const fmpz_t a, const fmpz_t p) //Legender symbol
{
	fmpz_t c;
	fmpz_init(c);
	fmpz_mod(c,a,p);
	int l = fmpz_jacobi(c,p);
	fmpz_clear(c);
	return l;
}

int Leg(slong a, ulong p) //Legender symbol
{
	fmpz_t c,q;
	fmpz_init_set_ui(q,p);
	fmpz_init_set_si(c,a);
	fmpz_mod(c,c,q);
	int l = fmpz_jacobi(c,q);
	fmpz_clear(c);
	fmpz_clear(q);
	return l;
}

void solveMatrix(nmod_mat_t mat, ulong I, ulong J, int *vec)
{
	nmod_mat_rref(mat);
	ulong minI = 0;
	for(minI = 0; minI < I; minI++)
	{
		mp_limb_t b = nmod_mat_entry(mat,minI,minI);
		if(b == 0) break;
	}
	vec[minI] = 1;
	for(ulong i = minI+1; i < J; i++)
		vec[i] = 0;
	for(ulong i = 0; i < minI; i++)
		vec[i] = nmod_mat_entry(mat,i,minI);
}

//void nmod_poly_powmod_ui_binexp(res,poly,e,f)

bool legal(const fmpz_poly_t fx, ulong p)
{
	ulong d = fmpz_poly_degree(fx);
	fmpz_t c;
	nmod_poly_t f;
	fmpz_init(c);
	nmod_poly_init(f,p);
	for(slong i = 0; i <= d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,fx,i);
		nmod_poly_set_coeff_ui(f,i,fmpz_mod_ui(c,c,p));
	}
	int r = nmod_poly_is_irreducible(f);
	fmpz_clear(c);
	nmod_poly_clear(f);
	return r;
}

void selectNonResidual(nmod_poly_t px, const nmod_poly_t fx, ulong p, const fmpz_t e, ulong d)
{
	mpz_t ee;
	mpz_init(ee);
	fmpz_get_mpz(ee,e);
	flint_rand_t frt;
	flint_randinit(frt);
	nmod_poly_t mp;
	nmod_poly_init(mp,p);
	while(true)
	{
		nmod_poly_randtest_monic(px,frt,d);
		nmod_poly_powmod_mpz_binexp(mp,px,ee,fx);
		assert(nmod_poly_degree(mp)==0);
		if(nmod_poly_get_coeff_ui(mp,0) == p-1) break;
	}
	flint_randclear(frt);
	nmod_poly_clear(mp);
	mpz_clear(ee);
}

ulong computeOrder(const nmod_poly_t px, const nmod_poly_t fx)
{
	mp_limb_t p = nmod_poly_modulus(px);
	nmod_poly_t g;
	nmod_poly_init(g,p);
	nmod_poly_set(g,px);
	ulong i;
	for(i = 0; i < 10000; i++)
	{
		if(nmod_poly_is_one(g)) break;
		nmod_poly_powmod_ui_binexp(g,g,2,fx);
	}
	nmod_poly_clear(g);
	if(i >= 10000) return -1;
	return i;
}

ulong computeSquareRoot(const fmpz_poly_t sx, const fmpz_poly_t f, slong p, const fmpz_t m, const fmpz_t Nm)
{
	assert(p > 1);

	slong d = fmpz_poly_degree(f);
	fmpz_t c,fp,q,q1,q2,s,s1,Nmp;
	mpz_t eq1,eq2,es,es1;
	mpz_init(eq1);
	mpz_init(eq2);
	mpz_init(es);
	mpz_init(es1);
	nmod_poly_t Fx, Sx, Sxp,lambda,omega,zeta,eta,pzeta1,pzeta2,omega2,NN;
	fmpz_init(c);
	fmpz_init_set_ui(fp,p);
	fmpz_init(q);
	fmpz_init(q1);
	fmpz_init(q2);
	fmpz_init(s);
	fmpz_init(s1);
	fmpz_init(Nmp);
	nmod_poly_init(Fx,p);
	nmod_poly_init(Sx,p);
	nmod_poly_init(Sxp,p);
	nmod_poly_init(lambda,p);
	nmod_poly_init(omega,p);
	nmod_poly_init(zeta,p);
	nmod_poly_init(eta,p);
	nmod_poly_init(pzeta1,p);
	nmod_poly_init(pzeta2,p);
	nmod_poly_init(omega2,p);
	nmod_poly_init(NN,p);

	for(int i = 0; i <= d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,f,i);
		nmod_poly_set_coeff_ui(Fx,i,fmpz_mod_ui(c,c,p));
	}
	for(int i = 0; i < d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,sx,i);
		nmod_poly_set_coeff_ui(Sx,i,fmpz_mod_ui(c,c,p));
	}
	fmpz_pow_ui(q,fp,d);
	fmpz_sub_ui(q,q,1);
	fmpz_fdiv_q_ui(q1,q,2);
	fmpz_fdiv_q_ui(q2,q,p-1);
	fmpz_get_mpz(eq1,q1);
	fmpz_get_mpz(eq2,q2);
	/*Check Sx is quadratic residual*/
	nmod_poly_powmod_mpz_binexp(Sxp,Sx,eq1,Fx);
	assert(nmod_poly_is_one(Sxp));
	/********************************/
	slong r = 1;
	fmpz_set(s,q1);
	while(true)
	{
		if(!fmpz_divisible_si(s,2)) break;
		fmpz_fdiv_q_ui(s,s,2);
		r++;
	}
	fmpz_get_mpz(es,s);
	fmpz_add_ui(s1,s,1);
	fmpz_fdiv_q_ui(s1,s1,2);
	fmpz_get_mpz(es1,s1);
	nmod_poly_powmod_mpz_binexp(lambda,Sx,es,Fx);
	nmod_poly_powmod_mpz_binexp(omega,Sx,es1,Fx);
	selectNonResidual(eta,Fx,p,q1,d);
	nmod_poly_powmod_mpz_binexp(zeta,eta,es,Fx);

	while(!nmod_poly_is_one(lambda))
	{
		slong m = computeOrder(lambda,Fx);
		assert(r > m);
		nmod_poly_set(pzeta1,zeta);
		for(slong i = 0; i < r-m-1; i++)
			nmod_poly_powmod_ui_binexp(pzeta1,pzeta1,2,Fx);
		nmod_poly_powmod_ui_binexp(pzeta2,pzeta1,2,Fx);
		nmod_poly_mulmod(lambda,lambda,pzeta2,Fx);
		nmod_poly_mulmod(omega,omega,pzeta1,Fx);

		nmod_poly_t lambdadelta;
		nmod_poly_init(lambdadelta,p);
		nmod_poly_powmod_ui_binexp(omega2,omega,2,Fx);
		nmod_poly_mulmod(lambdadelta,lambda,Sx,Fx);
		assert(nmod_poly_equal(omega2,lambdadelta));
		nmod_poly_clear(lambdadelta);
	}
	/*Check that omega_n^2 = delta*/
	nmod_poly_powmod_ui_binexp(omega2,omega,2,Fx);
	assert(nmod_poly_equal(omega2,Sx));
	/*If the norm of omega is not equivalent to norm, negate it*/
	nmod_poly_powmod_mpz_binexp(NN,omega,eq2,Fx);
	assert(nmod_poly_degree(NN)==0);
	slong pp = nmod_poly_get_coeff_ui(NN,0);
	slong Nmmp = fmpz_mod_ui(Nmp,Nm,p);
	if(pp != Nmmp)
	{
		nmod_poly_neg(omega,omega);
		nmod_poly_powmod_mpz_binexp(NN,omega,eq2,Fx);
		assert(nmod_poly_degree(NN)==0);
		pp = nmod_poly_get_coeff_ui(NN,0);
		assert(pp == Nmmp);
	}
	/*Substitute m in*/
	slong M = nmod_poly_evaluate_nmod(omega,fmpz_mod_ui(Nmp,m,p));

	mpz_init(eq1);
	mpz_init(eq2);
	mpz_init(es);
	mpz_init(es1);
	fmpz_clear(c);
	fmpz_clear(fp);
	fmpz_clear(q);
	fmpz_clear(q1);
	fmpz_clear(q2);
	fmpz_clear(s);
	fmpz_clear(s1);
	fmpz_clear(Nmp);
	nmod_poly_clear(Fx);
	nmod_poly_clear(Sx);
	nmod_poly_clear(Sxp);
	nmod_poly_clear(lambda);
	nmod_poly_clear(omega);
	nmod_poly_clear(zeta);
	nmod_poly_clear(eta);
	nmod_poly_clear(pzeta1);
	nmod_poly_clear(pzeta2);
	nmod_poly_clear(omega2);
	nmod_poly_clear(NN);
	return M;
}

double doublesum(double array[], ulong l, ulong r)
{
	if(r-l == 0) return 0;
	if(r-l == 1) return array[l];
	return doublesum(array,l,(l+r)/2) + doublesum(array,(l+r)/2,r);
}

ulong boundForSmoothness(ulong d, const fmpz_t n)
{
	assert(d > 0);
	double dlogd = d*log(d);
	assert(n > 0);
	double temp = 1.0/d * fmpz_dlog(n);
	assert(temp > 0);
	double e = dlogd + sqrt(dlogd*dlogd + 4*temp*log(temp));
	return 0.5*exp(0.5*e);
}

void prepareRationalBase(ulong *RB, double *lRB, ulong &nRB, ulong bound)
{
	ulong p;
	n_primes_t iter;
	n_primes_init(iter);
	for(nRB = 0; (p = n_primes_next(iter)) <= bound; nRB++)
	{
		RB[nRB] = p;
		lRB[nRB] = log(RB[nRB]);
	}
	//slong r = rand() % (nRB/10);
	//nRB -= r;
	assert(nRB <= MaxPrimeBufSize);
	n_primes_clear(iter);
}

void prepareAlgebraicBase(MyPair *AB, double *lAB, ulong &nAB, ulong size, fmpz_poly_t f)
{
	n_primes_t iter;
	n_primes_init(iter);
	for(ulong i = 0; i < size; i++)
	{
		ulong p = n_primes_next(iter);
		ulong roots[16]; ulong nroot;
		rootsMod(f,p,roots,nroot);
		for(ulong j = 0; j < nroot; j++)
		{
			AB[nAB] = MyPair(roots[j],p);
			nAB++;
		}
	}
	n_primes_clear(iter);
	assert(nAB <= MaxPrimeBufSize);
	for(ulong i = 0; i < nAB; i++)
		lAB[i] = log(AB[i].p);
}

void prepareQuadraticBase(MyPair *QB, ulong &nQB, ulong min, ulong max, fmpz_poly_t &f)
{
	ulong p;
	n_primes_t iter;
	n_primes_init(iter);
	n_primes_jump_after(iter,min);
	for(ulong i = 0; (p=n_primes_next(iter)) <= max; i++)
	{
		ulong roots[16]; ulong nroot;
		rootsMod(f,p,roots,nroot);
		for(ulong j = 0; j < nroot; j++)
		{
			QB[nQB] = MyPair(roots[j],p);
			nQB++;
		}
	}
	assert(nQB <= MaxPrimeBufSize);
	n_primes_clear(iter);
}

void printListOfNumbers(slong *A, ulong s, ulong N)
{
	for(ulong i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) cout << endl;
		cout << setw(10) << A[i];
	}
	cout << endl;
}

void printListOfNumbers(ulong *A, ulong s, ulong N)
{
	for(ulong i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) cout << endl;
		cout << setw(10) << A[i];
	}
	cout << endl;
}

void printListOfPairs(MyPair *A, ulong s, ulong N)
{
	for(ulong i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) cout << endl;
		cout << "(" << A[i].r << "," << A[i].p << ") ";
	}
	cout << endl;
}

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
		assert(fmpz_is_pm1(A));
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
		assert(fmpz_is_pm1(A));
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
			assert((a+b*s)%q != 0);
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

void select(MyPair *pairs, int *vec, ulong I, ulong &n)
{
	ulong loc = 0;
	for(ulong i = 0; i < I; i++)
	{
		if(vec[i]) pairs[loc++] = pairs[i];
	}
	n = loc;
	assert(n > 0);
}

void sqrtProductOfPairs(fmpz_t s, const MyPair *pairs, ulong num, const fmpz_t m)
{
	fmpz_t abm,bm,fa,r,e;
	fmpz_init(abm);
	fmpz_init(bm);
	fmpz_init(fa);
	fmpz_init(r);
	fmpz_init(e);
	fmpz_one(s);
	for(long i = 0; i < num; i++)
	{
		slong a = pairs[i].r;
		ulong b = pairs[i].p;
		fmpz_set_si(fa,a);
		fmpz_mul_ui(bm,m,b);
		fmpz_add(abm,bm,fa);
		fmpz_mul(s,s,abm);
	}
	fmpz_sqrtrem(r,e,s);
	/*Check that s is indeed a square number.*/
	assert(fmpz_is_zero(e));
	fmpz_set(s,r);

	fmpz_clear(abm);
	fmpz_clear(bm);
	fmpz_clear(fa);
	fmpz_clear(r);
	fmpz_clear(e);
}

void productOfPairs(fmpz_poly_t sx, const MyPair *abPairs, ulong num, const fmpz_poly_t f, fmpz_t Nm)
{
	fmpz_t r,e,nm,fa,fb;
	fmpz_poly_t gx,q;
	fmpz_init(r);
	fmpz_init(e);
	fmpz_init(nm);
	fmpz_init(fa);
	fmpz_init(fb);
	fmpz_poly_init(gx);
	fmpz_poly_init(q);
	fmpz_poly_one(sx);
	fmpz_one(Nm);

	for(ulong i = 0; i < num; i++)
	{
		slong a = abPairs[i].r;
		ulong b = abPairs[i].p;
		fmpz_set_si(fa,a);
		fmpz_set_ui(fb,b);
		fmpz_poly_set_coeff_si(gx,0,a);
		fmpz_poly_set_coeff_si(gx,1,b);
		fmpz_poly_mul(sx,sx,gx);
		fmpz_poly_divrem(q,sx,sx,f);
		norm(nm,f,fa,fb);
		fmpz_mul(Nm,Nm,nm);
	}
	fmpz_sqrtrem(r,e,Nm);
	/*Check that Nm is indeed a square number.*/
	assert(fmpz_is_zero(e));
	fmpz_set(Nm,r);

	fmpz_clear(r);
	fmpz_clear(e);
	fmpz_clear(nm);
	fmpz_clear(fa);
	fmpz_clear(fb);
	fmpz_poly_clear(gx);
	fmpz_poly_clear(q);
}

void estimateUpperBoundForX(fmpz_t res, const fmpz_poly_t delta, const fmpz_t m, ulong d)
{
	fmpz_zero(res);
	fmpz_t pom,c;
	fmpz_init_set_ui(pom,1);
	fmpz_init(c);
	for(ulong i = 0; i < d; i++)
	{
		fmpz_poly_get_coeff_fmpz(c,delta,i);
		fmpz_abs(c,c);
		fmpz_sqrt(c,c);
		fmpz_addmul(res,c,pom);
		fmpz_mul(pom,pom,m);
	}
	fmpz_mul_ui(res,res,100);
	fmpz_clear(pom);
	fmpz_clear(c);
}

bool selectPrimesCoverX(ulong *primes, ulong &nprimes, fmpz_t upperBound, ulong d, const fmpz_poly_t f)
{
	n_primes_t iter;
	n_primes_init(iter);
	n_primes_jump_after(iter,1000);
	nprimes = 0;
	while(fmpz_cmp_ui(upperBound,1) > 0)
	{
		ulong p = n_primes_next(iter);
		assert(p <= MaxPrime);
		assert(nprimes < MaxSelectedPrimes);
		if(legal(f,p))
		{
			primes[nprimes++] = p;
			fmpz_cdiv_q_ui(upperBound,upperBound,p);
			//cout << p << endl;
			//cout << upperBound << endl;
		}
	}
	n_primes_clear(iter);
	return fmpz_cmp_ui(upperBound,1)<=0;
}

void computeSquareRoots(ulong *XmodPi, ulong *primes, ulong nprimes, fmpz_poly_t delta,
						const fmpz_poly_t f, const fmpz_t m, const fmpz_t Nm)
{
	for(ulong i = 0; i < nprimes; i++)
	{
		ulong p = primes[i];
		XmodPi[i] = computeSquareRoot(delta,f,p,m,Nm);
	}
}

void computePinvs(ulong *Pinv, const ulong *primes, ulong nprimes)
{
	fmpz_t fp,fq,f;
	fmpz_init(fp);
	fmpz_init(fq);
	fmpz_init(f);
	for(ulong i = 0; i < nprimes; i++)
	{
		ulong p = primes[i];
		fmpz_set_ui(fp,p);
		Pinv[i] = 1;
		for(ulong j = 0; j < nprimes; j++)
		{
			if(j == i) continue;
			ulong q = primes[j];
			fmpz_set_ui(fq,q);
			fmpz_invmod(f,fq,fp);
			ulong qinv = fmpz_get_ui(f);
			Pinv[i] = (Pinv[i]*qinv)%p;
		}
	}
	fmpz_clear(fp);
	fmpz_clear(fq);
	fmpz_clear(f);
}

void computePinvsModn(fmpz_t *Pinvmodn, const ulong *primes, ulong nprimes, const fmpz_t n)
{
	fmpz_t fp,f;
	fmpz_init(fp);
	fmpz_init(f);
	for(ulong i = 0; i < nprimes; i++)
	{
		fmpz_init(Pinvmodn[i]);
		ulong p = primes[i];
		fmpz_set_ui(fp,p);
		fmpz_invmod(Pinvmodn[i],fp,n);
	}
	fmpz_clear(fp);
	fmpz_clear(f);
}

void computeAXoP(double *AXoP,const ulong *Pinv,const ulong *XmodPi,const ulong *primes,ulong nprimes)
{
	for(ulong i = 0; i < nprimes; i++)
		AXoP[i] = Pinv[i] * XmodPi[i] / (double) primes[i];
}

void productMod(fmpz_t res, ulong *a, ulong k, const fmpz_t n)
{
	fmpz_one(res);
	for(ulong i = 0; i < k; i++)
	{
		fmpz_mul_ui(res,res,a[i]);
		fmpz_mod(res,res,n);
	}
}

void sumOfAXPmodN(fmpz_t res, const ulong *Pinv, const ulong *XmodPi, const fmpz_t *pinvmodn, const fmpz_t Pmodn,
				  const ulong *primes, ulong nprimes, const fmpz_t n)
{
	fmpz_t s;
	fmpz_init(s);
	fmpz_zero(res);
	for(ulong i = 0; i < nprimes; i++)
	{
		fmpz_set(s,pinvmodn[i]);
		fmpz_mul(s,s,Pmodn);
		fmpz_mul_ui(s,s,XmodPi[i]);
		fmpz_mul_ui(s,s,Pinv[i]);
		fmpz_add(res,res,s);
		fmpz_mod(res,res,n);
	}
	fmpz_clear(s);
}

void freeArrayOfFmpz(fmpz_t *array, ulong size)
{
	for(ulong i = 0; i < size; i++)
		fmpz_clear(array[i]);
}

bool NFS(const fmpz_t n)
{
	/*--------------------Select polynomial-----------------------------------*/
	fmpz_poly_t f,delta;
	fmpz_t m, y, Nm, upperBoundOfX,r,x,Pmodn,rPmodn,xxmodn,yymodn;
	ulong d;
	fmpz_poly_init(delta);
#ifdef PRINT_PROCESS
	cout << "Selecting polynomial..." << endl;
#endif
	selectPolynomial(n,f,m,d);

#ifdef PRINT_MDF
	printf("m = "); fmpz_print(m); printf("\n");
	printf("d = %lu",d); printf("\n");
	printf("f(x) = "); fmpz_poly_print(f); printf("\n");
#endif
	
    /*--choose the bound for smoothness---------------------------------------*/
#ifdef PRINT_PROCESS
	cout << "Choosing smoothness bound..." << endl;
#endif
	ulong smoothBound = boundForSmoothness(d,n);
#ifdef PRINT_SMOOTH_BOUND
	cout << "Smoothness bound is " << smoothBound << endl;
#endif
	
	/*-Prepare the rational base----------------------------------------------*/
	ulong RB[MaxPrimeBufSize], nRB = 0;
	double lRB[MaxPrimeBufSize];
#ifdef PRINT_PROCESS
	cout << "Preparing the rational base..." << endl;
#endif
	prepareRationalBase(RB,lRB,nRB,smoothBound);
#ifdef PRINT_RATIONAL_BASE
	cout << "Rational base: " << endl;
	printListOfNumbers(RB,nRB,10);
#endif

	/*-Prepare the algebraic base---------------------------------------------*/
	MyPair AB[MaxPrimeBufSize];
	double lAB[MaxPrimeBufSize];
	ulong nAB = 0;
#ifdef PRITN_PROCESS
	cout << "Preparing the algebraic base..." << endl;
#endif
	prepareAlgebraicBase(AB,lAB,nAB,nRB,f);
#ifdef PRINT_ALGEBRAIC_BASE
	cout << "Algebraic base: " << endl;
	printListOfPairs(AB,nAB,10);
#endif

	/*-Prepare the quadratic base---------------------------------------------*/
	MyPair QB[MaxPrimeBufSize];
	ulong nQB = 0;
#ifdef PRINT_PROCESS
	cout << "Preparing the quadratic base..." << endl;
#endif
	prepareQuadraticBase(QB,nQB,smoothBound,smoothBound+20*log(smoothBound),f);
#ifdef PRINT_QUADRATIC_BASE
	cout << "Quadratic base: " << endl;
	printListOfPairs(QB,nQB,10);
#endif

	/*----------Sieve---------------------------------------------------------*/
	MyPair abPairs[2*MaxPrimeBufSize+1];
	ulong num = 2+nRB+nAB+nQB; /*Number of (a,b) pairs to search*/
	slong A = smoothBound*20, B = A/3;
#ifdef PRINT_PROCESS
	cout << "Sieving for " << num << " (a,b) pairs in region [" << -A << "~" << A << "]-[1," << B << "]"  << endl;
#endif
	//sieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs, num, smoothBound*20, m);
	latticeSieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs, num, A, B, m);
#ifdef PRINT_SELECTED_ABPAIRS
	cout << "Selected smooth (a,b) pairs: " << endl;
	printListOfPairs(abPairs,num,10);
#endif

	/*---------Form matrix----------------------------------------------------*/
	ulong I = num-1, J = num;
	nmod_mat_t matrix;
#ifdef PRINT_PROCESS
	cout << "Forming the matrix..." << endl;
#endif
	formMatrix(matrix,I,J,m,f,abPairs,RB,nRB,AB,nAB,QB,nQB);
#ifdef PRINT_MATRIX
	cout << "The matrix is: " << endl;
	nmod_mat_print_pretty(matrix);
#endif

	/*-------Solve the linear system------------------------------------------*/
	int *vec = new int[J];
#ifdef PRINT_PROCESS
	cout << "Solving the " << I << " by " << J << " linear system..." << endl;
#endif
	solveMatrix(matrix,I,J,vec);
	select(abPairs,vec,J,num); /*Select the pairs corresponding to 1 in vec*/
	nmod_mat_clear(matrix);
	delete[] vec;

#ifdef PRINT_SELECTED_SQUARE_ABPAIRS
	cout << "The selected (a,b) pairs whose product is square in both Z and Z[X]:" << endl;
	printListOfPairs(abPairs,num,10);
#endif

	/*---------Calculate prod(a+bm)-------------------------------------------*/
#ifdef PRINT_PROCESS
	cout << "Computing prod(a+bm)..." << endl;
#endif
	sqrtProductOfPairs(y,abPairs,num,m);
#ifdef PRINT_PROD_ABM
	fmpz_print(y);cout << endl;
#endif
	fmpz_mod(y,y,n);

	/*---------Calculate prod(a+b theta)--------------------------------------*/
	/*Compute the product of the norm of (a,b) pairs, used to select
	  beta or -beta when computing square root of delta mod p*/
#ifdef PRINT_PROCESS
	cout << "Computing prod(a+b theta)/f(theta)..." << endl;
#endif
	productOfPairs(delta,abPairs,num,f,Nm);









	/*=================Calculate x = phi(beta)================================*/
	/*This is a big project, divide it into several parts*/
#ifdef PRINT_PROCESS
	cout << "Computing phi(beta) mod n..." << endl;
#endif
	/*1. Estimate an upper bound for x~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	estimateUpperBoundForX(upperBoundOfX,delta,m,d);
#ifdef PRINT_UPPDER_BOUND
	cout << "Upper bound for X is ";
	fmpz_print(upperBoundOfX); cout << endl;
#endif
	/************************************************************************/



	/*2. Select p_i that Prod p_i > x~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	ulong primes[MaxSelectedPrimes]; ulong nprimes;
#ifdef PRINT_PROCESS
	cout << "----Selecting primes p_i such that prod p_i > x..." << endl;
#endif
	if(!selectPrimesCoverX(primes,nprimes,upperBoundOfX,d,f)) return false;
#ifdef PRINT_PRIMES
	cout << "--------Selected primes: " << endl;
	printListOfNumbers(primes,nprimes,10);
#endif
	/************************************************************************/



	/*3. Compute x_i = x mod p_i~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	ulong XmodPi[MaxSelectedPrimes];
#ifdef PRINT_PROCESS
	printf("----Computing x_i = x mod p_i...\n");
#endif
	computeSquareRoots(XmodPi,primes,nprimes,delta,f,m,Nm);
#ifdef PRINT_XI
	printListOfNumbers(XmodPi,nprimes,10);
#endif
	/************************************************************************/



	/**************************/
	/*4. Compute x mod n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	ulong Pinv[MaxSelectedPrimes]; /*Inverse of P_i mod p_i, where P_i = P/p_i, P = Prod p_i*/
	double AXoP[MaxSelectedPrimes];/*a_i*x_i/p_i, where a_i = Pinv[i]*/
	fmpz_t Pinvmodn[MaxSelectedPrimes];
#ifdef PRINT_PROCESS
	cout << "----Computing x from (x mod p_i) by Chinese remainder..." << endl;
#endif
	computePinvs(Pinv,primes,nprimes);
	computeAXoP(AXoP,Pinv,XmodPi,primes,nprimes);
	computePinvsModn(Pinvmodn,primes,nprimes,n);
	/*Let z be the result of Chinese remainder, then x = z mod P,
	  and z = x + rP, so r = (z-x)/P, since x < P, r = Floor(z/P)*/
	fmpz_set_d(r,doublesum(AXoP,0,nprimes));
	productMod(Pmodn,primes,nprimes,n);
	sumOfAXPmodN(x,Pinv,XmodPi,Pinvmodn,Pmodn,primes,nprimes,n);
	fmpz_mul(rPmodn,r,Pmodn);
	fmpz_mod(rPmodn,rPmodn,n);
	fmpz_sub(x,x,rPmodn);
	if(fmpz_cmp_ui(x,0)<0) fmpz_add(x,x,n);
	/*There might be cases where x < 0, then the x obtained above is not
	 * the real one, subtract P from it and mod n again.*/
	fmpz_mul(xxmodn,x,x);
	fmpz_mod(xxmodn,xxmodn,n);
	fmpz_mul(yymodn,y,y);
	fmpz_mod(yymodn,yymodn,n);
	if(!fmpz_equal(xxmodn,yymodn))
	{
		fmpz_sub(x,x,Pmodn);
		fmpz_mul(xxmodn,x,x);
		fmpz_mod(xxmodn,xxmodn,n);
	}
	/************************************************************************/
	fmpz_clear(m);
	fmpz_clear(Nm);
	fmpz_clear(upperBoundOfX);
	fmpz_clear(rPmodn);
	fmpz_poly_clear(f);
	fmpz_poly_clear(delta);
	freeArrayOfFmpz(Pinvmodn,nprimes);



	fmpz_t xpy, xmy, f1, f2;
	fmpz_init(xpy);
	fmpz_init(xmy);
	fmpz_init(f1);
	fmpz_init(f2);
	fmpz_add(xpy,x,y);
	fmpz_sub(xmy,x,y);
	fmpz_gcd(f1,xpy,n);
	fmpz_gcd(f2,xmy,n);
	/*Finally, we get our x mod n and y mod n. Time to sum up.*/
#if(!ONLY_RESULT)
	printf("x mod n = "); fmpz_print(x); printf("\n");
	printf("y mod n = "); fmpz_print(y); printf("\n");
	/*Check square of x and y*/
	printf("x^2 mod n = "); fmpz_print(xxmodn); printf("\n");
	printf("y^2 mod n = "); fmpz_print(yymodn); printf("\n");
	printf("x + y = "); fmpz_print(xpy); printf("\n");
	printf("x - y = "); fmpz_print(xmy); printf("\n");
	printf("GCD(x+y,n) = "); fmpz_print(f1); printf("\n");
	printf("GCD(x-y,n) = "); fmpz_print(f2); printf("\n");
#endif
	fmpz_clear(xpy);
	fmpz_clear(xmy);
	fmpz_clear(y);
	fmpz_clear(xxmodn);
	fmpz_clear(yymodn);
	/*Return true if any of f1 and f2 is a proper factor of n*/
	if(fmpz_cmp_ui(f1,1)>0 && fmpz_cmp(f1,n)<0)
	{
		fmpz_t nof1;
		fmpz_init(nof1);
		fmpz_fdiv_q(nof1,n,f1);
		fmpz_print(n);
		printf(" = ");
		fmpz_print(f1);
		printf(" * ");
		fmpz_print(nof1);
#if(!ONLY_RESULT)
		printf("\n");
#endif
		fmpz_clear(nof1);
		fmpz_clear(f1);
		fmpz_clear(f2);
		return true;
	}
	if(fmpz_cmp_ui(f2,1)>0 && fmpz_cmp(f2,n)<0)
	{
		fmpz_t nof2;
		fmpz_init(nof2);
		fmpz_fdiv_q(nof2,n,f2);
		fmpz_print(n);
		printf(" = ");
		fmpz_print(f2);
		printf(" * ");
		fmpz_print(nof2);
#if(!ONLY_RESULT)
		printf("\n");
#endif
		fmpz_clear(nof2);
		fmpz_clear(f1);
		fmpz_clear(f2);
		return true;
	}
	return false;
}

void testNFS(ulong e)
{
	time_t t = clock();
	fmpz_t n;
	fmpz_init_set_ui(n,2);
	fmpz_pow_ui(n,n,e);
	fmpz_add_ui(n,n,1);
	printf("2^%lu+1 = \t",e);
	for(int i = 0; i < 50; i++)
		if(NFS(n)) break;
	printf("\t");
	t = clock() - t;
	printf("%ld\n",t);
	fmpz_clear(n);
}

void init()
{
	srand((unsigned)time(0));
}

int main(int argc, char *argv[])
{
	init();
#if(!ONLY_RESULT)
	fmpz_t n; fmpz_init(n);
	char *ptr = NULL;
	char num[] = "268435457";
	if(argc > 1) ptr = argv[1];
	else ptr = num;
	fmpz_set_str(n,ptr,10);
	printf("n = "); fmpz_print(n); printf("\n");
	int Tries = 50;
	for(int i = 0; i < Tries; i++)
	{
		printf("--------------------------------------------------------------------------------\n");
		printf("Trying for the %d time...\n",i+1);
		if(NFS(n)) break;
		printf("--------------------------------------------------------------------------------\n");
	}
	fmpz_clear(n);
#else
	testNFS(10);
	testNFS(23);
	testNFS(30);
	testNFS(41);
	testNFS(50);
	testNFS(61);
	testNFS(71);
	testNFS(82);
	testNFS(90);
	testNFS(101);
#endif
	return 0;
}
