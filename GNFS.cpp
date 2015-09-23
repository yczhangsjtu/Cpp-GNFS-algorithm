#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <limits.h>
#include <sstream>

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/matrix.h>
#include <NTL/ZZX.h>

#define Int   ZZ
#define Intp  ZZ_p
#define IntX  ZZX
#define IntpX ZZ_pX

#define DEBUG 1

#if(DEBUG)
#define PRINT_MDF
#define PRINT_PROCESS
//#define PRINT_RATIONAL_BASE
//#define PRINT_ALGEBRAIC_BASE
//#define PRINT_QUADRATIC_BASE
//#define PRINT_SELECTED_ABPAIRS
//#define PRINT_MATRIX
//#define PRINT_SELECTED_SQUARE_ABPAIRS
//#define PRINT_PRIMES
#endif

using namespace std;
using namespace NTL;

const int MaxPrimeBufSize = 4096;
const int MaxDimOfMatrix = 4096;
const int MaxB = 1024;

typedef struct MyPair
{
	Int r,p;
	MyPair(Int a, Int b){r=a;p=b;}
	MyPair(){r=0;p=0;}
} MyPair;

bool operator<(const MyPair &a, const MyPair &b)
{
	if(a.r < b.r) return true;
	if(a.r > b.r) return false;
	return a.p < b.p;
}

long tol(Int k)
{
	long K;
	conv(K,k);
	return K;
}

long tol(Intp k)
{
	long K;
	conv(K,k);
	return K;
}

Int findRootByDiv(Int n, Int l, Int r, long d)
{
	assert(power(l,d)<=n);
	if(l >= r-1) return l;
	Int m = (l+r)/2;
	Int mtod = power(m,d);
	if(mtod > n) return findRootByDiv(n,l,m,d);
	if(mtod < n) return findRootByDiv(n,m,r,d);
	return m;
}

Int findRootByDiv(Int n, long d)
{
	Int m(2), mtod = power(Int(2),d), twotod = power(Int(2),d);
	while(mtod <= n) {mtod *= twotod; m *= 2;}
	return findRootByDiv(n,m/2,m,d);
}

/*
 * Given integer n, return f, m, d
 * f: Integer array representing the coefficients of polynomial
 * m: the number selected
 * d: degree of the polynomial
 */
void selectPolynomial(Int n, IntX &f, Int &m, Int &d)
{
	assert(n > 0);
	assert(log(n) > 0);
	d = pow(3*log(n)/log(log(n)),0.333);
	if(d % 2 == 0) d++;
	while(true)
	{
		if(n > LONG_MAX)
			m = findRootByDiv(n,tol(d));
		else
			m = pow(tol(n),1.0/tol(d));
		if(m > 20)
		{
			do
				m -= Int(rand()) % (m/8);
			while(n % m == 0);
		}
		assert(power(m,tol(d)) < n);
		f.SetLength(tol(d)+1);
		for(int i = 0; i <= d; i++)
		{
			SetCoeff(f,i,n % m);
			n /= m;
		}
		if(f[tol(d)] == 1 && f[0] != 0) break;
	}
}

bool randombool(double p)
{
	return rand() < RAND_MAX * p;
}

Int norm(IntX &f, Int a, Int b)
{
	Int d(deg(f));
	if(a == 0)
	{
		Int p(1);
		for(int i = 0; i < d; i++)
			p *= -b;
		return p;
	}
	Int p(1);
	Int s(0);
	for(int i = 0; i < d; i++) p *= a;
	for(int i = tol(d); i >= 0; i--)
	{
		s += p * f[i];
		p = -p / a * b;
	}
	return s;
}

bool isSmooth(Int num, const Int *ps, Int np)
{
	for(long i = 0; i < np; i++)
	{
		Int p = ps[i];
		while(num > 0 && num % p == 0)
			num /= p;
	}
	return num == 1;
}

bool isSmooth(Int num, const MyPair *ps, Int np)
{
	for(long i = 0; i < np; i++)
	{
		Int p = ps[i].p;
		while(num > 0 && num % p == 0)
			num /= p;
	}
	return num == 1;
}

/*
 * sieve_array: The array to accumulate
 * RB: Factor base
 * nRB: Size of factor base
 * sieve_len: sieve length
 */
void rationalSieve(double *sieve_array, Int sieve_len, const Int *RB, const double *lRB, Int nRB, Int A, Int bm)
{
	/*Zerolize the sieve array*/
	for(long i = 0; i < sieve_len; i++)
	{
		if(A+bm+i == 0) continue;
		sieve_array[i] = -log(abs(A+bm+i));
	}
	
	for(long i = 0; i < nRB; i++)
	{
		Int p = RB[i];
		if(p < 0) cout << i << ' ' << RB[i] << ' ' << p << endl;
		Int f = A+bm;
		Int loc = f % p == 0? Int(0): p - (f % p);
		
		while(loc < sieve_len)
		{
			sieve_array[tol(loc)] += lRB[i];
			loc += p;
		}
	}
}

void algebraicSieve(double *nf_sieve_array, IntX &f, const MyPair *AB, const double *lAB, Int nAB, Int sieve_len, Int A, Int b)
{
	/*Zerolize the sieve array*/
	for(int i = 0; i < sieve_len; i++)
	{
		Int nm = abs(norm(f,A+i,b));
		if(nm==0) continue;
		assert(nm > 0);
		nf_sieve_array[i] = -log(nm);
	}
	
	for(long i = 0; i < nAB; i++)
	{
		Int p = AB[i].p;
		Int r = AB[i].r;
		Int g = A+b*r;
		Int loc = g % p == 0? Int(0): p - (g % p);
		
		while(loc < sieve_len)
		{
			assert(p > 0);
			nf_sieve_array[tol(loc)] += lAB[i];
			loc += p;
		}
	}
}

void sieve(IntX &f, const Int *RB, const double *lRB, Int nRB,
		   const MyPair *AB, const double *lAB, Int nAB, MyPair *abPairs, Int num, Int N, Int m)
{
	Int loc(0);
	double *r_sieve_array = new double[tol(2*N+1)];
	double *a_sieve_array = new double[tol(2*N+1)];
	for(long b = 1; b < MaxB; b++)
	{
		Int bm = b * m;
		rationalSieve(r_sieve_array, 2*N+1, RB, lRB, nRB, -N, bm);
		algebraicSieve(a_sieve_array, f, AB, lAB, nAB, 2*N+1, -N, Int(b));
		for(long i = 0; i < 2*N+1; i++)
		{
			Int a = i - N;
			if(r_sieve_array[i] >= -5.0 && a_sieve_array[i] >= -5.0 && GCD(a,Int(b))==1)
			{
				if(randombool(0.9))
				{
					if(isSmooth(abs(a+bm),RB,nRB) && isSmooth(abs(norm(f,a,(Int)b)),AB,nAB))
					{
						abPairs[tol(loc)] = MyPair(a, Int(b));
						loc++;
						if(loc >= num) break;
#ifdef PRINT_PROCESS
						cout << "\r" << loc*100/num << "%"; cout.flush();
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
}

void primeSieve(Int *result, Int *num, Int bound)
{
	bool *buf = new bool[tol(bound+1)];
	for(long i = 2; i <= bound; i++)
		buf[i] = true;
	Int sbound = SqrRoot(bound);
	for(long i = 2; i <= sbound; i++)
	{
		if(buf[i] != 1) continue;
		for(int j = 2*i; j <= bound; j += i)
			buf[j] = false;
	}
	*num = 0;
	for(long i = 2; i <= bound; i++)
	{
		if(buf[i])
		{
			long l = tol(*num);
			result[l] = i;
			(*num)++;
		}
	}
	delete []buf;
}

void rootsMod(IntX &f, Int p, Int *roots, Int *nroot)
{
	Int d(deg(f));
	Int q(p);
	assert(q > 1);
	Intp::init(q);
	IntpX g(INIT_SIZE,tol(d+1));
	for(int i = 0; i <= d; i++)
	{
		Intp c(INIT_VAL,f[i]);
		SetCoeff(g,i,c);
	}

	Vec< Pair<IntpX, long> > factors;
	CanZass(factors,g);
	Int n(factors.length());
	*nroot = 0;
	for(long i = 0; i < n; i++)
	{
		if(deg(factors[i].a) == 1)
		{
			Intp a,b;
			GetCoeff(a,factors[i].a,0);
			GetCoeff(b,factors[i].a,1);
			if(b != 0)
			{
				roots[tol(*nroot)] = Int(tol(-a/b));
				(*nroot)++;
			}
		}
	}
}

Int powermod(Int a, Int e, Int p)
{
	if(e == 0) return Int(1);
	if(e == 1) return a;
	Int pm = powermod(a,e/2,p);
	if(e % 2 == 0) return pm*pm%p;
	return pm*pm*a%p;
}

Int Leg(Int a, Int p) //Legender symbol
{
	if(a % p == 0) return Int(0);
	Int q = (p-1)/2;
	Int s = powermod(a,q,p);
	if(s == p-1) s = -1;
	if(s == 1-p) s = 1;
	assert(abs(s)==1);
	return s;
}

void Swap(Int &x, Int &y)
{
	Int tmp = x;
	x = y;
	y = tmp;
}

void solveMatrix(char **matrix, Int I, Int J, int *vec)
{
	Int piv[MaxDimOfMatrix], qiv[MaxDimOfMatrix];
	for(long i = 0; i < I; i++)
		piv[i] = i;
	for(long j = 0; j < J; j++)
		qiv[j] = j;
	Int minI(0);
	long j = 0;
	for(j = 0; j < J; j++)
	{
		if(minI >= I) break;
		bool allzero = true;
		for(long i = tol(minI); i < I; i++)
		{
			if(matrix[tol(piv[i])][tol(qiv[j])])
			{
				Swap(piv[tol(minI)],piv[i]);
				allzero = false;
				break;
			}
		}
		if(allzero)
		{
			bool allzero2 = true;
			for(long ii = tol(minI); ii < I; ii++)
			{
				for(long jj = j + 1; jj < J; jj++)
				{
					if(matrix[tol(piv[ii])][tol(qiv[jj])])
					{
						Swap(piv[tol(minI)],piv[ii]);
						Swap(qiv[j],qiv[jj]);
						allzero2 = false;
						break;
					}
				}
			}
			if(allzero2) break;
		}
		for(long ii = 0; ii < I; ii++)
		{
			if(ii == tol(minI) || matrix[tol(piv[ii])][tol(qiv[j])] == 0) continue;
			for(long jj = 0; jj < J; jj++)
				matrix[tol(piv[ii])][tol(qiv[jj])] ^= matrix[tol(piv[tol(minI)])][tol(qiv[jj])];
		}
		minI++;
#ifdef PRINT_PROCESS
		cout << "\r" << Int(j)*100/J << "%"; cout.flush();
#endif
	}

	for(long jj = j; jj < J; jj++)
		vec[tol(qiv[jj])] = 0;
	vec[tol(qiv[j])] = 1;
	for(long i = 0; i < minI; i++)
		vec[tol(qiv[i])] = matrix[tol(piv[i])][tol(qiv[j])];
}

void modpower(IntpX &px, IntpX sx, IntpX fx, Int e)
{
	if(e == 0)
	{
		px = 1;
		return;
	}
	IntpX mp;
	modpower(mp,sx,fx,e/2);
	if(e % 2 == 1)
		px = (mp * mp * sx) % fx;
	else
		px = (mp * mp) % fx;
}

bool legal(IntX &fx, Int p, Int d)
{
	assert(p > 1);
	Intp::init(p);
	IntpX Fx(INIT_SIZE,tol(d+1));
	for(long i = 0; i <= d; i++)
	{
		Intp c(INIT_VAL,fx[i]);
		SetCoeff(Fx,i,c);
	}
	Int zp(p), q = power(zp,tol(d));
	for(int i = 0; i < 20; i++)
	{
		IntpX px, mp;
		random(px,tol(d+1));
		modpower(mp,px,Fx,(q-1)/2);
		if(mp != tol(p-1) && mp != 1) return false;
	}
	return true;
}

void selectNonResidual(IntpX &px, IntpX fx, Int p, Int q, Int d)
{
	IntpX mp;
	while(true)
	{
		random(px,tol(d));
		modpower(mp,px,fx,(q-1)/2);
		if(mp == tol(p-1)) break;
	}
}

Int computeOrder(IntpX px, IntpX fx)
{
	for(long i = 0; i < 10000; i++)
	{
		if(px == 1) return Int(i);
		px = (px * px) % fx;
	}
	return Int(-1);
}

Int computeSquareRoot(IntX sx, IntX &f, Int p, Int m, Int Nm)
{
	Int d(deg(f));
	assert(p > 1);
	Intp::init(p);
	IntpX Fx(INIT_SIZE,tol(d+1));
	for(int i = 0; i <= d; i++)
	{
		Intp c(INIT_VAL,f[i]);
		SetCoeff(Fx,i,c);
	}
	IntpX Sx(INIT_SIZE,tol(d));
	for(int i = 0; i < d; i++)
	{
		Intp c(INIT_VAL, sx[i] % p);
		SetCoeff(Sx,i,c);
	}
	Int q = power(p,tol(d));
	Int q1 = q - 1;
	/*Check Sx is quadratic residual*/
	IntpX Sxp;
	modpower(Sxp,Sx,Fx,q1/2);
	assert(Sxp == 1);
	/********************************/
	Int r(0); Int s = q1;
	while(true)
	{
		if(s % 2 != 0) break;
		r++;
		s /= 2;
	}
	IntpX lambda, omega, zeta, eta;
	modpower(lambda,Sx,Fx,s);
	modpower(omega,Sx,Fx,(s+1)/2);
	selectNonResidual(eta,Fx,p,q,d);
	modpower(zeta,eta,Fx,s);

	while(lambda != 1)
	{
		Int m = computeOrder(lambda,Fx);
		assert(r > m);
		IntpX pzeta1 = zeta, pzeta2;
		for(long i = 0; i < r-m-1; i++)
			pzeta1 = pzeta1 * pzeta1;
		pzeta2 = pzeta1 * pzeta1;
		lambda = lambda * pzeta2;
		omega = omega * pzeta1;
	}
	/*Check that omega_n^2 = delta*/
	IntpX omega2 = (omega * omega) % Fx;
	assert(Sx == omega2);
	/*If the norm of omega is not equivalent to norm, negate it*/
	IntpX NN; Intp pp; Int ppp;
	modpower(NN,omega,Fx,q1/(p-1));
	GetCoeff(pp,NN,0);
	conv(ppp,pp);
	if(ppp != Nm % p)
	{
		omega = -omega;
		modpower(NN,omega,Fx,q1/(p-1));
		GetCoeff(pp,NN,0);
		conv(ppp,pp);
		assert(ppp == Nm % p);
	}
	/*Substitute m in*/
	Intp M;
	eval(M,omega,Intp(INIT_VAL,m));
	Int mm;
	conv(mm,M);
	return mm;
}

double doublesum(double array[], Int l, Int r)
{
	if(r-l == 0) return 0;
	if(r-l == 1) return array[tol(l)];
	return doublesum(array,l,(l+r)/2) + doublesum(array,(l+r)/2,r);
}

Int boundForSmoothness(long d, Int n)
{
	assert(d > 0);
	double dlogd = d*log(d);
	assert(n > 0);
	double temp = 1.0/d * log(n);
	assert(temp > 0);
	double e = dlogd + sqrt(dlogd*dlogd + 4*temp*log(temp));
	return Int(INIT_VAL,exp(0.5*e));
}

void prepareRationalBase(Int *RB, double *lRB, Int &nRB, Int bound)
{
	primeSieve(RB, &nRB, bound);
	assert(nRB <= MaxPrimeBufSize);
	for(int i = 0; i < nRB; i++)
		lRB[i] = log(RB[i]);
}

void prepareAlgebraicBase(MyPair *AB, double *lAB, Int &nAB, Int *primes, Int size, IntX &f)
{
	for(long i = 0; i < size; i++)
	{
		Int p = primes[i];
		Int roots[16], nroot;
		rootsMod(f,p,roots,&nroot);
		for(long j = 0; j < nroot; j++)
		{
			AB[tol(nAB)] = MyPair(roots[j],p);
			nAB++;
		}
	}
	assert(nAB <= MaxPrimeBufSize);
	for(int i = 0; i < nAB; i++)
		lAB[i] = log(AB[i].p);
}

void prepareQuadraticBase(MyPair *QB, Int &nQB, Int min, Int max, IntX &f)
{
	Int pQB[MaxPrimeBufSize]; Int np;
	primeSieve(pQB, &np, max);
	for(long i = 0; i < np; i++)
	{
		Int p = pQB[i];
		if(p <= min) continue;
		Int roots[16], nroot;
		rootsMod(f,p,roots,&nroot);
		for(long j = 0; j < nroot; j++)
		{
			QB[tol(nQB)] = MyPair(roots[j],p);
			nQB++;
		}
	}
	assert(nQB <= MaxPrimeBufSize);
}

void printListOfNumbers(Int *A, long s, long N)
{
	for(long i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) cout << endl;
		cout << setw(10) << tol(A[i]);
	}
	cout << endl;
}

void printListOfPairs(MyPair *A, long s, long N)
{
	for(long i = 0; i < s; i++)
	{
		if(N && i && i % N == 0) cout << endl;
		cout << "(" << A[i].r << "," << A[i].p << ") ";
	}
	cout << endl;
}

char **allocMatrix(long I, long J)
{
	char **matrix = new char*[I];
	for(long i = 0; i < I; i++)
		matrix[i] = new char[J];
	return matrix;
}

void formMatrix(char **matrix, Int I, Int J, Int m, ZZX &f, MyPair *abPairs,
				Int *RB, Int nRB, MyPair* AB, Int nAB, MyPair* QB, Int nQB)
{
	for(long j = 0; j < J; j++)
	{
		Int a = abPairs[j].r;
		Int b = abPairs[j].p;
		Int A = a+b*m;
		matrix[0][j] = A >= 0? 0: 1;
		for(long i = 0; i < nRB; i++)
		{
			Int p = RB[i];
			Int e(0);
			while(A != 0 && A % p ==0)
			{
				A /= p;
				e++;
			}
			matrix[i+1][j] = e % 2;
		}
		assert(abs(A) == 1);
		A = norm(f,a,b);
		for(long i = 0; i < nAB; i++)
		{
			Int p = AB[i].p;
			Int r = AB[i].r;
			Int e(0);
			while(A != 0 && A % p ==0 && (a+b*r)%p == 0)
			{
				A /= p;
				e++;
			}
			matrix[tol(i+1+nRB)][j] = e % 2;
		}
		assert(abs(A) == 1);
		for(long i = 0; i < nQB; i++)
		{
			Int s = QB[i].r;
			Int q = QB[i].p;
			Int l = Leg(a+b*s,q);
			if((a+b*s)%q==0)
			{
				cout << q << " divides " << a+b*s << endl;
				cout << a << ' ' << b << endl;
				A = norm(f,a,b);
				cout << A << endl;
				cout << A % q << endl;
			}
			assert((a+b*s)%q != 0);
			matrix[tol(i+1+nRB+nAB)][j] = (l == 1? 0: 1);
		}
	}
}

void printMatrix(char **matrix, long I, long J)
{
	for(long i = 0; i < I; i++)
	{
		for(long j = 0; j < J; j++)
			printf("%d ",matrix[i][j]);
		printf("\n");
	}
	printf("\n");
}

void freeMatrix(char **matrix, long I)
{
	for(long i = 0; i < I; i++)
		delete[] matrix[i];
	delete[] matrix;
}

void select(MyPair *pairs, int *vec, long I, Int &n)
{
	long loc = 0;
	for(long i = 0; i < I; i++)
	{
		if(vec[i]) pairs[loc++] = pairs[i];
	}
	n = Int(loc);
}

Int sqrtProductOfPairs(MyPair *pairs, long num, Int m)
{
	Int s(1);
	for(long i = 0; i < num; i++)
	{
		Int a = pairs[i].r;
		Int b = pairs[i].p;
		Int t(a+b*m);
		s *= t;
	}
	Int r = SqrRoot(s);
	/*Check that s is indeed a square number.*/
	assert(r*r == s);
	return r;
}

IntX productOfPairs(MyPair *abPairs, long num, ZZX &f, Int &Nm)
{
	IntX sx(1);
	Nm = Int(1);
	for(long i = 0; i < num; i++)
	{
		Int a = abPairs[i].r;
		Int b = abPairs[i].p;
		IntX gx(INIT_SIZE,2);
		SetCoeff(gx,0,a);
		SetCoeff(gx,1,b);
		MulMod(sx,sx,gx,f);
		Nm *= norm(f,a,b);
	}
	/*Check that Nm is square*/
	assert(Nm > 0);
	Int NN = SqrRoot(Nm);
	assert(NN * NN == Nm);
	Nm = NN;
	return sx;
}

Int estimateUpperBoundForX(ZZX &delta, Int m, Int d)
{
	Int S(0);
	Int pom(1); /*power of m*/
	for(long i = 0; i < d; i++)
	{
		Int c = delta[i];
		if(c < 0) c = -c;
		c = SqrRoot(c);
		S += c * pom;
		pom *= m;
	}
	return S * 100;
}

bool selectPrimesCoverX(Int *primes, Int &nprimes, Int upperBound, Int d, ZZX &f)
{
	Int ps[20000], nps, bound(200000);
	if(upperBound < 10000) bound = Int(10000);
	else bound = findRootByDiv(upperBound,10);
	if(bound < 10000) bound = Int(10000);
	if(bound > 200000) bound = Int(200000);
	nprimes = Int(0);
	primeSieve(ps,&nps,bound);
	while(upperBound > 1)
	{
		nps--;
		if(nps < 0) return false;
		Int p = ps[tol(nps)];
		if(legal(f,p,d))
		{
			primes[tol(nprimes)] = p;
			nprimes++;
			upperBound = upperBound / p + 1;
		}
	}
	return true;
}

void computeSquareRoots(Int *XmodPi, Int *primes, Int nprimes, IntX &delta,
						IntX &f, Int m, Int Nm)
{
	for(long i = 0; i < nprimes; i++)
	{
		Int p = primes[i];
		XmodPi[i] = computeSquareRoot(delta,f,p,m,Nm);
	}
}

void computePinvs(Int *Pinv, Int *primes, Int nprimes)
{
	for(long i = 0; i < nprimes; i++)
	{
		Int p = primes[i];
		Pinv[i] = Int(1);
		for(long j = 0; j < nprimes; j++)
		{
			if(j == i) continue;
			Int pinv = InvMod(primes[j]%p,p);
			Pinv[i] = (Pinv[i]*pinv)%p;
		}
	}
}

void computeAXoP(double *AXoP,Int *Pinv,Int *XmodPi,Int *primes,Int nprimes)
{
	for(int i = 0; i < nprimes; i++)
		AXoP[i] = tol(Pinv[i] * XmodPi[i]) / (double) tol(primes[i]);
}

Int productMod(Int *a, Int k, Int n)
{
	Int s(1);
	for(long i = 0; i < k; i++)
		MulMod(s,s,a[i],n);
	return s;
}

Int sumOfAXPmodN(Int *Pinv, Int *XmodPi, Int Pmodn, Int *primes, Int nprimes, Int n)
{
	Int x(0);
	for(long i = 0; i < nprimes; i++)
		x = (x + Pinv[i] * XmodPi[i] * Pmodn * InvMod(primes[i]%n,n)) % n;
	return x;
}

Int sumOfAXP(Int *Pinv, Int *XmodPi, Int P, Int *primes, Int nprimes)
{
	Int x(0);
	for(long i = 0; i < nprimes; i++)
		x += Pinv[i] * XmodPi[i] * P / primes[i];
	return x;
}

bool NFS(Int n)
{
	/*--------------------Select polynomial-----------------------------------*/
	IntX f;
	Int m, d;
#ifdef PRINT_PROCESS
	cout << "Selecting polynomial..." << endl;
#endif
	selectPolynomial(n,f,m,d);
#ifdef PRINT_MDF
	cout << "m = " << m << endl;
	cout << "d = " << d << endl;
	cout << "f(x) = " << f << endl;
#endif
	
    /*--choose the bound for smoothness---------------------------------------*/
#ifdef PRINT_PROCESS
	cout << "Choosing smoothness bound..." << endl;
#endif
	Int smoothBound = boundForSmoothness(tol(d),n);
#ifdef PRINT_SMOOTH_BOUND
	cout << "Smoothness bound is " << smoothBound << endl;
#endif
	
	/*-Prepare the rational base----------------------------------------------*/
	Int RB[MaxPrimeBufSize], nRB(0);
	double lRB[MaxPrimeBufSize];
#ifdef PRINT_PROCESS
	cout << "Preparing the rational base..." << endl;
#endif
	prepareRationalBase(RB,lRB,nRB,smoothBound);
#ifdef PRINT_RATIONAL_BASE
	cout << "Rational base: " << endl;
	printListOfNumbers(RB,tol(nRB),10);
#endif

	/*-Prepare the algebraic base---------------------------------------------*/
	MyPair AB[MaxPrimeBufSize];
	double lAB[MaxPrimeBufSize];
	Int nAB(0);
#ifdef PRITN_PROCESS
	cout << "Preparing the algebraic base..." << endl;
#endif
	prepareAlgebraicBase(AB,lAB,nAB,RB,nRB,f);
#ifdef PRINT_ALGEBRAIC_BASE
	cout << "Algebraic base: " << endl;
	printListOfPairs(AB,tol(nAB),10);
#endif

	/*-Prepare the quadratic base---------------------------------------------*/
	MyPair QB[MaxPrimeBufSize]; Int nQB(0);
#ifdef PRINT_PROCESS
	cout << "Preparing the quadratic base..." << endl;
#endif
	prepareQuadraticBase(QB,nQB,smoothBound,smoothBound+Int(INIT_VAL,20*log(smoothBound)),f);
#ifdef PRINT_QUADRATIC_BASE
	cout << "Quadratic base: " << endl;
	printListOfPairs(QB,tol(nQB),10);
#endif

	/*----------Sieve---------------------------------------------------------*/
	MyPair abPairs[2*MaxPrimeBufSize+1];
	Int num = 2+nRB+nAB+nQB; /*Number of (a,b) pairs to search*/
#ifdef PRINT_PROCESS
	cout << "Sieving for " << num << " (a,b) pairs..." << endl;
#endif
	sieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs, num, smoothBound*20, m);
#ifdef PRINT_SELECTED_ABPAIRS
	cout << "Selected smooth (a,b) pairs: " << endl;
	printListOfPairs(abPairs,tol(num),10);
#endif

	/*---------Form matrix----------------------------------------------------*/
	Int I = num-1, J = num;
	char **matrix = allocMatrix(tol(I),tol(J));
#ifdef PRINT_PROCESS
	cout << "Forming the matrix..." << endl;
#endif
	formMatrix(matrix,I,J,m,f,abPairs,RB,nRB,AB,nAB,QB,nQB);
#ifdef PRINT_MATRIX
	cout << "The matrix is: " << endl;
	printMatrix(matrix,tol(I),tol(J));
#endif

	/*-------Solve the linear system------------------------------------------*/
	int *vec = new int[tol(J)];
#ifdef PRINT_PROCESS
	cout << "Solving the " << I << " by " << J << " linear system..." << endl;
#endif
	solveMatrix(matrix,I,J,vec);
	freeMatrix(matrix,tol(I));
	select(abPairs,vec,tol(J),num); /*Select the pairs corresponding to 1 in vec*/
	delete[] vec;

#ifdef PRINT_SELECTED_SQUARE_ABPAIRS
	cout << "The selected (a,b) pairs whose product is square in both Z and Z[X]:" << endl;
	printListOfPairs(abPairs,tol(num),10);
#endif

	/*---------Calculate prod(a+bm)-------------------------------------------*/
#ifdef PRINT_PROCESS
	cout << "Computing prod(a+bm)..." << endl;
#endif
	Int y = sqrtProductOfPairs(abPairs,tol(num),m) % n;

	/*---------Calculate prod(a+b theta)--------------------------------------*/
	Int Nm(1); /*Compute the product of the norm of (a,b) pairs, used to select
				beta or -beta when computing square root of delta mod p*/
#ifdef PRINT_PROCESS
	cout << "Computing prod(a+b theta)/f(theta)..." << endl;
#endif
	IntX delta = productOfPairs(abPairs,tol(num),f,Nm);









	/*=================Calculate x = phi(beta)================================*/
	/*This is a big project, divide it into several parts*/
#ifdef PRINT_PROCESS
	cout << "Computing phi(beta) mod n..." << endl;
#endif
	/*1. Estimate an upper bound for x~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	Int upperBoundOfX = estimateUpperBoundForX(delta,m,d);
	/************************************************************************/



	/*2. Select p_i that Prod p_i > x~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	Int primes[10000], nprimes;
#ifdef PRINT_PROCESS
	cout << "----Selecting primes p_i such that prod p_i > x..." << endl;
#endif
	if(!selectPrimesCoverX(primes,nprimes,upperBoundOfX,d,f)) return false;
#ifdef PRINT_PRIMES
	cout << "--------Selected primes: " << endl;
	printListOfNumbers(primes,tol(nprimes),10);
#endif
	/************************************************************************/



	/*3. Compute x_i = x mod p_i~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	Int XmodPi[10000];
#ifdef PRINT_PROCESS
	printf("----Computing x_i = x mod p_i...\n");
#endif
	computeSquareRoots(XmodPi,primes,nprimes,delta,f,m,Nm);
	/************************************************************************/



	/*4. Compute x mod n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	Int Pinv[10000]; /*Inverse of P_i mod p_i, where P_i = P/p_i, P = Prod p_i*/
	double AXoP[tol(nprimes)];/*a_i*x_i/p_i, where a_i = Pinv[i]*/
	Int r; /*Let z be the result of Chinese remainder, then x = z mod P,
			and z = x + rP, so r = (z-x)/P, since x < P, r = Floor(z/P)*/
	Int x(0); /*This in fact is x mod n*/
	Int Pmodn(1); /*P mod n*/
#ifdef PRINT_PROCESS
	cout << "----Computing x from (x mod p_i) by Chinese remainder..." << endl;
#endif
	computePinvs(Pinv,primes,nprimes);
	computeAXoP(AXoP,Pinv,XmodPi,primes,nprimes);
	r = Int(INIT_VAL, doublesum(AXoP,Int(0),nprimes));
	Pmodn = productMod(primes,nprimes,n);
	x = (sumOfAXPmodN(Pinv,XmodPi,Pmodn,primes,nprimes,n) - r*Pmodn) % n;

	Int P(1);
	for(int i = 0; i < nprimes; i++)
		P *= primes[i];

	if(x < 0) x += n;
	/*There might be cases where x < 0, then the x obtained above is not
	 * the real one, subtract P from it and mod n again.*/
	if(x*x%n != y*y%n) x=(x-Pmodn)%n;
	/************************************************************************/



	/*Finally, we get our x mod n and y mod n. Time to sum up.*/
	cout << "x mod n = " << x << endl;
	cout << "y mod n = " << y << endl;
	/*Check square of x and y*/
	cout << "x^2 mod n = " << x * x % n << endl;
	cout << "y^2 mod n = " << y * y % n << endl;
	assert(x*x%n == y*y%n);
	cout << "x + y = " << x+y << endl;
	cout << "x - y = " << x-y << endl;
	Int f1 = GCD(x+y,n);
	Int f2 = GCD(x-y,n);
	cout << "GCD(x+y,n) = " << f1 << endl;
	cout << "GCD(x-y,n) = " << f2 << endl;
	/*Return true if any of f1 and f2 is a proper factor of n*/
	if(f1 > 1 && f1 < n)
	{
		cout << n << " = " << f1 << " * " << n/f1 << endl;
		return true;
	}
	if(f2 > 1 && f2 < n)
	{
		cout << n << " = " << f2 << " * " << n/f2 << endl;
		return true;
	}
	return false;
}

void init()
{
	srand((unsigned)time(0));
}

int main(int argc, char *argv[])
{
	init(); Int n;
	if(argc > 1)
	{
		istringstream is(argv[1]);
		is >> n;
	}
	else
	{
		char num[] = "1125899906842625";
		istringstream is(num);
		is >> n;
	}
	cout << "n = " << n << endl;
	int Tries = 50;
	for(int i = 0; i < Tries; i++)
	{
		printf("--------------------------------------------------------------------------------\n");
		printf("Trying for the %d time...\n",i+1);
		if(NFS(n)) break;
		printf("--------------------------------------------------------------------------------\n");
	}
	return 0;
}
