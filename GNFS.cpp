#include <cstdio>
#include <iostream>
#include <cmath>
#include <set>
#include <cassert>
#include <ctime>
#include <cstdlib>

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/matrix.h>
#include <NTL/ZZX.h>

#define Int   ZZ
#define Intp  ZZ_p
#define IntX  ZZX
#define IntpX ZZ_pX

using namespace std;
using namespace NTL;

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

/*
 * Given integer n, return f, m, d
 * f: Integer array representing the coefficients of polynomial
 * m: the number selected
 * d: degree of the polynomial
 */
void selectPolynomial(Int n, IntX &f, Int *m, Int *d)
{
	*d = pow(3*log(n)/log(log(n)),0.333)+1;
	*m = pow(tol(n),1.0/tol(*d));
	f.SetLength(tol(*d)+1);
	for(int i = 0; i <= *d; i++)
	{
		Int tmp = n % *m;
		SetCoeff(f,i,tmp);
		n /= *m;
	}
}

/*
 * sieve_array: The array to accumulate
 * RB: Factor base
 * nRB: Size of factor base
 * sieve_len: sieve length
 */
void rationalSieve(Int *sieve_array, Int sieve_len, Int *RB, Int nRB, Int A, Int bm)
{
	/*Zerolize the sieve array*/
	for(long i = 0; i < sieve_len; i++)
		sieve_array[i] = A+bm+i;
	
	/*printf("Before sieve.\n");
	for(int i = 0; i < sieve_len; i++)
		printf("%ld ",sieve_array[i]);
	printf("\n");*/
	for(long i = 0; i < nRB; i++)
	{
		Int p = RB[i];
		Int f = A+bm;
		Int loc = f % p == 0? Int(0): p - (f % p);
		
		while(loc < sieve_len)
		{
			while(sieve_array[tol(loc)] != 0 && sieve_array[tol(loc)] % p ==0)
				sieve_array[tol(loc)] /= p;
			loc += p;
		}
	}
	/*printf("After sieve.\n");
	for(int i = 0; i < sieve_len; i++)
		printf("%ld ",sieve_array[i]);
	printf("\n");*/
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

void algebraicSieve(Int *nf_sieve_array, IntX &f, MyPair *AB, Int nAB, Int sieve_len, Int A, Int b)
{
	/*Zerolize the sieve array*/
	for(int i = 0; i < sieve_len; i++)
		nf_sieve_array[i] = norm(f,A+i,b);
	
	/*printf("Before sieve.\n");
	for(int i = 0; i < sieve_len; i++)
		printf("%ld ",nf_sieve_array[i]);
	printf("\n");*/
	for(long i = 0; i < nAB; i++)
	{
		Int p = AB[i].p;
		Int r = AB[i].r;
		Int g = A+b*r;
		Int loc = g % p == 0? Int(0): p - (g % p);
		
		while(loc < sieve_len)
		{
			while(nf_sieve_array[tol(loc)] != 0 && (g+loc)%p==0 && nf_sieve_array[tol(loc)] % p ==0)
				nf_sieve_array[tol(loc)] /= p;
			loc += p;
		}
	}
	/*printf("After sieve.\n");
	for(int i = 0; i < sieve_len; i++)
		printf("%ld ",nf_sieve_array[i]);
	printf("\n");*/
}

bool randombool(double p)
{
	return rand() < RAND_MAX * p;
}

void sieve(IntX &f, Int *RB, Int nRB, MyPair *AB, Int nAB, MyPair *result, Int num, Int N, Int m)
{
	Int loc(0);
	Int *r_sieve_array = new Int[tol(2*N+1)];
	Int *a_sieve_array = new Int[tol(2*N+1)];
	for(long b = 1; true; b++)
	{
		//printf("Sieve with b = %d\n",b);
		Int bm = b * m;
		//printf("Rational sieve.\n");
		rationalSieve(r_sieve_array, 2*N+1, RB, nRB, -N, bm);
		//printf("Algebraic sieve.\n");
		algebraicSieve(a_sieve_array, f, AB, nAB, 2*N+1, -N, Int(b));
		//printf("Searching for smooth pairs.\n");
		for(long i = 0; i < 2*N+1; i++)
		{
			if(abs(r_sieve_array[i]) == 1 && abs(a_sieve_array[i]) == 1 && GCD(i-N,Int(b))==1)
			{
				if(randombool(0.8))
				{
					result[tol(loc)] = MyPair(i-N, Int(b));
					loc++;
					if(loc >= num) break;
				}
			}
		}
		if(loc >= num) break;
	}
	num = loc;
	delete []r_sieve_array;
	delete []a_sieve_array;
}

void primeSieve(Int *result, Int *num, Int bound)
{
	Int *buf = new Int[tol(bound+1)];
	for(long i = 2; i <= bound; i++)
		buf[i] = Int(1);
	Int sbound = SqrRoot(bound);
	for(long i = 2; i <= sbound; i++)
	{
		if(buf[i] != 1) continue;
		for(int j = 2*i; j <= bound; j += i)
			buf[j] = 0;
	}
	*num = 0;
	for(long i = 2; i <= bound; i++)
	{
		if(tol(buf[i]))
		{
			result[tol(*num)] = i;
			(*num)++;
		}
	}
	delete []buf;
}

void rootsMod(IntX &f, Int p, Int *roots, Int *nroot)
{
	Int d(deg(f));
	Int q(p);
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

void solveMatrix(int **matrix, Int I, Int J, int *vec)
{
	Int piv[1024], qiv[1024];
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
	}

#if 0
	for(int i = 0; i < I; i++)
	{
		for(int j = 0; j < J; j++)
			printf("%d ",matrix[piv[i]][qiv[j]]);
		printf("\n");
	}
	printf("%d %d\n",minI,j);
	printf("\n");
#endif

	for(long jj = j; jj < J; jj++)
		vec[tol(qiv[jj])] = 1;
	for(long i = 0; i < minI; i++)
		vec[tol(qiv[i])] = 0;
	for(long jj = j; jj < J; jj++)
	{
		for(long i = 0; i < minI; i++)
			vec[tol(qiv[i])] ^= matrix[tol(piv[i])][tol(qiv[jj])];
	}
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
	Intp::init(Int(p));
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
	Int zp(p);
	Intp::init(zp);
	IntpX Fx(INIT_SIZE,tol(d+1));
	for(int i = 0; i <= d; i++)
	{
		Intp c(INIT_VAL,f[i]);
		SetCoeff(Fx,i,c);
	}
	IntpX Sx(INIT_SIZE,tol(d));
	for(int i = 0; i < d; i++)
	{
		Intp c(INIT_VAL, sx[i] % zp);
		SetCoeff(Sx,i,c);
	}
	Int q = power(zp,tol(d));
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
#if 0
	cout << endl;
	cout << "p = " << p << endl;
	cout << "n = 2^" << r << "*" << s << endl;
	cout << "eta = " << eta << endl;
	cout << "zeta = eta^s = " << zeta << endl;
	cout << "lambda_0 = " << lambda << endl;
	cout << "omega_0 = " << omega << endl;
#endif

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
#if 0
	cout << "lambda_n = " << lambda << endl;
	cout << "omega_n = " << omega << endl;
#endif
	/*Check that omega_n^2 = delta*/
	IntpX omega2 = (omega * omega) % Fx;
	assert(Sx == omega2);
	/*If the norm of omega is not equivalent to norm, negate it*/
	IntpX NN; Intp pp; Int ppp;
	modpower(NN,omega,Fx,q1/(p-1));
	GetCoeff(pp,NN,0);
	conv(ppp,pp);
	if(ppp != Nm % zp)
	{
		omega = -omega;
		modpower(NN,omega,Fx,q1/(p-1));
		GetCoeff(pp,NN,0);
		conv(ppp,pp);
		assert(ppp == Nm % zp);
	}
	/*Substitute m in*/
	Intp M;
	eval(M,omega,Intp(INIT_VAL,m));
	// cout << M << endl;
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

bool NFS(Int n)
{
	srand((unsigned)time(0));
	IntX f;
	Int m, d;
	selectPolynomial(n,f,&m,&d);
	printf("m=%ld\n",tol(m));
	printf("d=%ld\n",tol(d));
	printf("f(x)=");
	cout << f << endl;
	
    /*  choose y - the bound for smoothness  */
    double dlogd = tol(d)*log(tol(d));
    double temp = 1.0/tol(d) * log( tol(n) );
    double e = dlogd + sqrt( dlogd*dlogd + 4*temp*log(temp) );
    Int y((long)(1.5*exp( 0.5*e )));
    printf( "smoothness bound is %ld.\n", tol(y));
	
	/* Prepare the rational base */
	printf("Preparing the rational base.\n");
	Int RB[256], nRB(0);
	primeSieve(RB, &nRB, y);
	assert(nRB <= 256);
	printf("Rational base.\n");
	for(long i = 0; i < nRB; i++)
	{
		if(i && i % 10 == 0) printf("\n");
		printf("%6ld ",tol(RB[i]));
	}
	printf("\n");
	
	/* Prepare the algebraic base */
	printf("Preparing the algebraic base.\n");
	MyPair AB[256];
	Int nAB(0);
	for(long i = 0; i < nRB; i++)
	{
		Int p = RB[i];
		Int roots[16], nroot;
		rootsMod(f,p,roots,&nroot);
		for(long j = 0; j < nroot; j++)
		{
			AB[tol(nAB)] = MyPair(roots[j],p);
			nAB++;
		}
	}
	assert(nAB <= 256);
	printf("Algebraic base.\n");
	for(long i = 0; i < nAB; i++)
	{
		if(i && i % 10 == 0) printf("\n");
		printf("(%ld,%ld) ",tol(AB[i].r),tol(AB[i].p));
	}
	printf("\n");

	/* Prepare the quadratic base */
	MyPair QB[256]; Int nQB(0);
	Int pQB[256]; Int np;
	primeSieve(pQB, &np, 2*y);
	for(long i = 0; i < np; i++)
	{
		Int p = pQB[i];
		if(p <= y) continue;
		Int roots[16], nroot;
		rootsMod(f,p,roots,&nroot);
		for(long j = 0; j < nroot; j++)
		{
			QB[tol(nQB)] = MyPair(roots[j],p);
			nQB++;
		}
	}
	assert(nQB <= 256);
	printf("Quadratic base.\n");
	for(long i = 0; i < nQB; i++)
	{
		if(i && i % 10 == 0) printf("\n");
		printf("(%ld,%ld) ",tol(QB[i].r),tol(QB[i].p));
	}
	printf("\n");

	/*Sieve*/
	MyPair result[513];
	Int num = 2+nRB+nAB+nQB;
	printf("Sieving...\n");
	sieve(f, RB, nRB, AB, nAB, result, num, y*5, m);
	printf("Selected smooth (a,b) pairs\n");
	for(long i = 0; i < num; i++)
	{
		if(i && i % 10 == 0) printf("\n");
		printf("(%ld,%ld) ",tol(result[i].r),tol(result[i].p));
	}
	printf("\n");

	/*Form matrix*/
	printf("Forming the matrix...\n");
	Int I = num-1, J = num;
	int **matrix = new int*[tol(I)];
	for(long i = 0; i < I; i++)
		matrix[i] = new int[tol(J)];
	for(long j = 0; j < J; j++)
	{
		Int a = result[j].r;
		Int b = result[j].p;
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
			assert((a+b*s)%q != 0);
			matrix[tol(i+1+nRB+nAB)][j] = (l == 1? 0: 1);
		}
	}

	/*Print the matrix*/
#if 0
	for(int i = 0; i < I; i++)
	{
		for(int j = 0; j < J; j++)
			printf("%d ",matrix[i][j]);
		printf("\n");
	}
	printf("\n");
#endif
	int *vec = new int[tol(J)];
	solveMatrix(matrix,I,J,vec);

	for(long i = 0; i < I; i++)
		delete[] matrix[i];
	delete[] matrix;
	delete[] vec;

	/*Print the selected*/
	printf("The selected (a,b) pairs whose product is square in both Z and Z[X]\n");
	Int count(0);
	for(long i = 0; i < J; i++)
		if(vec[i])
		{
			printf("(%ld,%ld) ",tol(result[i].r),tol(result[i].p));
			count++;
			if(count % 10 == 0) printf("\n");
		}
	printf("\n");

	/*Calculate prod(a+bm)*/
	printf("Computing prod(a+bm)...\n");
	Int s(1);
	for(long i = 0; i < J; i++)
	{
		Int a = result[i].r;
		Int b = result[i].p;
		Int t(a+b*m);
		if(vec[i]) s *= t;
	}
	Int r = SqrRoot(s);
	/*Check that s is indeed a square number.*/
	assert(r*r == s);

	/*Calculate prod(a+b theta)*/
	printf("Computing prod(a+b theta)/f(theta)\n");
	Int Nm(1);
	IntX sx(1);
	Int phibeta2(0);
	for(long i = 0; i < J; i++)
	{
		Int a = result[i].r;
		Int b = result[i].p;
		IntX gx(INIT_SIZE,2);
		SetCoeff(gx,0,a);
		SetCoeff(gx,1,b);
		if(vec[i])
		{
			sx *= gx;
			sx %= f;
			Nm *= norm(f,a,b);
		}
	}
	printf("delta = Prod(a+b theta) = ");
	cout << sx << endl;

	/*Compute phi(beta^2) for checking*/
	Int pm(1);
	for(long i = 0; i < d; i++)
	{
		phibeta2 += sx[i] * pm;
		pm *= m;
	}
	/*Check that Nm is square*/
	Int NN = SqrRoot(Nm);
	assert(NN * NN == Nm);
	Nm = NN;
	/*Calculate x = phi(beta)*/
	/*1. Estimate an uppder bound for x*/
	printf("Computing phi(beta) mod n...\n");
	printf("----Selecting primes p_i such that prod p_i > x\n");
	Int S(0);
	Int pom(1);
	for(long i = 0; i < d; i++)
	{
		Int c = sx[i];
		if(c < 0) c = -c;
		c = SqrRoot(c) * 100;
		S += c * pom;
		pom *= m;
	}
	/*2. Select p_i that Prod p_i > x*/
	Int bound(100000);
	Int primes[10000], nprimes;
	primeSieve(primes,&nprimes,bound);
	Int ps[10000], nps(0);
	while(S > 1)
	{
		nprimes--;
		Int p = primes[tol(nprimes)];
		if(legal(f,p,d))
		{
			ps[tol(nps)] = p;
			nps++;
			S = S / p + 1;
		}
	}
#if 1
	printf("--------Selected primes: ");
	for(long i = 0; i < nps; i++)
		printf("%ld ",tol(ps[i]));
	printf("\n");
#endif
	/*3. Compute x_i = x mod p_i*/
	printf("----Computing x_i = x mod p_i\n");
	Int X[10000];
	for(long i = 0; i < nps; i++)
	{
		Int p = ps[i];
		X[i] = computeSquareRoot(sx,f,p,m,Nm);
		/*Check x^2 mod pi*/
		// cout << (Int(X[i])*Int(X[i])) % Int(p) << endl;
		// cout << phibeta2 % Int(p) << endl;
		// assert(Int(X[i])*Int(X[i]) % Int(p) == phibeta2 % Int(p));
	}
	/*4. Compute x mod n*/
	/*Compute the real answer*/
#if 0
	{
		Int z(0);
		Int P(1);
		for(int i = 0; i < nps; i++)
			P *= ps[i];
		Int Pmodn = P % Int(n);
		for(int i = 0; i < nps; i++)
		{
			Int ai = InvMod((P/ps[i])%Int(ps[i]),Int(ps[i]));
			z += ai * X[i] * (P/ps[i]);
		}
		Int x = z % P;
		cout << "Real P is " << P << endl;
		cout << "Real P mod n is " << P % Int(n) << endl;
		cout << "Real z is " << z << endl;
		cout << "Real r is " << z/P << endl;
		cout << "z - rP is " << z - (z/P)*P << endl;
		cout << "Real x is " << x << endl;
		cout << "Real x mod n is " << x % Int(n) << endl;
		cout << "Real x^2 mod n is " << (x*x) % n << endl;
		for(int i = 0; i < nps; i++)
			assert(x % Int(ps[i]) == X[i]);
		cout << "phi(beta^2) mod n = " << phibeta2 % Int(n) << endl;
	}
#endif
	/****************************/
	printf("----Computing x from (x mod p_i) by Chinese remainder.\n");
	Int Pinv[10000];
	for(long i = 0; i < nps; i++)
	{
		Int p = ps[i];
		Int Piinv(1);
		for(long j = 0; j < nps; j++)
		{
			if(j == i) continue;
			Int pp = ps[j];
			Int pinv = InvMod(pp%p,p);
			Piinv = (Piinv*pinv)%p;
		}
		Pinv[i] = Piinv;
	}

	double array[tol(nps)];
	for(int i = 0; i < nps; i++)
		array[i] = tol(Pinv[i] * X[i]) / (double) tol(ps[i]);
	Int rr((int)doublesum(array,Int(0),nps));
	// cout << "r = " << rr << endl;

	Int XX(0);
	Int Pmodn(1);
	for(long i = 0; i < nps; i++)
		Pmodn = (Pmodn * ps[i]) % n;
	// cout << "P mod n = " << Pmodn << endl;

	for(long i = 0; i < nps; i++)
		XX = (XX + Int(Pinv[i]) * Int(X[i]) * Int(Pmodn) * Int(InvMod(ps[i]%n,n))) % Int(n);
	XX = XX - rr * Pmodn;
	XX %= Int(n);
	if(XX < 0) XX += n;

	Int(YY) = Int(r)%Int(n);
	if(XX*XX%n != YY*YY%n) XX=(XX-Pmodn)%n;
	cout << "x mod n = " << XX << endl;
	cout << "y mod n = " << YY << endl;
	/*Check square of x and y*/
	cout << "x^2 mod n = " << XX * XX % n << endl;
	cout << "y^2 mod n = " << YY * YY % n << endl;
	cout << "x + y = " << XX+YY << endl;
	cout << "x - y = " << XX-YY << endl;
	Int f1 = GCD(XX+YY,Int(n));
	Int f2 = GCD(XX-YY,Int(n));
	cout << "GCD(x+y,n) = " << f1 << endl;
	cout << "GCD(x-y,n) = " << f2 << endl;
	return (f1 > 1 && f1 < Int(n)) || (f2 > 1 && f2 < Int(n));
}

int main(int argc, char *argv[])
{
	Int n(132163);
	if(argc > 1) n = Int(atoi(argv[1]));
	int Tries = 50;
	for(int i = 0; i < Tries; i++)
	{
		printf("--------------------------------------------------------------------------------\n");
		printf("Trying for the %d time...",i+1);
		if(NFS(n)) break;
		printf("--------------------------------------------------------------------------------\n");
	}
	return 0;
}
