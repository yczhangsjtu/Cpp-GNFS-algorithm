#include <cstdlib>
#include "GNFS.h"

/**
 *	Compute the square root of Prod(a+bm) for all the pairs (a,b) in an array
 */
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

/**
 *	Compute product of polynomials (a+bx) modular f for all the pairs (a,b) in the array.
 */
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

/**
 *	Compute the square root of a polynomial sx modular polynomial f in Fp[x].
 */
ulong computeSquareRoot(const fmpz_poly_t sx, const fmpz_poly_t f, slong p, const fmpz_t m, const fmpz_t Nm)
{
#ifdef DEBUG
	assert(p > 1);
#endif

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
#ifdef DEBUG
	assert(nmod_poly_is_one(Sxp));
#endif
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

/**
 *	We use Chinese remainder theorem to compute x = phi(beta), first we need to find
 *	a set of primes p, such that their products is larger than x.
 */

/**
 *	First estimate an upperbound for x.
 */
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

/**
 *	Then select the primes by repeatedly dividing the upperbound by larger and
 *	larger primes.
 */
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
		if(irreducibleMod(f,p))
		{
			primes[nprimes++] = p;
			fmpz_cdiv_q_ui(upperBound,upperBound,p);
		}
	}
	n_primes_clear(iter);
	return fmpz_cmp_ui(upperBound,1)<=0;
}

/**
 *	For each selected prime p, compute square root of polynomial delta mod f mod p,
 *	then substitute m in the result polynomial, we will get a square root of
 *	phi(beta^2) mod p, which are the x mod p[i]'s used in the CRT.
 */
void computeSquareRoots(ulong *XmodPi, ulong *primes, ulong nprimes, fmpz_poly_t delta,
						const fmpz_poly_t f, const fmpz_t m, const fmpz_t Nm)
{
	for(ulong i = 0; i < nprimes; i++)
	{
		ulong p = primes[i];
		XmodPi[i] = computeSquareRoot(delta,f,p,m,Nm);
	}
}

/**
 *	A step in the CRT: compute the inverse of P[i] modular p[i], where P[i] is P/p[i],
 *	and P is the product of all the p[i]'s.
 */
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

/**
 *	For each p[i], compute (double) a[i]x[i]/p[i], where a[i] is the inverse of P[i] mod p[i].
 */
void computeAXoP(double *AXoP,const ulong *Pinv,const ulong *XmodPi,const ulong *primes,ulong nprimes)
{
	for(ulong i = 0; i < nprimes; i++)
		AXoP[i] = Pinv[i] * XmodPi[i] / (double) primes[i];
}

/**
 *	Summation of a[i]x[i]p[i].
 */
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

	fmpz_t n, m, x, y, r, upperBoundOfX, Pmodn, rPmodn, xxmodn, yymodn, Nm;
	fmpz_t xpy, xmy, f1, f2;
	fmpz_poly_t f, delta;
	ulong num = 0, d, nprimes;
	ulong primes[MaxSelectedPrimes], XmodPi[MaxSelectedPrimes], Pinv[MaxSelectedPrimes];
	double AXoP[MaxSelectedPrimes];
	fmpz_t Pinvmodn[MaxSelectedPrimes];
	MyPair abPairs[2*MaxPrimeBufSize+1];

	fmpz_init(n);
	fmpz_init(m);
	fmpz_init(x);
	fmpz_init(y);
	fmpz_init(r);
	fmpz_init(Pmodn);
	fmpz_init(rPmodn);
	fmpz_init(xxmodn);
	fmpz_init(yymodn);
	fmpz_init(Nm);
	fmpz_init(upperBoundOfX);
	fmpz_init(xpy);
	fmpz_init(xmy);
	fmpz_init(f1);
	fmpz_init(f2);
	fmpz_poly_init(f);
	fmpz_poly_init(delta);
	fmpz_fread(input,n);
	fmpz_fread(input,m);
	fmpz_poly_fread(input,f);
	d = fmpz_poly_degree(f);
	fscanf(input,"%lu",&num);
	for(slong i = 0; i < num; i++)
		fscanf(input,"%ld%ld",&abPairs[i].r,&abPairs[i].p);
	
	sqrtProductOfPairs(y,abPairs,num,m);
	fmpz_mod(y,y,n);
	productOfPairs(delta,abPairs,num,f,Nm);
	estimateUpperBoundForX(upperBoundOfX,delta,m,d);

	if(!selectPrimesCoverX(primes,nprimes,upperBoundOfX,d,f)) exit(-1);
	computeSquareRoots(XmodPi,primes,nprimes,delta,f,m,Nm);
	computePinvs(Pinv,primes,nprimes);
	computeAXoP(AXoP,Pinv,XmodPi,primes,nprimes);
	computePinvsModn(Pinvmodn,primes,nprimes,n);
	fmpz_set_d(r,doublesum(AXoP,0,nprimes));
	productMod(Pmodn,primes,nprimes,n);
	sumOfAXPmodN(x,Pinv,XmodPi,Pinvmodn,Pmodn,primes,nprimes,n);
	fmpz_mul(rPmodn,r,Pmodn);
	fmpz_mod(rPmodn,rPmodn,n);
	fmpz_sub(x,x,rPmodn);
	if(fmpz_cmp_ui(x,0)<0) fmpz_add(x,x,n);
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

	fmpz_init(xpy);
	fmpz_init(xmy);
	fmpz_init(f1);
	fmpz_init(f2);
	fmpz_add(xpy,x,y);
	fmpz_sub(xmy,x,y);
	fmpz_gcd(f1,xpy,n);
	fmpz_gcd(f2,xmy,n);
#if(DEBUG)
	/*Finally, we get our x mod n and y mod n. Time to sum up.*/
	fprintf(output,"x mod n = "); fmpz_fprint(output,x); fprintf(output,"\n");
	fprintf(output,"y mod n = "); fmpz_fprint(output,y); fprintf(output,"\n");
	/*Check square of x and y*/
	fprintf(output,"x^2 mod n = "); fmpz_fprint(output,xxmodn); fprintf(output,"\n");
	fprintf(output,"y^2 mod n = "); fmpz_fprint(output,yymodn); fprintf(output,"\n");
	fprintf(output,"x + y = "); fmpz_fprint(output,xpy); fprintf(output,"\n");
	fprintf(output,"x - y = "); fmpz_fprint(output,xmy); fprintf(output,"\n");
	fprintf(output,"GCD(x+y,n) = "); fmpz_fprint(output,f1); fprintf(output,"\n");
	fprintf(output,"GCD(x-y,n) = "); fmpz_fprint(output,f2); fprintf(output,"\n");
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
		fmpz_fprint(output,n);
		fprintf(output," = ");
		fmpz_fprint(output,f1);
		fprintf(output," * ");
		fmpz_fprint(output,nof1);
		fprintf(output,"\n");
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
		fmpz_fprint(output,n);
		fprintf(output," = ");
		fmpz_fprint(output,f2);
		fprintf(output," * ");
		fmpz_fprint(output,nof2);
		fprintf(output,"\n");
		fmpz_clear(nof2);
		fmpz_clear(f1);
		fmpz_clear(f2);
		return true;
	}
	return false;

	fmpz_clear(n);
	fmpz_clear(m);
	fmpz_clear(x);
	fmpz_clear(y);
	fmpz_clear(r);
	fmpz_clear(Pmodn);
	fmpz_clear(rPmodn);
	fmpz_clear(xxmodn);
	fmpz_clear(yymodn);
	fmpz_clear(Nm);
	fmpz_clear(upperBoundOfX);
	fmpz_clear(xpy);
	fmpz_clear(xmy);
	fmpz_clear(f1);
	fmpz_clear(f2);
	fmpz_poly_clear(f);
	fmpz_poly_clear(delta);

	fclose(input);
	fclose(output);
	return 0;
}
