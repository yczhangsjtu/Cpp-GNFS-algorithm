#include "GNFS.h"

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
	ulong N = smoothBound*20;
	time_t start = clock();
#ifdef PRINT_PROCESS
	cout << "Sieving for " << num << " (a,b) pairs in region [-" << N << "," << N << "]-[1," << MaxB << "]"  << endl;
#endif
	sieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs, num, N, m);
	//latticeSieve(f, RB, lRB, nRB, AB, lAB, nAB, abPairs, num, A, B, m);
#ifdef PRINT_PROCESS
	cout << "Time used by the sieving process: " << (clock() - start)/1000 << "ms" << endl;
#endif
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
	printf("x mod n = "); fmpz_print(x); printf("\n");
	printf("y mod n = "); fmpz_print(y); printf("\n");
	/*Check square of x and y*/
	printf("x^2 mod n = "); fmpz_print(xxmodn); printf("\n");
	printf("y^2 mod n = "); fmpz_print(yymodn); printf("\n");
	printf("x + y = "); fmpz_print(xpy); printf("\n");
	printf("x - y = "); fmpz_print(xmy); printf("\n");
	printf("GCD(x+y,n) = "); fmpz_print(f1); printf("\n");
	printf("GCD(x-y,n) = "); fmpz_print(f2); printf("\n");
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
		printf("\n");
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
		printf("\n");
		fmpz_clear(nof2);
		fmpz_clear(f1);
		fmpz_clear(f2);
		return true;
	}
	return false;
}

