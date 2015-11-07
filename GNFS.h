#ifndef __GNFS_H__
#define __GNFS_H__

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

#define RANDOMIZE 0
#define DEBUG 0
#define PRINT_PROCESS 1
#define PRINT_SIEVE_PROCESS 1
//#define SLOW_PRINT_SIEVE_PROCESS 1

#include "algebraic.h"
#include "factorbase.h"
#include "linear.h"
#include "poly.h"
#include "polyselect.h"
#include "rational.h"
#include "sieve.h"
#include "sqrt.h"
#include "util.h"

#if(PRINT_PROCESS)
#define PRINT_MDF
//#define PRINT_SMOOTH_BOUND
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

/*******************************************************************************
 *	Global Constants
 ******************************************************************************/
const int MaxPrimeBufSize = 30000;
const int MaxPrime = 20000000;
const int MaxSelectedPrimes = 100000;
const int MaxB = 20240;
const ulong MaxT = 200000;
const int smoothfactor = 5;
const double threshold = 5.0;
const double boundaryFactor = 1.0;

#endif
