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

#define RANDOMIZE 0
#define DEBUG 0
#define PRINT_PROCESS 1
#define PRINT_SIEVE_PROCESS 1
//#define SLOW_PRINT_SIEVE_PROCESS 1
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

/*******************************************************************************
 *	Global Constants
 ******************************************************************************/
const int MaxPrimeBufSize = 30000;
const int DefaultMaxPrime = 20000000;
const double DefaultThreshold = 5.0;
const int DefaultAfactor = 20;

#endif
