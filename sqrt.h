#ifndef __SQRT_H__
#define __SQRT_H__

#include "mypair.h"

void sqrtProductOfPairs(fmpz_t s, const MyPair *pairs, ulong num, const fmpz_t m);
void productOfPairs(fmpz_poly_t sx, const MyPair *abPairs, ulong num, const fmpz_poly_t f, fmpz_t Nm);
ulong computeSquareRoot(const fmpz_poly_t sx, const fmpz_poly_t f, slong p, const fmpz_t m, const fmpz_t Nm);
void estimateUpperBoundForX(fmpz_t res, const fmpz_poly_t delta, const fmpz_t m, ulong d);
bool selectPrimesCoverX(ulong *primes, ulong &nprimes, fmpz_t upperBound, ulong d, const fmpz_poly_t f);
void computeSquareRoots(ulong *XmodPi, ulong *primes, ulong nprimes, fmpz_poly_t delta,
						const fmpz_poly_t f, const fmpz_t m, const fmpz_t Nm);
void computePinvs(ulong *Pinv, const ulong *primes, ulong nprimes);
void computePinvsModn(fmpz_t *Pinvmodn, const ulong *primes, ulong nprimes, const fmpz_t n);
void computeAXoP(double *AXoP,const ulong *Pinv,const ulong *XmodPi,const ulong *primes,ulong nprimes);
void sumOfAXPmodN(fmpz_t res, const ulong *Pinv, const ulong *XmodPi, const fmpz_t *pinvmodn, const fmpz_t Pmodn,
				  const ulong *primes, ulong nprimes, const fmpz_t n);

#endif
