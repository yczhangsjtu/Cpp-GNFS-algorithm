#ifndef __POLY_H__
#define _POLY_H__

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod_poly.h>
#include "mypair.h"

extern int MaxPrime;
void norm(fmpz_t nm, const fmpz_poly_t f, const fmpz_t a, const fmpz_t b);
bool isSmooth(const fmpz_t num, const ulong *ps, ulong np);
bool isSmooth(const fmpz_t num, const MyPair *ps, ulong np);
void rootsMod(const fmpz_poly_t f, ulong p, ulong *roots, ulong &nroot);
int Leg(slong a, ulong p);
bool irreducibleMod(const fmpz_poly_t fx, ulong p);
void getMaxCoeff(fmpz_t m, const fmpz_poly_t f);
bool testPolynomial(const fmpz_poly_t f);
void selectNonResidual(nmod_poly_t px, const nmod_poly_t fx, ulong p, const fmpz_t e, ulong d);
ulong computeOrder(const nmod_poly_t px, const nmod_poly_t fx);

#endif
