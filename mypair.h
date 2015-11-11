#ifndef __MYPAIR_H__
#define __MYPAIR_H__

#include <flint/fmpz.h>
/**
 *	Pair of integer
 */
typedef struct MyPair
{
	slong r;
	slong p;
	MyPair(slong a, slong b){r=a; p=b;}
	MyPair(){r=0;p=0;}
} MyPair;

inline bool operator<(MyPair p1, MyPair p2)
{
	if(p1.r < p2.r) return true;
	if(p1.r > p2.r) return false;
	return p1.p < p2.p;
}

inline bool operator==(MyPair p1, MyPair p2)
{
	return p1.r == p2.r && p1.p == p2.p;
}

#endif
