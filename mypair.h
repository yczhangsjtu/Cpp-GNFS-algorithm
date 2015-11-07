#ifndef __MYPAIR_H__
#define __MYPAIR_H__

#include "GNFS.h"

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

#endif
