#include "GNFS.h"
#include "GNFS-lattice.h"

#if 0
void gaussianLatticeReduce(MyPair &u, MyPair &v)
{
	slong lu = length(u), lv = length(v);
	do{
		if(lu > lv)
		{
			MyPair tmp = u;
			u = v;
			v = tmp;
			lu = length(u);
			lv = length(v);
		}
		if(lu == 0) return;
		slong k = (u.p*v.p+u.r*v.r)/lu;
		slong k1 = k + 1;
		slong l = (v.p-k*u.p)*(v.p-k*u.p)+(v.r-k*u.r)*(v.r-k*u.r);
		slong l1 = (v.p-k1*u.p)*(v.p-k1*u.p)+(v.r-k1*u.r)*(v.r-k1*u.r);
		if(l > l1)
		{
			v.p -= k1*u.p;
			v.r -= k1*u.r;
		}
		else
		{
			v.p -= k*u.p;
			v.r -= k*u.r;
		}
		lv = length(v);
	} while(lu > lv);
}
#else
void gaussianLatticeReduce(MyPair &u, MyPair &v)
{
	v.r %= u.r;
}
#endif

/**
 *	Given C and D, and two vectors s, t, find the bound E1,E2 for s, F1,F2 for t, such that
 *	the region {es+ft | E1<=e<=E2, F1<=f<=F2} covers the region {(c,d) | -C<=c<=C, -D<=d<=D}
 *	as 'good' as possible.
 *
 *	Here the algorithm is simply taken as: stretch the vectors s,t until reaching the bound
 */
void getBound(slong &E1, slong &E2, slong &F1, slong &F2, slong C, slong D, MyPair s, MyPair t)
{
	if(s.r == 0) E1 = D/abs(s.p);
	else if(s.p == 0) E1 = C/abs(s.r);
	else
	{
		E1 = C/abs(s.r);
		if(D/abs(s.p) < E1) E1 = D/abs(s.p);
	}
	E1 *= boundaryFactor;
	E2 = E1;
	E1 = -E1;
	if(t.r == 0) F2 = D/abs(t.p);
	else if(t.p == 0) F2 = C/abs(t.r);
	else
	{
		F2 = C/abs(t.r);
		if(D/abs(t.p) < F2) F2 = D/abs(t.p);
	}
	F2 = F2 * boundaryFactor;
	F1 = 1;
#ifdef DEBUG
	assert(E2 > 0 && F2 > 0);
#endif
}

