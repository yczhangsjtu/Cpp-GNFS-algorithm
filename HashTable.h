#ifndef __HASHTABLE_H__
#define __HASHTABLE_H__

#include "GNFS.h"

/**
 *	Pair used in the HashTable mapping integer pairs to doubles
 */
typedef struct PairDouble
{
	MyPair p;
	double v;
	PairDouble(MyPair a, double b){p=a;v=b;}
	PairDouble():p(){v=0.0;}
} PairDouble;

inline bool operator==(const MyPair &a, const MyPair &b)
{
	return a.r == b.r && a.p == b.p;
}

inline slong length(MyPair u)
{
	return u.p*u.p+u.r*u.r;
}


/*******************************************************************************
 *	Definition of HashTable class
 ******************************************************************************/
class HashTable{
public:
	ulong _T;
	HashTable(ulong T){
		if(T > MaxT) T = MaxT;
		_T = T;
		_data = new list<PairDouble>[T];
	}
	~HashTable(){
		if(size()/_T >= 3)
			cerr << endl << "Warning: Load factor = "
				 << (double)size()/_T << endl;
		delete []_data;
	}
	list<PairDouble> * _data;
	list<MyPair> _list;
	inline ulong size(){return _list.size();}
	inline slong hash(const MyPair &p)
	{
		return (p.r+p.p+p.r*(_T/2))%_T;
	}
	void set(const MyPair &p, double v)
	{
		/* Map pair p to value v, if not exist, insert it into the table */
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
		{
			if(iter->p==p)
			{
				iter->v=v;
				return;
			}
		}
		_data[s].push_back(PairDouble(p,v));
		_list.push_back(p);
	}
	inline void insert(const MyPair &p)
	{
		/* Insert pair p into table, used when the table is used as a set and
		 * the value is unimportant. */
		set(p,0.0);
	}
	void add(const MyPair &p, double v)
	{
		/* Add v to the value corresponding to p, if p does not exist, insert it
		 * and set the value to v. */
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
			if(iter->p==p)
			{
				iter->v+=v;
				return;
			}
		_data[s].push_back(PairDouble(p,v));
		_list.push_back(p);
	}
	double get(const MyPair &p)
	{
		/* Get the value of pair p. If not exist, return 0. */
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
			if(iter->p==p) return iter->v;
		return 0.0;
	}
	double find(const MyPair &p)
	{
		/* Return if p is contained in the table. */
		slong s = hash(p);
		list<PairDouble>::iterator iter = _data[s].begin();
		for(;iter!=_data[s].end();iter++)
			if(iter->p==p) return true;
		return false;
	}
};

#endif
