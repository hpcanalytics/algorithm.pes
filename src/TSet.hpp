/* 
 * File:   TSet.hpp [Modified version of "Set.hpp" supporting template data type]
 * Author: maleq
 * Modification made by: Md Hasanuzzaman Bhuiyan
 * Created on November 11, 2008, 4:35 PM
 * Modified on February 14, 2012
 */

#ifndef _SET_HPP
#define	_SET_HPP 

#include "utility.hpp"
				 
typedef  long int  Vsize; 

/*
inline int  Compare(const void *a, const void *b) 
{
	return *(ULI *)a < *(ULI *)b;
}
*/

template <typename DType>
class Set {
public:
	Vsize   size;
	Vsize   maxsize;
	DType 	*set;

	Set() {
		size = 0; 
		maxsize = 0;
		set = NULL;
	}
	
	~Set() {
		if (set!=NULL) delete [] set;
	}

	Set(Vsize s) {
		size = 0; 
		maxsize = s; 
		set = new DType[s];
	}
	
	void  init(Vsize s) {
		if (s <= 0) return;
		size = 0;
		maxsize = s;
		set = new DType[s];
	}
	
	void clear() {
		size = 0;
	}
	
	void destroy() {
		if (set!=NULL) delete [] set;
		set = NULL;
		size = maxsize = 0;
	}
	
	Vsize Size() const {
		return size;
	}
	
	Vsize getMaxSize() const {
		return maxsize;
	}
	
	void setSize(Vsize s)
	{
		size = s;
	}
	
	Vsize incrementSize()
	{
		return ++size;
	}
	
	Vsize decrementSize()
	{
		if(size > 0)
			size--;
		else
			cout << "\nError! Trying to decrement the value of the unsigned long size = 0\n";
		return size;
	}
	
	// check if element is a member of the set -- linear search in unsorted case
	bool member(DType element)
	{
		return LSearch(set, size, element);
	}
	
	// check if element is a member of the set -- if yes, return the index of it's appeareance, else return -1 -- linear search in unsorted case
	Vsize getPosition(DType element)
	{
		return LSearchPos(set, size, element);
	}
	
		// Insert an element; NO boundary check - it is caller's responsibility
	void insert(DType element) {
		set[size] = element; 
		size++;
	}

		// Insert an element; but, dynamically increase the size of no empty space
	void Dinsert(DType element) { //zhao zhao: dynamic expanding set size
		if(size == maxsize)
		{
			if(maxsize != 0)
			{
				if(maxsize >= 1024)
					maxsize = (Vsize)ceil(1.5L * size);
				else
					maxsize = (Vsize)ceil(2.0L * size);
				DType* tmp = set;
				set = new DType[maxsize];
				memcpy(set, tmp, sizeof(DType)*size);
				delete [] tmp;		
			}
			else
			{
				maxsize = 2;
				set = new DType[maxsize];
			}
		}
		set[size] = element; 
		size++;
	}
	
	
	
	DType &operator[](Vsize i) {
		return set[i];
	}
	
	int operator()(DType v) {		// membership function
		return BSearchD(set, 0, int(size)-1, v);  
	}
	
	int  Compare(const void *a, const void *b)
	{
		return *(DType *)a < *(DType *)b;
	}
	
	void  sort() {
		qsort(set, size, sizeof(DType), Compare);
	}
	
	void  finalize() {
        // qsort(set, size, sizeof(DType), Compare);
		if (size < maxsize) {
			DType *tmpset = set;
			set = new DType[size];
			for (Vsize i=0; i<size; i++) 
				set[i] = tmpset[i];
			delete [] tmpset;
		}
	}
	
	// both A and B must be sorted in descending order before calling intersect
	void intersect(Vsize i, Set &A, Set &B) {
        Vsize  k;
		size = 0;
		for (k=0; k<B.size && i<A.size; i++) {
			while (k<B.size && B.set[k] > A.set[i])  k++;
			if (k<B.size && B.set[k] == A.set[i]) insert(A.set[i]);
		}
	}
			
	// both A and B must be sorted in descending order before calling intersect
	void intersect(Set &A, Set &B) {
        Vsize  i, k;
		size = 0;
		
		for (i=0, k=0; k<B.size && i<A.size; i++) {
			while (k<B.size && B.set[k] > A.set[i])  k++;
			if (k<B.size && B.set[k] == A.set[i]) insert(A.set[i]);
		}
	}
	
	// both A and B must be sorted in descending order before calling intersect
	void intersect(Set &B) {
        Vsize  i, k;
		Vsize tsize = size;
		size = 0;
		
		for (i=0, k=0; k<B.Size() && i<tsize; i++) {
			while (k<B.Size() && B[k] > set[i])  k++;
			if (k<B.Size() && B[k] == set[i]) insert(set[i]);
		}
	}

	// both A and B must be sorted in descending order before calling set_union
	void set_union(Set &A, Set &B) {
        Vsize  i, k;
		size = 0;
		
		for (i=0, k=0; k<B.Size() && i<A.Size(); i++) {
			while (k<B.Size() && B[k] > A[i])  insert(B[k++]);
			if (k>=B.Size() || B[k] != A[i])  insert(A[i]);
		}
		while (i<A.Size())  insert(A[i++]);
		while (k<B.Size())  insert(B[k++]);
	}
	
	// both A and B must be sorted in descending order before calling set_union
	void set_union(Set &B) {
        Vsize  i, k;
		DType *tset = set;
		Vsize tsize = size; 
		
		size = 0;
		init(tsize+B.Size());
		
		for (i=0, k=0; k<B.Size() && i<tsize; i++) {
			while (k<B.Size() && B[k] > tset[i])  insert(B[k++]);
			if (k>=B.Size() || B[k] != tset[i])  insert(tset[i]);
		}
		while (i<tsize)  insert(tset[i++]);
		while (k<B.Size())  insert(B[k++]);
		
		if (tset) delete [] tset;
	}
		
	void print() {
		cout<< "Size "<< size <<": "; 
		for (int i=0; i<size; i++) 
			cout << set[i] << " "; 
		cout<<endl;
		
		cout<< "MaxSize "<< maxsize <<": "; 
		for (int i=0; i<maxsize; i++) 
			cout << set[i] << " "; 
		cout<<endl;
	}

	void intersect_idx(Set &A, Set &B) {
        Vsize  i, k;
		size = 0;
		
		for (i=0, k=0; k<B.size && i<A.size; i++) {
			while (k<B.size && B[k] > A[i])  k++;
			if (k<B.size && B[k] == A[i]) insert(i);
		}
	}

};

#endif	/* _SET_HPP */

