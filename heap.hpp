#ifndef _HEAP_HPP_
#define _HEAP_HPP_

#include <vector>
#include <iostream>
#include <cstdlib>
#include <stdio.h>

#ifndef HEAPDONE
#define HEAPDONE -32323
#endif


#ifndef _HEAPSTRUCT_
#define _HEAPSTRUCT_
struct helt
{
  int i;
  double d;
  double aux[3];
};
#endif

using std::vector;
using std::cout;
using std::endl;

/// Much better to use the built-in priority queue, but we retain our own implementation a sorted heap here for posterity.
class Heap
{
private:
  // disable automatic constructors
  explicit Heap( ); // no empty constructor
  explicit Heap(const Heap &input); // no copy constructor
  // Private data
  vector<struct helt> h;
  int numelt;
  int capacity;
  int blocksize;
  const bool useAux;

  // Private member functions
  void flip(int ix1, int ix2);
  int par(const int ix);
  int lch(const int ix);
  int rch(const int ix);

public:
  
  explicit Heap(const int initcapacity);
  explicit Heap(const int initcapacity, const bool _useAux);
  inline void addToHeap(const int ii, const double dd);
  void addToHeap(const int ii, const double dd, const double *cp);
  void showHeap() const;
  struct helt popFromHeap();
  inline int showCapacity() const;
  inline int numberOfElements() const;
}; 

inline int Heap::par(const int ix)
{ return((ix-1)/2); }

inline int Heap::lch(const int ix)
{ return(2*ix+1); }

inline int Heap::rch(const int ix)
{ return(2*ix+2); }

inline void Heap::addToHeap(const int ii, const double dd)
{ addToHeap(ii,dd,NULL); }

inline int Heap::showCapacity() const
{ return(capacity); }

inline int Heap::numberOfElements() const
{ return(numelt); }

inline void showHelt(const struct helt &h);

inline void showHelt(const struct helt &h)
{ cout << "(" << h.i << "," << h.d << ";" << h.aux[0] << "," << h.aux[1] << ")" << endl; }


#endif
