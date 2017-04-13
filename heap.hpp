// kill HEAPDONE
// add bool returnu

#ifndef _HEAP_HPP_
#define _HEAP_HPP_

#include <vector>
#include <utility>
#include <iostream>
#include <cstdlib>
#include <stdio.h>

using idx_t = std::size_t;

/// Heap class, sorted on S, with auxilary data stored in T
/// Much better to use the built-in priority queue, but we retain our own implementation a sorted heap here for posterity.
template<class S, class T>
class Heap
{
private:
  // Private data
  std::vector<std::pair<S, T> > h;
  idx_t numelt;
  idx_t capacity;
  idx_t blocksize;

  // Private member functions
  void flip(idx_t const ix1, idx_t const ix2);
  idx_t par(idx_t const ix);
  idx_t lch(idx_t const ix);
  idx_t rch(idx_t const ix);

public:
  Heap( ) = delete; // no empty constructor
  explicit Heap(Heap const &input) = delete; // no copy constructor
  
  explicit Heap(idx_t const initcapacity);
  Heap(idx_t const initcapacity, idx_t const _blocksize);
  void addToHeap(T const &t, S const &s);
  void addToHeap(std::pair<S,T> const &st);
  void showHeap() const;
  bool popFromHeap(std::pair<S, T>& st);
  inline idx_t showCapacity() const;
  inline idx_t numberOfElements() const;
}; 

/// Compute parent index for ix
/// \param[in] ix : Child index
/// \return         Parent index
template<class S, class T>
inline idx_t Heap<S,T>::par(idx_t const ix)
{ return((ix-1)/2); }

/// Compute left child index for ix
/// \param[in] ix : Parent index
/// \return         Left child index
template<class S, class T>
inline idx_t Heap<S,T>::lch(idx_t const ix)
{ return(2*ix+1); }

/// Compute right child index for ix
/// \param[in] ix : Parent index
/// \return         Right child index
template<class S, class T>
inline idx_t Heap<S,T>::rch(idx_t const ix)
{ return(2*ix+2); }

/// Get heap capacity
/// \return Heap<S,T> capacity
template<class S, class T>
inline idx_t Heap<S,T>::showCapacity() const
{ return(capacity); }

/// Get number of elements in heap
/// \return Number of elements
template<class S, class T>
inline idx_t Heap<S,T>::numberOfElements() const
{ return(numelt); }

#endif
