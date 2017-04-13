#include "heap.hpp"
#include <cassert>
// maintains heap sorted on S, T is auxilary data

/// Constructor for empty heap
/// \param[in] initcapacity : Initial capacity and blocksize
template<class S, class T>
Heap<S,T>::Heap(idx_t const initcapacity) :
  Heap(initcapacity, initcapacity)
{ }

/// Constructor for empty heap
/// \param[in] initcapacity : Initial capacity
/// \param[in] _blocksize   : Block size
template<class S, class T>
Heap<S,T>::Heap(idx_t const initcapacity, idx_t const _blocksize) :
  h(initcapacity),
  numelt(0),
  capacity(initcapacity),
  blocksize(_blocksize)
{ }

/// Swap two elements
/// \param[in] ix1 : First index location to swap
/// \param[in] ix2 : Second index location to swap
template<class S, class T>
void Heap<S,T>::flip(idx_t ix1, idx_t ix2)
{
  assert((ix1 < capacity) && (ix2 < capacity));
  std::swap(h[ix1], h[ix2]);
}

/// Add a new element to the help
/// \param[in] t : Data component of element to add to heap
/// \param[in] s : Sort component of element to add to heap
template<class S, class T>
void Heap<S,T>::addToHeap(T const &t, S const &s)
{
  addToHeap(std::make_pair<S,T>(s,t));
}

/// Add a new element to the help
/// \param[in] st : Element to add to heap
template<class S, class T>
void Heap<S,T>::addToHeap(std::pair<S,T> const &st)
{ // min-heap: minimum element at top
  // Step 0: add more space to heap if needed
  if(numelt == capacity)
  {
    capacity += blocksize;
    h.resize(capacity);
  }
  // Step 1: add element (i,d) to bottom of heap
  h[numelt] = st;
  
  // Step 2: upheap
  idx_t pos = numelt;
  while(pos > 0)
    if(h[pos].first < h[par(pos)].first) {
      flip(pos,par(pos));
      pos = par(pos);
    }
    else
      break;
  // Step 3: increment numelt
  numelt++;
}

/// Pop the top element from the heap and return it
/// \param[out] : Top element from the heap if return is true
/// \return     : True if pop succeeds
template<class S, class T>
bool Heap<S,T>::popFromHeap(std::pair<S,T> &st)
{ // pop top element
  if(numelt == 0)
    return false;
  
  st = h[0];
  numelt--;

  h[0] = h[numelt];
  idx_t pos = 0;
  while(true) {
    if(rch(pos) < numelt) {
      if(h[lch(pos)].d <= h[rch(pos)].d) {
	if(h[lch(pos)].d < h[pos].d) {
	  flip(lch(pos),pos);
	  pos = lch(pos);
	}
	else
	  return true;
      } else {
	if(h[rch(pos)].d < h[pos].d) {
	  flip(rch(pos),pos);
	  pos = rch(pos);
	} else
	  return true;
      }
    } else {
      if(lch(pos) < numelt) {
	if(h[lch(pos)].d < h[pos].d) {
	  flip(lch(pos),pos);
	  pos = lch(pos);
	}
	else
	  return true;
      }
      else
	return true;
    }
  }
}

/// Print the heap
template<class S, class T>
void Heap<S,T>::showHeap() const
{
  for(int ii=0; ii<numelt; ++ii)
    std::cout << "(" << h[ii].first << "," << h[ii].second << std::endl;
}
