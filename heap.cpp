#include "heap.hpp"
// maintains heap sorted by the double value (the int value is an index)

Heap::Heap(const int initcapacity) :
  useAux(false)
{
  h.resize(initcapacity);
  capacity = initcapacity;
  blocksize = initcapacity;
  numelt = 0;
}

Heap::Heap(const int initcapacity, const bool _useAux) :
  useAux(_useAux)
{
  h.resize(initcapacity);
  capacity = initcapacity;
  blocksize = initcapacity;
  numelt = 0;
}

void Heap::flip(int ix1, int ix2)
{
  struct helt tmp;
  if((ix1 >= 0) && (ix1 <= numelt) && (ix2 >= 0) && (ix2 <= numelt))
  {
    tmp.d = h[ix1].d;
    tmp.i = h[ix1].i;
    h[ix1].d = h[ix2].d;
    h[ix1].i = h[ix2].i;
    h[ix2].d = tmp.d;
    h[ix2].i = tmp.i;
    if(useAux)
    {
      tmp.aux[0] = h[ix1].aux[0];
      tmp.aux[1] = h[ix1].aux[1];
      tmp.aux[2] = h[ix1].aux[2];
      h[ix1].aux[0] = h[ix2].aux[0];
      h[ix1].aux[1] = h[ix2].aux[1];
      h[ix1].aux[2] = h[ix2].aux[2];
      h[ix2].aux[0] = tmp.aux[0];
      h[ix2].aux[1] = tmp.aux[1];
      h[ix2].aux[2] = tmp.aux[2];
    }
  }
  else
  {
    cout << "Trying to flip out-of-bounds elements in heap. Aborting ..." << endl;
    abort();
  }
}

void Heap::addToHeap(const int ii, const double dd, const double *cp)
{ // min-heap: minimum element at top
  // Step 0: add more space to heap if needed
  if(numelt == capacity)
  {
    capacity += blocksize;
    h.resize(capacity);
  }
  // Step 1: add element (i,d) to bottom of heap
  h[numelt].d = dd;
  h[numelt].i = ii;
  if(useAux)
  {
    h[numelt].aux[0] = cp[0];
    h[numelt].aux[1] = cp[1];
    h[numelt].aux[2] = cp[2];
  }
  
  // Step 2: upheap
  int done = 0;
  int pos = numelt;
  while(!done && pos > 0)
    if(h[pos].d < h[par(pos)].d)
    {
      flip(pos,par(pos));
      pos = par(pos);
    }
    else
      done = 1;
  // Step 3: increment numelt
  numelt++;
}

struct helt Heap::popFromHeap()
{ // pop top element
  struct helt lval;
  // Step 1: Store return value
  if(numelt == 0)
  {
    lval.i = HEAPDONE;
    lval.d = 0.0f;
    if(useAux)
    {
      lval.aux[0] = 0.0f;
      lval.aux[1] = 0.0f;
      lval.aux[2] = 0.0f;
    }
    return(lval);
  }
  lval.d = h[0].d;
  lval.i = h[0].i;
  if(useAux)
  {
    lval.aux[0] = h[0].aux[0];
    lval.aux[1] = h[0].aux[1];
    lval.aux[2] = h[0].aux[2];
  }
  numelt--;

  h[0].d = h[numelt].d;
  h[0].i = h[numelt].i;
  if(useAux)
  {
    h[0].aux[0] = h[numelt].aux[0];
    h[0].aux[1] = h[numelt].aux[1];
    h[0].aux[2] = h[numelt].aux[2];
  }
  int done = 0;
  int pos = 0;
  while(!done)
  {
    if(rch(pos) < numelt)
    {
      if(h[lch(pos)].d <= h[rch(pos)].d)
      {
	if(h[lch(pos)].d < h[pos].d)
	{
	  flip(lch(pos),pos);
	  pos = lch(pos);
	}
	else
	  done = 1;
      }
      else
      {
	if(h[rch(pos)].d < h[pos].d)
	{
	  flip(rch(pos),pos);
	  pos = rch(pos);
	}
	else
	  done = 1;
      }
    }
    else
    {
      if(lch(pos) < numelt)
      {
	if(h[lch(pos)].d < h[pos].d)
	{
	  flip(lch(pos),pos);
	  pos = lch(pos);
	}
	else
	  done = 1;
      }
      else
	done = 1;
    }
  }
  return(lval);
}

void Heap::showHeap() const
{
  for(int ii=0; ii<numelt; ++ii)
  {
    printf(" i=%d, d = %f\n",h[ii].i,h[ii].d);
  }
}
