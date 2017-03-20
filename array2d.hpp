#ifndef _ARRAY2D_HPP_
#define _ARRAY2D_HPP_

// To-do:
// [1] Doxygen

#include<vector>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<cassert>
#include<cstring>
#include<algorithm>

#include "defs.h"

using idx_t = std::size_t;

/// Array2D: template class for two-dimensional arrays with fixed spacing.
template <typename T> 
class Array2D
{

protected:
  std::vector<T> data; ///< contains (m x n) = N elements of type T
  idx_t const m;       ///< number of elements in the y-direction
  idx_t const n;       ///< number of elements in the x-direction
  idx_t const N;       ///< total number of elements
  double const dx;     ///< grid spacing in x
  double const dy;     ///< grid spacing in y

  idx_t sub2ind(idx_t const ii, idx_t const jj) const;

public: 
  Array2D<T>() = delete; // don't use empty constructor
  explicit Array2D<T>(idx_t const nn);
  Array2D<T>(idx_t const mm, idx_t const nn);
  Array2D<T>(T* const indata, idx_t const mm, idx_t const nn);
  Array2D<T>(idx_t const mm, idx_t const nn, double const _dx, double const _dy);
  Array2D<T>(T* const indata, idx_t const mm, idx_t const nn, double const _dx, double const _dy);
  Array2D<T>(Array2D<T> const &Input, int const CopyFlag);
  Array2D<T>(Array2D<T> const &Input); 
  T get(idx_t const idx) const;
  T get(idx_t const ii, idx_t const jj) const;
  idx_t getm() const;
  idx_t getn() const;
  idx_t getN() const;
  double getdx() const;
  double getdy() const;
  double lenx() const;
  double leny() const;
  idx_t size() const;
  void put(T const value, idx_t const idx);
  void put(T const value, idx_t const ii, idx_t const jj);
  void putadd(T const value, idx_t const idx);
  void putadd(T const value, idx_t const ii, idx_t const jj);

  T maxval() const;
  T minval() const;
  T minabsval() const;
  void fillWithValue(T const value);
  T fourNborMin(idx_t const idx) const;
  T fourNborMax(idx_t const idx) const;
  std::vector<T> const& returnData() const;
  double getX(idx_t const idx) const;
  double getY(idx_t const idx) const;
  idx_t getXidx(idx_t const idx) const;
  idx_t getYidx(idx_t const idx) const;
  idx_t xp(idx_t const idx) const;
  idx_t xm(idx_t const idx) const;
  idx_t yp(idx_t const idx) const;
  idx_t ym(idx_t const idx) const;
  T getxp(idx_t const idx) const;
  T getxm(idx_t const idx) const;
  T getyp(idx_t const idx) const;
  T getym(idx_t const idx) const;

  idx_t countVal(T const value) const;
  
  // algebraic operations
  T& operator[] ( idx_t I );
  T const& operator[] ( idx_t I ) const;
  Array2D<T>& operator= ( const Array2D<T>& Vec );
  Array2D<T>& operator*= ( T const Value );
  Array2D<T>& operator/= ( T const Value );
  Array2D<T>& operator+= ( const Array2D<T> &Vec );
  Array2D<T>& operator-= ( const Array2D<T> &Vec );
  void addMultiple ( const Array2D<T>& Vec, T Factor );
};

/// Full constructor
/// \param[in] indata : Pointer to array of input data (to copy). Assumed to contain (mm x nn) elements.
/// \param[in] mm     : Number of elements per row
/// \param[in] nn     : Number of elements per column  
/// \param[in] _dx    : Element spacing in x
/// \param[in] _dy    : Element spacing in y
template <typename T>
Array2D<T>::Array2D(T* const indata, idx_t const mm, idx_t const nn, double const _dx, double const _dy) :
  m(mm),
  n(nn),
  N(mm*nn),
  dx(_dx),
  dy(_dy)
{
  if(N <= 0)
  {
    std::cout << "N = " << N << ", should be positive. Aborting ..." << std::endl;
    abort();
  }
  if (indata != nullptr)
    data = std::vector<T>(indata, indata+N);
  else
    data.assign(N, T(0));
}

/// Constructor
/// \param[in] nn : Number of elements per row and per column
template <typename T> 
Array2D<T>::Array2D(idx_t const nn) :
  Array2D(nullptr, nn, nn)
{ }

/// Constructor
/// \param[in] mm : Number of elements per row
/// \param[in] nn : Number of elements per column
template <typename T> 
Array2D<T>::Array2D(idx_t const mm, idx_t const nn) :
  Array2D(nullptr, mm, nn)
{ }

/// Constructor
/// \param[in] indata : Pointer to array of input data (to copy). Assumed to contain (mm x nn) elements.
/// \param[in] mm     : Number of elements per row
/// \param[in] nn     : Number of elements per column  
template <typename T>
Array2D<T>::Array2D(T* const indata, idx_t const mm, idx_t const nn) :
  Array2D(indata, mm, nn, 1./static_cast<double>(nn), 1./static_cast<double>(mm))
{ }

/// Constructor
/// \param[in] mm     : Number of elements per row
/// \param[in] nn     : Number of elements per column  
/// \param[in] _dx    : Element spacing in x
/// \param[in] _dy    : Element spacing in y
template <typename T> 
Array2D<T>::Array2D(idx_t const mm, idx_t const nn, double const _dx, double const _dy) :
  Array2D(nullptr, mm, nn, _dx, _dy)
{ }

/// Copy constructor, deep or structure
/// \param[in] Input    : Array2D<T> object to copy
/// \param[in] CopyFlag : Deep copy if 0, structure copy otherwise
template <typename T>
Array2D<T>::Array2D(Array2D<T> const &Input, int const CopyFlag) :
  Array2D(Input.m, Input.n, Input.dx, Input.dy)
{
  if (CopyFlag == 0) // deep copy
    data = Input.data;
  // else if (CopyFlag == 1) // structure copy, done by call to Array2D constructor
}

/// Copy constructor
/// \param[in] Input    : Array2D<T> object to copy
template <typename T>
Array2D<T>::Array2D(Array2D<T> const &Input) :
  Array2D(Input, 0) // default to deep copy
{ } 

/// Convert (ii,jj) row-column index to flat index.
/// \param[in] ii : row index
/// \param[in] jj : column index
/// \return         flat index
template <typename T>
idx_t Array2D<T>::sub2ind(idx_t const ii, idx_t const jj) const
{
  assert((ii < m) && (jj < n));
  return(ii+m*jj);
}

template <typename T>
idx_t Array2D<T>::xp(idx_t const idx) const
{
  return((idx)<(N-m) ? (idx +m) : (idx+m-N));
}

template <typename T>
idx_t Array2D<T>::xm(idx_t const idx) const
{
  return((idx)>(m-1) ? (idx -m) : (idx-m+N));
}

template <typename T>
idx_t Array2D<T>::yp(idx_t const idx) const
{
  return( (idx%m)==(m-1) ? (idx-m+1) : (idx+1));
}

template <typename T>
idx_t Array2D<T>::ym(idx_t const idx) const
{
  return( (idx%m)==(0) ? (idx+m-1) : (idx-1));
}

template <typename T>
T Array2D<T>::getxp(idx_t const idx) const
{
  return(data[xp(idx)]);
}

template <typename T>
T Array2D<T>::getxm(idx_t const idx) const
{
  return(data[xm(idx)]);
}

template <typename T>
T Array2D<T>::getyp(idx_t const idx) const
{
  return(data[yp(idx)]);
}

template <typename T>
T Array2D<T>::getym(idx_t const idx) const
{
  return(data[ym(idx)]);
}

template <typename T>
idx_t Array2D<T>::getm() const
{
  return(m);
}

template <typename T>
idx_t Array2D<T>::getn() const
{
  return(n);
}

template <typename T>
idx_t Array2D<T>::getN() const
{
  return(N);
}

template <typename T>
double Array2D<T>::getdx() const
{
  return(dx);
}

template <typename T>
double Array2D<T>::getdy() const
{
  return(dy);
}

template <typename T>
double Array2D<T>::lenx() const
{
  return(static_cast<double>(n) * dx);
}

template <typename T>
double Array2D<T>::leny() const
{
  return(static_cast<double>(m) * dy);
}

template <typename T>
idx_t Array2D<T>::size() const
{
  return(N);
}

/// Flat-index getter
/// \param[in] idx : Flat index to get array value at
/// \return          Value of array at idx
template <typename T>
T Array2D<T>::get(idx_t const idx) const
{
  return(data[idx]);
}

/// 2D-index getter
/// \param[in] ii : Row index to get array value at
/// \param[in] jj : Column index to get array value at
/// \return         Value of array at idx
template <typename T>
T Array2D<T>::get(idx_t const ii, idx_t const jj) const
{
  return(data[sub2ind(ii,jj)]);
}

template <typename T>
void Array2D<T>::put(T const value, idx_t const idx)
{
  data[idx] = value;
}

template <typename T>
void Array2D<T>::put(T const value, idx_t const ii, idx_t const jj)
{
  data[sub2ind(ii,jj)] = value;
}

template <typename T>
void Array2D<T>::putadd(T const value, idx_t const idx)
{
  data[idx] += value;
}

template <typename T>
void Array2D<T>::putadd(T const value, idx_t const ii, idx_t const jj)
{
  data[sub2ind(ii,jj)] += value;
}

template <typename T>
T Array2D<T>::maxval() const
{
  assert(N>0);
  return *std::max_element(data.begin(), data.end());
}

template <typename T>
T Array2D<T>::minval() const
{
  assert(N>0);
  return *std::min_element(data.begin(), data.end());
}

template <typename T>
T Array2D<T>::minabsval() const
{
  assert(N>0);
  return *std::min_element(data.begin(), data.end(), [](T const &v1, T const &v2)->bool { return std::abs(v1) < std::abs(v2); });
}

template <typename T>
void Array2D<T>::fillWithValue(T const value)
{
  std::fill(data.begin(), data.end(), value);
}

template <typename T>
T Array2D<T>::fourNborMin(idx_t const idx) const
{
  return(std::min(std::min(data[xp(idx)],data[xm(idx)]),std::min(data[yp(idx)],data[ym(idx)])));
}

template <typename T>
T Array2D<T>::fourNborMax(idx_t const idx) const
{
  return(std::max(std::max(data[xp(idx)],data[xm(idx)]),std::max(data[yp(idx)],data[ym(idx)])));
}

template <typename T>
std::vector<T> const& Array2D<T>::returnData() const
{ // returns a const ref to data std::vector
  return(data);
}

template <typename T>
double Array2D<T>::getX(idx_t const idx) const
{
  return(static_cast<double>(getXidx(idx)) * dx);
}

template <typename T>
double Array2D<T>::getY(idx_t const idx) const
{
  return(static_cast<double>(getYidx(idx)) * dy);
}

template <typename T>
idx_t Array2D<T>::getXidx(idx_t const idx) const
{
  return(idx/m);
}

template <typename T>
idx_t Array2D<T>::getYidx(idx_t const idx) const
{
  return(idx%m);
}

template <typename T>
idx_t Array2D<T>::countVal(T const value) const
{
  return std::count_if(data.begin(), data.end(), [&value](T const &t)->bool{return value == t;});
}

template <typename T>
T& Array2D<T>::operator[] ( idx_t I )
{
  return data[I];
}

template <typename T>
T const& Array2D<T>::operator[] ( idx_t I ) const
{
  return data[I];
}

template <typename T>
Array2D<T>& Array2D<T>::operator= ( Array2D<T> const& Vec )
{
  assert(m == Vec.getm());
  assert(n == Vec.getn());
  data = Vec.data;
  return *this;
}

template <typename T>
Array2D<T>& Array2D<T>::operator*= ( T const Value )
{
  std::transform(data.begin(), data.end(), data.begin(),
                 [&Value](T const &t)->T { return Value*t; });
  return *this;
}

template <typename T>
Array2D<T>& Array2D<T>::operator/= ( T const Value )
{
  assert(Value != T(0));
  std::transform(data.begin(), data.end(), data.begin(),
                 [&Value](T const &t)->T { return t/Value; });
  return *this;
}

template <typename T>
Array2D<T>& Array2D<T>::operator+= ( const Array2D<T> &Vec )
{
  assert(Vec.getN() == getN());
  std::transform(data.begin(), data.end(), Vec.returnData().begin(), data.begin(),
                 [](T const &inp1, T const &inp2)->T { return inp1+inp2; });
  return *this;
}

template <typename T>
Array2D<T>& Array2D<T>::operator-= ( const Array2D<T> &Vec )
{
  assert(Vec.getN() == getN());
  std::transform(data.begin(), data.end(), Vec.returnData().begin(), data.begin(),
                 [](T const &inp1, T const &inp2)->T { return inp1-inp2; });
  return *this;
}

template <typename T>
void Array2D<T>::addMultiple ( const Array2D<T>& Vec, T Factor )
{
  assert(getN() == Vec.getN());
  std::transform(data.begin(), data.end(), Vec.returnData().begin(), data.begin(),
                 [&Factor](T const &t1, T const &t2)->T { return t1 + t2*Factor; });
}

#endif
