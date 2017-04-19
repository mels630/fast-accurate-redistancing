#ifndef _ARRAY2D_HPP_
#define _ARRAY2D_HPP_

#include "defs.h"

#include<vector>
#include<cstdlib>
#include<iostream>
#include<cmath>
#include<cassert>
#include<algorithm>

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
  std::vector<T> &returnData();
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
};

/// Full constructor
/// \param[in] indata : Pointer to array of input data (to copy). Assumed to contain (mm x nn) elements.
/// \param[in] mm     : Number of elements per row
/// \param[in] nn     : Number of elements per column  
/// \param[in] _dx    : Element spacing in x
/// \param[in] _dy    : Element spacing in y
template <typename T>
Array2D<T>::Array2D(T* const indata, idx_t const mm, idx_t const nn, double const _dx, double const _dy) :
  data(),
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

/// Convert (ii,jj) row-column index to flat index
/// \param[in] ii : Row index
/// \param[in] jj : Column index
/// \return         Flat index
template <typename T>
idx_t Array2D<T>::sub2ind(idx_t const ii, idx_t const jj) const
{
  assert((ii < m) && (jj < n));
  return(ii+m*jj);
}

/// Return flat index of adjacent location in positive x direction
/// \param[in] idx : Current flat index
/// \return           Flat index of adjacent location in positive x direction
template <typename T>
idx_t Array2D<T>::xp(idx_t const idx) const
{
  return((idx)<(N-m) ? (idx +m) : (idx+m-N));
}

/// Return flat index of adjacent location in negative x direction
/// \param[in] idx : Current flat index
/// \return           Flat index of adjacent location in negative x direction
template <typename T>
idx_t Array2D<T>::xm(idx_t const idx) const
{
  return((idx)>(m-1) ? (idx -m) : (idx-m+N));
}

/// Return flat index of adjacent location in positive y direction
/// \param[in] idx : Current flat index
/// \return           Flat index of adjacent location in positive y direction
template <typename T>
idx_t Array2D<T>::yp(idx_t const idx) const
{
  return( (idx%m)==(m-1) ? (idx-m+1) : (idx+1));
}

/// Return flat index of adjacent location in negative y direction
/// \param[in] idx : Current flat index
/// \return           Flat index of adjacent location in negative y direction
template <typename T>
idx_t Array2D<T>::ym(idx_t const idx) const
{
  return( (idx%m)==(0) ? (idx+m-1) : (idx-1));
}

/// Return array value at adjacent location in positive x direction
/// \param[in] idx : Current flat index
/// \return           Array value at adjacent location in positive x direction
template <typename T>
T Array2D<T>::getxp(idx_t const idx) const
{
  return(data[xp(idx)]);
}

/// Return array value at adjacent location in negative x direction
/// \param[in] idx : Current flat index
/// \return           Array value at adjacent location in negative x direction
template <typename T>
T Array2D<T>::getxm(idx_t const idx) const
{
  return(data[xm(idx)]);
}

/// Return array value at adjacent location in positive y direction
/// \param[in] idx : Current flat index
/// \return           Array value at adjacent location in positive y direction
template <typename T>
T Array2D<T>::getyp(idx_t const idx) const
{
  return(data[yp(idx)]);
}

/// Return array value at adjacent location in negative y direction
/// \param[in] idx : Current flat index
/// \return           Array value at adjacent location in negative y direction
template <typename T>
T Array2D<T>::getym(idx_t const idx) const
{
  return(data[ym(idx)]);
}

/// Return number of rows in the array
/// \return Number of rows in the array
template <typename T>
idx_t Array2D<T>::getm() const
{
  return(m);
}

/// Return number of columns in the array
/// \return Number of columns in the array
template <typename T>
idx_t Array2D<T>::getn() const
{
  return(n);
}

/// Return total number of elements in the array
/// \return Total number of elements in the array
template <typename T>
idx_t Array2D<T>::getN() const
{
  return(N);
}

/// Return pixel x spacing
/// \return Pixel x spacing
template <typename T>
double Array2D<T>::getdx() const
{
  return(dx);
}

/// Return pixel y spacing
/// \return Pixel y spacing
template <typename T>
double Array2D<T>::getdy() const
{
  return(dy);
}

/// Return total length of the array in the x direction (number of elements * pixel x spacing)
/// \return Total length of array in the x direction
template <typename T>
double Array2D<T>::lenx() const
{
  return(static_cast<double>(n) * dx);
}

/// Return total length of the array in the y direction (number of elements * pixel y spacing)
/// \return Total length of array in the y direction
template <typename T>
double Array2D<T>::leny() const
{
  return(static_cast<double>(m) * dy);
}

/// Return total number of elements in the array
/// \return Total number of elements in the array
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
/// \return         Value of array at location
template <typename T>
T Array2D<T>::get(idx_t const ii, idx_t const jj) const
{
  return(data[sub2ind(ii,jj)]);
}

/// Flat-index putter
/// \param[in] value : Value to place in the array
/// \param[in] idx   : Flat index to set array value at
template <typename T>
void Array2D<T>::put(T const value, idx_t const idx)
{
  data[idx] = value;
}

/// 2D-index putter
/// \param[in] value : Value to place in the array
/// \param[in] ii    : Row index to set array value at
/// \param[in] jj    : Column index to set array value at
template <typename T>
void Array2D<T>::put(T const value, idx_t const ii, idx_t const jj)
{
  data[sub2ind(ii,jj)] = value;
}

/// Flat-index incrementer
/// \param[in] value : Value to add to location in the array
/// \param[in] idx   : Flat index to increment array value at
template <typename T>
void Array2D<T>::putadd(T const value, idx_t const idx)
{
  data[idx] += value;
}

/// 2D-index incrementer
/// \param[in] value : Value to add to location in the array
/// \param[in] ii    : Row index to increment array value at
/// \param[in] jj    : Column index to increment array value at
template <typename T>
void Array2D<T>::putadd(T const value, idx_t const ii, idx_t const jj)
{
  data[sub2ind(ii,jj)] += value;
}

/// Maximum value getter
/// \return Maximum value in the array
template <typename T>
T Array2D<T>::maxval() const
{
  assert(N>0);
  return *std::max_element(data.begin(), data.end());
}

/// Minimum value getter
/// \return Minimum value in the array
template <typename T>
T Array2D<T>::minval() const
{
  assert(N>0);
  return *std::min_element(data.begin(), data.end());
}

/// Minimum absolute value getter
/// \return Minimum absolute value in the array
template <typename T>
T Array2D<T>::minabsval() const
{
  assert(N>0);
  return *std::min_element(data.begin(), data.end(), [](T const &v1, T const &v2)->bool { return std::abs(v1) < std::abs(v2); });
}

/// Fill the array with value
/// \param[in] value : Value to fill array with
template <typename T>
void Array2D<T>::fillWithValue(T const value)
{
  std::fill(data.begin(), data.end(), value);
}

/// Return the value of the minimum of the four neighbors of flat index idx
/// \param[in] idx : Flat index to compute four-neighbor minimum at
template <typename T>
T Array2D<T>::fourNborMin(idx_t const idx) const
{
  return(std::min(std::min(data[xp(idx)],data[xm(idx)]),std::min(data[yp(idx)],data[ym(idx)])));
}

/// Return the value of the maximum of the four neighbors of flat index idx
/// \param[in] idx : Flat index to compute four-neighbor maximum at
template <typename T>
T Array2D<T>::fourNborMax(idx_t const idx) const
{
  return(std::max(std::max(data[xp(idx)],data[xm(idx)]),std::max(data[yp(idx)],data[ym(idx)])));
}

/// Get const reference to vector of data in array
/// \return Const reference to vector of data in array
template <typename T>
std::vector<T> const& Array2D<T>::returnData() const
{ // returns a const ref to data std::vector
  return(data);
}

/// Get reference to vector of data in array
/// \return Reference to vector of data in array
template <typename T>
std::vector<T>& Array2D<T>::returnData()
{ // returns a ref to data std::vector
  return(data);
}

/// Get physical x location of flat index in array (assuming origin at (0,0) coordinate)
/// \param[in] idx : Flat index
/// \return          Physical x location
template <typename T>
double Array2D<T>::getX(idx_t const idx) const
{
  return(static_cast<double>(getXidx(idx)) * dx);
}

/// Get physical y location of flat index in array (assuming origin at (0,0) coordinate)
/// \param[in] idx : Flat index
/// \return Physical y location
template <typename T>
double Array2D<T>::getY(idx_t const idx) const
{
  return(static_cast<double>(getYidx(idx)) * dy);
}

/// Get x-index from flat index
/// \param[in] idx : Flat index
/// \return          X-index
template <typename T>
idx_t Array2D<T>::getXidx(idx_t const idx) const
{
  return(idx/m);
}

/// Get y-index from flat index
/// \param[in] idx : Flat index
/// \return          Y-index
template <typename T>
idx_t Array2D<T>::getYidx(idx_t const idx) const
{
  return(idx%m);
}

/// Return number of elements in array matching value
/// \param[in] value : Value to match
/// \return            Count of matching values
template <typename T>
idx_t Array2D<T>::countVal(T const value) const
{
  return std::count_if(data.begin(), data.end(), [&value](T const &t)->bool{return value == t;});
}

/// Flat index access operator
/// \param[in] I: Flat index
/// \return       Reference to element at this flat index
template <typename T>
T& Array2D<T>::operator[] ( idx_t I )
{
  return data[I];
}

/// Flat index access const operator
/// \param[in] I: Flat index
/// \return       Const reference to element at this flat index
template <typename T>
T const& Array2D<T>::operator[] ( idx_t I ) const
{
  return data[I];
}

/// Assignment operator
/// \param[in] Vec : Array to assign from
/// \return          Reference to the assigned object
template <typename T>
Array2D<T>& Array2D<T>::operator= ( Array2D<T> const& Vec )
{
  assert(m == Vec.getm());
  assert(n == Vec.getn());
  data = Vec.data;
  return *this;
}

/// In-place constant multiplication operator
/// \param[in] Value : Value to mutiply each element of the array by
/// \return            Reference to this
template <typename T>
Array2D<T>& Array2D<T>::operator*= ( T const Value )
{
  std::transform(data.begin(), data.end(), data.begin(),
                 [&Value](T const &t)->T { return Value*t; });
  return *this;
}

/// In-place constant division operator
/// \param[in] Value : Value to divide each element of the array by
/// \return            Reference to this
template <typename T>
Array2D<T>& Array2D<T>::operator/= ( T const Value )
{
  assert(Value != T(0));
  std::transform(data.begin(), data.end(), data.begin(),
                 [&Value](T const &t)->T { return t/Value; });
  return *this;
}

/// Elementwise in-place addition operator
/// \param[in] Vec : Array to add to this array
/// \return          Reference to this
template <typename T>
Array2D<T>& Array2D<T>::operator+= ( const Array2D<T> &Vec )
{
  assert(Vec.getN() == getN());
  std::transform(data.begin(), data.end(), Vec.returnData().begin(), data.begin(),
                 [](T const &inp1, T const &inp2)->T { return inp1+inp2; });
  return *this;
}

/// Elementwise in-place subtraction operator
/// \param[in] Vec : Array to subtract from this array
/// \return          Reference to this
template <typename T>
Array2D<T>& Array2D<T>::operator-= ( const Array2D<T> &Vec )
{
  assert(Vec.getN() == getN());
  std::transform(data.begin(), data.end(), Vec.returnData().begin(), data.begin(),
                 [](T const &inp1, T const &inp2)->T { return inp1-inp2; });
  return *this;
}

#endif
