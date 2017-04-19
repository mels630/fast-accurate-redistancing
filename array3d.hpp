#ifndef _ARRAY3D_HPP_
#define _ARRAY3D_HPP_

#include "defs.h"

#include<vector>
#include<cstdlib>
#include<iostream>
#include<cmath>
#include<cassert>
#include<algorithm>

using idx_t = std::size_t;

/// Array3D: template class for three-dimensional arrays with fixed spacing.
template <typename T> 
class Array3D
{
protected:
  std::vector<T> data; ///< contains (m x n x p) = N elements of type T
  const idx_t m;       ///< number of elements in the y-direction
  const idx_t n;       ///< number of elements in the x-direction
  const idx_t k;       ///< number of elements in the z-direction
  const idx_t N;       ///< total number of elements
  const double dx;     ///< grid spacing in x
  const double dy;     ///< grid spacing in y
  const double dz;     ///< grid spacing in z

  idx_t sub2ind(const idx_t ii, const idx_t jj, const idx_t kk) const;

public: 
  Array3D<T>() = delete; // don't use empty constructor
  explicit Array3D<T>(const idx_t nn);
  Array3D<T>(const idx_t mm, const idx_t nn, const idx_t kk);
  Array3D<T>(const T* indata, const idx_t mm, const idx_t nn, const idx_t kk);
  Array3D<T>(const idx_t mm, const idx_t nn, const idx_t kk, const double _dx, const double _dy, const double _dz);
  Array3D<T>(const T* indata, const idx_t mm, const idx_t nn, const idx_t kk, const double _dx, const double _dy, const double _dz);
  Array3D<T>(const Array3D<T> &Input, const int CopyFlag);
  Array3D<T>(const Array3D<T> &Input); 
  T get(const idx_t idx) const;
  T get(const idx_t ii, const idx_t jj, const idx_t kk) const;
  idx_t getm() const;
  idx_t getn() const;
  idx_t getk() const;
  idx_t getN() const;
  double getdx() const;
  double getdy() const;
  double getdz() const;
  double lenx() const;
  double leny() const;
  double lenz() const;
  idx_t size() const;
  void put(const T value, const idx_t idx);
  void put(const T value, const idx_t ii, const idx_t jj, const idx_t kk);
  void putadd(const T value, const idx_t idx);
  void putadd(const T value, const idx_t ii, const idx_t jj, const idx_t kk);

  T maxval() const;
  T minval() const;
  T minabsval() const;
  T maxabsval() const;
  void fillWithValue(const T value);
  T sixNborMin(const idx_t idx) const;
  T sixNborMax(const idx_t idx) const;
  std::vector<T> const &returnData() const;
  std::vector<T>& returnData();
  double getX(const idx_t idx) const;
  double getY(const idx_t idx) const;
  double getZ(const idx_t idx) const;
  idx_t getXidx(const idx_t idx) const;
  idx_t getYidx(const idx_t idx) const;
  idx_t getZidx(const idx_t idx) const;
  idx_t xp(const idx_t idx) const;
  idx_t xm(const idx_t idx) const;
  idx_t yp(const idx_t idx) const;
  idx_t ym(const idx_t idx) const;
  idx_t zp(const idx_t idx) const;
  idx_t zm(const idx_t idx) const;
  T getxp(const idx_t idx) const;
  T getxm(const idx_t idx) const;
  T getyp(const idx_t idx) const;
  T getym(const idx_t idx) const;
  T getzp(const idx_t idx) const;
  T getzm(const idx_t idx) const;

  idx_t countVal(const T value) const;

  // algebraic operations
  T& operator[] ( idx_t I );
  const T& operator[] ( idx_t I ) const;
  Array3D<T>& operator= ( const Array3D<T>& Vec );
  Array3D<T>& operator*= ( const T Value );
  Array3D<T>& operator/= ( const T Value );
  Array3D<T>& operator+= ( const Array3D<T> &Vec );
  Array3D<T>& operator-= ( const Array3D<T> &Vec );
};

/// Full constructor
/// \param[in] indata : Pointer to array of input data (to copy). Assumed to contain (mm x nn x kk) elements.
/// \param[in] mm     : Number of elements per row
/// \param[in] nn     : Number of elements per column  
/// \param[in] kk     : Number of elements per slice
/// \param[in] _dx    : Element spacing in x
/// \param[in] _dy    : Element spacing in y
/// \param[in] _dz    : Element spacing in z
template <typename T>
Array3D<T>::Array3D(const T* indata, const idx_t mm, const idx_t nn, const idx_t kk, const double _dx, const double _dy, const double _dz) :
  data(),
  m(mm),
  n(nn),
  k(kk),
  N(mm*nn*kk),
  dx(_dx),
  dy(_dy),
  dz(_dz)
{ // assume indata contains at least (exactly) N entries
  assert( (m>0) && (n>0) && (k>0) && (dx>0.) && (dy>0.) && (dz>0.) );
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
/// \param[in] mm     : Number of elements per row
/// \param[in] nn     : Number of elements per column  
/// \param[in] kk     : Number of elements per slice
/// \param[in] _dx    : Element spacing in x
/// \param[in] _dy    : Element spacing in y
/// \param[in] _dz    : Element spacing in z
template <typename T> 
Array3D<T>::Array3D(const idx_t mm, const idx_t nn, const idx_t kk, const double _dx, const double _dy, const double _dz) :
  Array3D(nullptr, mm, nn, kk, _dx, _dy, _dz)
{ }

/// Constructor
/// \param[in] mm     : Number of elements per row
/// \param[in] nn     : Number of elements per column  
/// \param[in] kk     : Number of elements per slice
template <typename T> 
Array3D<T>::Array3D(const idx_t mm, const idx_t nn, const idx_t kk) :
  Array3D(mm, nn, kk, 1./static_cast<double>(nn), 1./static_cast<double>(mm), 1./static_cast<double>(kk))
{ } 

/// Constructor
/// \param[in] nn     : Number of elements in each dimension
template <typename T> 
Array3D<T>::Array3D(const idx_t nn) :
  Array3D(nn, nn, nn)
{ }

/// Constructor
/// \param[in] indata : Pointer to array of input data (to copy). Assumed to contain (mm x nn x kk) elements.
/// \param[in] mm     : Number of elements per row
/// \param[in] nn     : Number of elements per column  
/// \param[in] kk     : Number of elements per slice
template <typename T>
Array3D<T>::Array3D(const T* indata, const idx_t mm, const idx_t nn, const idx_t kk) :
  Array3D(indata, mm, nn, kk, 1./static_cast<double>(nn), 1./static_cast<double>(mm), 1./static_cast<double>(kk))
{ }

/// Copy constructor, deep or structure
/// \param[in] Input    : Array2D<T> object to copy
/// \param[in] CopyFlag : Deep copy if 0, structure copy otherwise
template <typename T>
Array3D<T>::Array3D( const Array3D<T> &Input, const int CopyFlag) :
  Array3D(Input.m, Input.n, Input.k, Input.dx, Input.dy, Input.dz)
{
  if (CopyFlag == 0) // deep copy
    data = Input.data;
  // else if (CopyFlag == 1) // structure copy, done by call to Array2D constructor
}

/// Copy constructor, deep or structure
/// \param[in] Input    : Array2D<T> object to copy
template <typename T>
Array3D<T>::Array3D( const Array3D<T> &Input) :
  Array3D(Input, 0) // default to deep copy
{
  data = Input.data;
}

/// Convert (ii,jj,kk) row-column-slice index to flat index
/// \param[in] ii : Row index
/// \param[in] jj : Column index
/// \param[in] kk : Slice index
/// \return         Flat index
template <typename T>
idx_t Array3D<T>::sub2ind(const idx_t ii, const idx_t jj, const idx_t kk) const
{
  assert((ii < m) && (jj < n) && (kk < k));
  return(ii+m*(jj+n*kk));
}

/// Flat-index getter
/// \param[in] idx : Flat index to get array value at
/// \return          Value of array at idx
template <typename T>
T Array3D<T>::get(const idx_t idx) const
{
  return(data[idx]);
}

/// 3D-index getter
/// \param[in] ii : Row index to get array value at
/// \param[in] jj : Column index to get array value at
/// \param[in] kk : Slice index to get array value at
/// \return         Value of array at location
template <typename T>
T Array3D<T>::get(const idx_t ii, const idx_t jj, const idx_t kk) const
{
  return(data[sub2ind(ii,jj,kk)]);
}

/// Return number of rows in the array
/// \return Number of rows in the array
template <typename T>
idx_t Array3D<T>::getm() const
{
  return(m);
}

/// Return number of columns in the array
/// \return Number of columns in the array
template <typename T>
idx_t Array3D<T>::getn() const
{
  return(n);
}

/// Return number of slices in the array
/// \return Number of slices in the array
template <typename T>
idx_t Array3D<T>::getk() const
{
  return(k);
}

/// Return total number of elements in the array
/// \return Total number of elements in the array
template <typename T>
idx_t Array3D<T>::getN() const
{
  return(N);
}

/// Return pixel x spacing
/// \return Pixel x spacing
template <typename T>
double Array3D<T>::getdx() const
{
  return(dx);
}

/// Return pixel y spacing
/// \return Pixel y spacing
template <typename T>
double Array3D<T>::getdy() const
{
  return(dy);
}

/// Return pixel z spacing
/// \return Pixel z spacing
template <typename T>
double Array3D<T>::getdz() const
{
  return(dz);
}

/// Return total length of the array in the x direction (number of elements * pixel x spacing)
/// \return Total length of array in the x direction
template <typename T>
double Array3D<T>::lenx() const
{
  return(static_cast<double>(n) * dx);
}

/// Return total length of the array in the y direction (number of elements * pixel y spacing)
/// \return Total length of array in the y direction
template <typename T>
double Array3D<T>::leny() const
{
  return(static_cast<double>(m) * dy);
}

/// Return total length of the array in the z direction (number of elements * pixel z spacing)
/// \return Total length of array in the z direction
template <typename T>
double Array3D<T>::lenz() const
{
  return(static_cast<double>(k) * dz);
}

/// Return total number of elements in the array
/// \return Total number of elements in the array
template <typename T>
idx_t Array3D<T>::size() const
{
  return(N);
}

/// Flat-index putter
/// \param[in] value : Value to place in the array
/// \param[in] idx   : Flat index to set array value at
template <typename T>
void Array3D<T>::put(const T value, const idx_t idx)
{
  data[idx] = value;
}

/// 2D-index putter
/// \param[in] value : Value to place in the array
/// \param[in] ii    : Row index to set array value at
/// \param[in] jj    : Column index to set array value at
/// \param[in] kk    : Slice index to set array value at
template <typename T>
void Array3D<T>::put(const T value, const idx_t ii, const idx_t jj, const idx_t kk)
{
  data[sub2ind(ii,jj,kk)] = value;
}

/// Flat-index incrementer
/// \param[in] value : Value to add to location in the array
/// \param[in] idx   : Flat index to increment array value at
template <typename T>
void Array3D<T>::putadd(const T value, const idx_t idx)
{
  data[idx] += value;
}

/// 2D-index incrementer
/// \param[in] value : Value to add to location in the array
/// \param[in] ii    : Row index to increment array value at
/// \param[in] jj    : Column index to increment array value at
/// \param[in] kk    : Slice index to increment array value at
template <typename T>
void Array3D<T>::putadd(const T value, const idx_t ii, const idx_t jj, const idx_t kk)
{
  data[sub2ind(ii,jj,kk)] += value;
}

/// Maximum value getter
/// \return Maximum value in the array
template <typename T>
T Array3D<T>::maxval() const
{
  assert(N>0);
  return *std::max_element(data.begin(), data.end());
}

/// Minimum value getter
/// \return Minimum value in the array
template <typename T>
T Array3D<T>::minval() const
{
  assert(N>0);
  return *std::min_element(data.begin(), data.end());
}

/// Minimum absolute value getter
/// \return Minimum absolute value in the array
template <typename T>
T Array3D<T>::minabsval() const
{
  assert(N>0);
  return *std::min_element(data.begin(), data.end(), [](T const &v1, T const &v2)->bool { return std::abs(v1) < std::abs(v2); });
}

/// Minimum absolute value getter
/// \return Maximum absolute value in the array
template <typename T>
T Array3D<T>::maxabsval() const
{
  assert(N>0);
  return *std::max_element(data.begin(), data.end(), [](T const &v1, T const &v2)->bool { return std::abs(v1) < std::abs(v2); });
}

/// Fill the array with value
/// \param[in] value : Value to fill array with
template <typename T>
void Array3D<T>::fillWithValue(const T value)
{
  std::fill(data.begin(), data.end(), value);
}

/// Return the value of the minimum of the six neighbors of flat index idx
/// \param[in] idx : Flat index to compute six-neighbor minimum at
template <typename T>
T Array3D<T>::sixNborMin(const idx_t idx) const
{
  return(mymin(mymin(mymin(data[xp(idx)],data[xm(idx)]),mymin(data[yp(idx)],data[ym(idx)])),mymin(data[zp(idx)],data[zm(idx)])));
}

/// Return the value of the maximum of the six neighbors of flat index idx
/// \param[in] idx : Flat index to compute six-neighbor maximum at
template <typename T>
T Array3D<T>::sixNborMax(const idx_t idx) const
{
  return(mymax(mymax(mymax(data[xp(idx)],data[xm(idx)]),mymax(data[yp(idx)],data[ym(idx)])),mymax(data[zp(idx)],data[zm(idx)])));
}
 
/// Get const reference to vector of data in array
/// \return Const reference to vector of data in array
template <typename T>
std::vector<T> const &Array3D<T>::returnData() const
{ // returns a const ref to data
  return(data);
}

/// Get reference to vector of data in array
/// \return Reference to vector of data in array
template <typename T>
std::vector<T> &Array3D<T>::returnData()
{ // returns a ref to data
  return(data);
}

/// Get physical x location of flat index in array (assuming origin at (0,0,0) coordinate)
/// \param[in] idx : Flat index
/// \return          Physical x location
template <typename T>
double Array3D<T>::getX(const idx_t idx) const
{
  return(static_cast<double>(getXidx(idx)) / static_cast<double>(n));
}

/// Get physical y location of flat index in array (assuming origin at (0,0,0) coordinate)
/// \param[in] idx : Flat index
/// \return          Physical y location
template <typename T>
double Array3D<T>::getY(const idx_t idx) const
{
  return(static_cast<double>(getYidx(idx)) / static_cast<double>(m));
}

/// Get physical z location of flat index in array (assuming origin at (0,0,0) coordinate)
/// \param[in] idx : Flat index
/// \return          Physical z location
template <typename T>
double Array3D<T>::getZ(const idx_t idx) const
{
  return(static_cast<double>(getZidx(idx)) / static_cast<double>(k));
}

/// Get x-index from flat index
/// \param[in] idx : Flat index
/// \return          X-index
template <typename T>
idx_t Array3D<T>::getXidx(const idx_t idx) const
{
  return((idx%(m*n))/m);
}

/// Get y-index from flat index
/// \param[in] idx : Flat index
/// \return          Y-index
template <typename T>
idx_t Array3D<T>::getYidx(const idx_t idx) const
{
  return((idx%(m*n))%m);
}

/// Get z-index from flat index
/// \param[in] idx : Flat index
/// \return          Z-index
template <typename T>
idx_t Array3D<T>::getZidx(const idx_t idx) const
{
  return(idx/(m*n));
}

/// Return flat index of adjacent location in positive x direction
/// \param[in] idx : Current flat index
/// \return           Flat index of adjacent location in positive x direction
template <typename T>
idx_t Array3D<T>::xp(const idx_t idx) const
{
  return( (idx%(m*n))<(n*m-m) ? (idx +m) : (idx+m-m*n) );
}

/// Return flat index of adjacent location in negative x direction
/// \param[in] idx : Current flat index
/// \return           Flat index of adjacent location in negative x direction
template <typename T>
idx_t Array3D<T>::xm(const idx_t idx) const
{
  return( (idx%(m*n))>(m-1) ? (idx -m) : (idx-m+m*n) );
}

/// Return flat index of adjacent location in positive y direction
/// \param[in] idx : Current flat index
/// \return           Flat index of adjacent location in positive y direction
template <typename T>
idx_t Array3D<T>::yp(const idx_t idx) const
{
  return( ((idx%(m*n))%m)==(m-1) ? (idx-m+1) : (idx+1) );
}

/// Return flat index of adjacent location in negative y direction
/// \param[in] idx : Current flat index
/// \return           Flat index of adjacent location in negative y direction
template <typename T>
idx_t Array3D<T>::ym(const idx_t idx) const
{
  return( ((idx%(m*n))%m)==(0) ? (idx+m-1) : (idx-1) );
}

/// Return flat index of adjacent location in positive z direction
/// \param[in] idx : Current flat index
/// \return           Flat index of adjacent location in positive z direction
template <typename T>
idx_t Array3D<T>::zp(const idx_t idx) const
{
  return( (idx < (m*n*(k-1))) ? (idx+m*n) : (idx+m*n-m*n*k) );
}

/// Return flat index of adjacent location in negative z direction
/// \param[in] idx : Current flat index
/// \return           Flat index of adjacent location in negative z direction
template <typename T>
idx_t Array3D<T>::zm(const idx_t idx) const
{
  return( (idx > (m*n-1)) ? (idx-m*n) : (idx-m*n+m*n*k) );
}

/// Return array value at adjacent location in positive x direction
/// \param[in] idx : Current flat index
/// \return           Array value at adjacent location in positive x direction
template <typename T>
T Array3D<T>::getxp(const idx_t idx) const
{
  return(data[xp(idx)]);
}

/// Return array value at adjacent location in negative x direction
/// \param[in] idx : Current flat index
/// \return           Array value at adjacent location in negative x direction
template <typename T>
T Array3D<T>::getxm(const idx_t idx) const
{
  return(data[xm(idx)]);
}

/// Return array value at adjacent location in positive y direction
/// \param[in] idx : Current flat index
/// \return           Array value at adjacent location in positive y direction
template <typename T>
T Array3D<T>::getyp(const idx_t idx) const
{
  return(data[yp(idx)]);
}

/// Return array value at adjacent location in negative y direction
/// \param[in] idx : Current flat index
/// \return           Array value at adjacent location in negative y direction
template <typename T>
T Array3D<T>::getym(const idx_t idx) const
{
  return(data[ym(idx)]);
}

/// Return array value at adjacent location in positive z direction
/// \param[in] idx : Current flat index
/// \return           Array value at adjacent location in positive z direction
template <typename T>
T Array3D<T>::getzp(const idx_t idx) const
{
  return(data[zp(idx)]);
}

/// Return array value at adjacent location in negative z direction
/// \param[in] idx : Current flat index
/// \return           Array value at adjacent location in negative z direction
template <typename T>
T Array3D<T>::getzm(const idx_t idx) const
{
  return(data[zm(idx)]);
}

/// Return number of elements in array matching value
/// \param[in] value : Value to match
/// \return            Count of matching values
template <typename T>
idx_t Array3D<T>::countVal(const T value) const
{
  return std::count_if(data.begin(), data.end(), [&value](T const &t)->bool{return value == t;});
}

/// Flat index access operator
/// \param[in] I: Flat index
/// \return       Reference to element at this flat index
template <typename T>
T& Array3D<T>::operator[] ( idx_t I )
{
  return data[I];
}

/// Flat index access const operator
/// \param[in] I: Flat index
/// \return       Const reference to element at this flat index
template <typename T>
const T& Array3D<T>::operator[] ( idx_t I ) const
{
  return(data[I]);
}

/// Assignment operator
/// \param[in] Vec : Array to assign from
/// \return          Reference to the assigned object
template <typename T>
Array3D<T>& Array3D<T>::operator= ( const Array3D<T>& Vec )
{
  assert(m == Vec.getm());
  assert(n == Vec.getn());
  assert(k == Vec.getk());
  data = Vec.data;
  return(*this);
}

/// In-place constant multiplication operator
/// \param[in] Value : Value to mutiply each element of the array by
/// \return            Reference to this
template <typename T>
Array3D<T>& Array3D<T>::operator*= ( const T Value )
{
  std::transform(data.begin(), data.end(), data.begin(),
                 [&Value](T const &t)->T { return Value*t; });
  return(*this);
}

/// In-place constant division operator
/// \param[in] Value : Value to divide each element of the array by
/// \return            Reference to this
template <typename T>
Array3D<T>& Array3D<T>::operator/= ( const T Value )
{
  assert(Value != T(0));
  std::transform(data.begin(), data.end(), data.begin(),
                 [&Value](T const &t)->T { return t/Value; });
  return(*this);
}

/// Elementwise in-place addition operator
/// \param[in] Vec : Array to add to this array
/// \return          Reference to this
template <typename T>
Array3D<T>& Array3D<T>::operator+= ( const Array3D<T> &Vec )
{
  assert(Vec.getN() == getN());
  std::transform(data.begin(), data.end(), Vec.returnData().begin(), data.begin(),
                 [](T const &inp1, T const &inp2)->T { return inp1+inp2; });
  return(*this);
}

/// Elementwise in-place subtraction operator
/// \param[in] Vec : Array to subtract from this array
/// \return          Reference to this
template <typename T>
Array3D<T>& Array3D<T>::operator-= ( const Array3D<T> &Vec )
{
  assert(Vec.getN() == getN());
  std::transform(data.begin(), data.end(), Vec.returnData().begin(), data.begin(),
                 [](T const &inp1, T const &inp2)->T { return inp1-inp2; });
  return(*this);
}

#endif
