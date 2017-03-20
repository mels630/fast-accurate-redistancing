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
#include<stack>
#include<cassert>
#include<cstring>
#include<algorithm>
#include<memory>

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

  inline idx_t sub2ind(idx_t const ii, idx_t const jj) const;
  void transposeData();

public: 
  Array2D<T>() = delete; // don't use empty constructor
  explicit Array2D<T>(idx_t const nn);
  Array2D<T>(idx_t const mm, idx_t const nn);
  Array2D<T>(T* const indata, idx_t const mm, idx_t const nn);
  Array2D<T>(idx_t const mm, idx_t const nn, double const _dx, double const _dy);
  Array2D<T>(T* const indata, idx_t const mm, idx_t const nn, double const _dx, double const _dy);
  Array2D<T>(Array2D<T> const &Input, int const CopyFlag);
  Array2D<T>(Array2D<T> const &Input); 
  inline T get(idx_t const idx) const;
  inline T get(idx_t const ii, idx_t const jj) const;
  inline idx_t getm() const;
  inline idx_t getn() const;
  inline idx_t getN() const;
  inline double getdx() const;
  inline double getdy() const;
  inline double lenx() const;
  inline double leny() const;
  inline idx_t size() const;
  inline void put(T const value, idx_t const idx);
  inline void put(T const value, idx_t const ii, idx_t const jj);
  inline void putadd(T const value, idx_t const idx);
  inline void putadd(T const value, idx_t const ii, idx_t const jj);
  void copyData(const Array2D<T> &input);
  inline T maxval() const;
  inline T minval() const;
  inline T minabsval() const;
  inline void fillWithValue(T const value);
  inline T fourNborMin(idx_t const idx) const;
  inline T fourNborMax(idx_t const idx) const;
  T* dataAddress();
  std::vector<T> const& returnData() const;
  std::shared_ptr<T> returnDataArray() const;
  Array2D<T> duplicateArray2D() const;
  inline double getX(idx_t const idx) const;
  inline double getY(idx_t const idx) const;
  inline idx_t getXidx(idx_t const idx) const;
  inline idx_t getYidx(idx_t const idx) const;
  inline idx_t xp(idx_t const idx) const;
  inline idx_t xm(idx_t const idx) const;
  inline idx_t yp(idx_t const idx) const;
  inline idx_t ym(idx_t const idx) const;
  inline T getxp(idx_t const idx) const;
  inline T getxm(idx_t const idx) const;
  inline T getyp(idx_t const idx) const;
  inline T getym(idx_t const idx) const;

  idx_t countVal(T const value) const;
  
  // algebraic operations
  inline T& operator[] ( idx_t I );
  inline T const& operator[] ( idx_t I ) const;
  inline Array2D<T>& operator= ( const Array2D<T>& Vec );
  inline Array2D<T>& operator*= ( T const Value );
  inline Array2D<T>& operator/= ( T const Value );
  inline Array2D<T>& operator+= ( const Array2D<T> &Vec );
  inline Array2D<T>& operator-= ( const Array2D<T> &Vec );
  inline void setAll( T const val );
  inline void addMultiple ( const Array2D<T>& Vec, T Factor );

  // load/save routines
  inline void save ( const char *FileName ) const;
  inline void load ( const char *FileName );
  inline void saveASCII ( const char *FileName ) const;
  inline void loadASCII ( const char *FileName );
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
inline idx_t Array2D<T>::sub2ind(idx_t const ii, idx_t const jj) const
{
  assert((ii < m) && (jj < n));
  return(ii+m*jj);
}

template <typename T>
inline idx_t Array2D<T>::xp(idx_t const idx) const
{
  return((idx)<(N-m) ? (idx +m) : (idx+m-N));
}

template <typename T>
inline idx_t Array2D<T>::xm(idx_t const idx) const
{
  return((idx)>(m-1) ? (idx -m) : (idx-m+N));
}

template <typename T>
inline idx_t Array2D<T>::yp(idx_t const idx) const
{
  return( (idx%m)==(m-1) ? (idx-m+1) : (idx+1));
}

template <typename T>
inline idx_t Array2D<T>::ym(idx_t const idx) const
{
  return( (idx%m)==(0) ? (idx+m-1) : (idx-1));
}

template <typename T>
inline T Array2D<T>::getxp(idx_t const idx) const
{
  return(data[xp(idx)]);
}

template <typename T>
inline T Array2D<T>::getxm(idx_t const idx) const
{
  return(data[xm(idx)]);
}

template <typename T>
inline T Array2D<T>::getyp(idx_t const idx) const
{
  return(data[yp(idx)]);
}

template <typename T>
inline T Array2D<T>::getym(idx_t const idx) const
{
  return(data[ym(idx)]);
}

template <typename T>
inline idx_t Array2D<T>::getm() const
{
  return(m);
}

template <typename T>
inline idx_t Array2D<T>::getn() const
{
  return(n);
}

template <typename T>
inline idx_t Array2D<T>::getN() const
{
  return(N);
}

template <typename T>
inline double Array2D<T>::getdx() const
{
  return(dx);
}

template <typename T>
inline double Array2D<T>::getdy() const
{
  return(dy);
}

template <typename T>
inline double Array2D<T>::lenx() const
{
  return(static_cast<double>(n) * dx);
}

template <typename T>
inline double Array2D<T>::leny() const
{
  return(static_cast<double>(m) * dy);
}

template <typename T>
inline idx_t Array2D<T>::size() const
{
  return(N);
}

/// Flat-index getter
/// \param[in] idx : Flat index to get array value at
/// \return          Value of array at idx
template <typename T>
inline T Array2D<T>::get(idx_t const idx) const
{
  return(data[idx]);
}

/// 2D-index getter
/// \param[in] ii : Row index to get array value at
/// \param[in] jj : Column index to get array value at
/// \return         Value of array at idx
template <typename T>
inline T Array2D<T>::get(idx_t const ii, idx_t const jj) const
{
  return(data[sub2ind(ii,jj)]);
}

template <typename T>
inline void Array2D<T>::put(T const value, idx_t const idx)
{
  data[idx] = value;
}

template <typename T>
inline void Array2D<T>::put(T const value, idx_t const ii, idx_t const jj)
{
  data[sub2ind(ii,jj)] = value;
}

template <typename T>
inline void Array2D<T>::putadd(T const value, idx_t const idx)
{
  data[idx] += value;
}

template <typename T>
inline void Array2D<T>::putadd(T const value, idx_t const ii, idx_t const jj)
{
  data[sub2ind(ii,jj)] += value;
}

template <typename T>
void Array2D<T>::copyData(const Array2D<T> &input)
{ // copy data after checking for agreement of sizes
  if((input.getm() != m) || (input.getn() != n))
    std::cout << "Data does not agree in size. Skipping call to copyData" << std::endl;
  else
    data = input.returnData();
}

template <typename T>
inline T Array2D<T>::maxval() const
{
  assert(N>0);
  return *std::max_element(data.begin(), data.end());
}

template <typename T>
inline T Array2D<T>::minval() const
{
  assert(N>0);
  return *std::min_element(data.begin(), data.end());
}

template <typename T>
inline T Array2D<T>::minabsval() const
{
  assert(N>0);
  return *std::min_element(data.begin(), data.end(), [](T const &v1, T const &v2)->bool { return std::abs(v1) < std::abs(v2); });
}

template <typename T>
inline void Array2D<T>::fillWithValue(T const value)
{
  std::fill(data.begin(), data.end(), value);
}

template <typename T>
inline T Array2D<T>::fourNborMin(idx_t const idx) const
{
  return(std::min(std::min(data[xp(idx)],data[xm(idx)]),std::min(data[yp(idx)],data[ym(idx)])));
}

template <typename T>
inline T Array2D<T>::fourNborMax(idx_t const idx) const
{
  return(std::max(std::max(data[xp(idx)],data[xm(idx)]),std::max(data[yp(idx)],data[ym(idx)])));
}

template <typename T>
T* Array2D<T>::dataAddress()
{ // returns address to first data element
  return(&(data.front()));
}

template <typename T>
std::vector<T> const& Array2D<T>::returnData() const
{ // returns a const ref to data std::vector
  return(data);
}

template <typename T>
std::shared_ptr<T> Array2D<T>::returnDataArray() const
{ // returns a pointer to a new (independent) array with data copied in.
  std::shared_ptr<T> pTData = std::make_shared<T>(getN());
  std::memcpy(pTData, &data, getN() * sizeof(T));
  return pTData;
}


template <typename T>
Array2D<T> Array2D<T>::duplicateArray2D() const
{ // deep copy of Array2D data
  return *this;
}

template <typename T>
void Array2D<T>::transposeData()
{
  assert(m==n);
  for(idx_t ii=0;ii<n;++ii)
    for(idx_t jj=0;jj<n;++jj)
      std::swap(data[ii+n*jj], data[jj+n*ii]);
}

template <typename T>
inline double Array2D<T>::getX(idx_t const idx) const
{
  return(static_cast<double>(getXidx(idx)) * dx);
}

template <typename T>
inline double Array2D<T>::getY(idx_t const idx) const
{
  return(static_cast<double>(getYidx(idx)) * dy);
}

template <typename T>
inline idx_t Array2D<T>::getXidx(idx_t const idx) const
{
  return(idx/m);
}

template <typename T>
inline idx_t Array2D<T>::getYidx(idx_t const idx) const
{
  return(idx%m);
}

template <typename T>
idx_t Array2D<T>::countVal(T const value) const
{
  return std::count_if(data.begin(), data.end(), [&value](T const &t)->bool{return value == t;});
}

template <typename T>
inline T& Array2D<T>::operator[] ( idx_t I )
{
  return data[I];
}

template <typename T>
inline T const& Array2D<T>::operator[] ( idx_t I ) const
{
  return data[I];
}

template <typename T>
inline Array2D<T>& Array2D<T>::operator= ( Array2D<T> const& Vec )
{
  assert(m == Vec.getm());
  assert(n == Vec.getn());
  data = Vec.data;
  return *this;
}

template <typename T>
inline Array2D<T>& Array2D<T>::operator*= ( T const Value )
{
  std::transform(data.begin(), data.end(), data.begin(),
                 [&Value](T const &t)->T { return Value*t; });
  return *this;
}

template <typename T>
inline Array2D<T>& Array2D<T>::operator/= ( T const Value )
{
  assert(Value != T(0));
  std::transform(data.begin(), data.end(), data.begin(),
                 [&Value](T const &t)->T { return t/Value; });
  return *this;
}

template <typename T>
inline Array2D<T>& Array2D<T>::operator+= ( const Array2D<T> &Vec )
{
  assert(Vec.getN() == getN());
  std::transform(data.begin(), data.end(), Vec.returnData(), data.begin(),
                 [](T const &inp1, T const &inp2)->T { return inp1+inp2; });
  return *this;
}

template <typename T>
inline Array2D<T>& Array2D<T>::operator-= ( const Array2D<T> &Vec )
{
  assert(Vec.getN() == getN());
  std::transform(data.begin(), data.end(), Vec.returnData(), data.begin(),
                 [](T const &inp1, T const &inp2)->T { return inp1-inp2; });
  return *this;
}

template <typename T>
inline void Array2D<T>::setAll( T const val )
{
  fillWithValue(val); // setAll is a redundant API
}

template <typename T>
inline void Array2D<T>::addMultiple ( const Array2D<T>& Vec, T Factor )
{
  assert(getN() == Vec.getN());
  std::transform(data.begin(), data.end(), Vec.returnDataArray(), data.begin(),
                 [&Factor](T const &t1, T const &t2)->T { return t1 + t2*Factor; });
}

template <typename T>
inline void Array2D<T>::save ( const char *FileName ) const {
  FILE *out;
  out = fopen ( FileName, "wb" );
  fwrite ( &data.front(), sizeof ( T ), N, out );
  fclose ( out );
}

template <typename T>
inline void Array2D<T>::load ( const char *FileName ) {
  FILE *in;
  in = fopen ( FileName, "rb" );
  fread ( &data.front(), sizeof ( T ), N, in );
  fclose ( in );
}

template <typename T>
inline void Array2D<T>::saveASCII ( const char *FileName ) const {
  std::ofstream file ( FileName );
  for ( idx_t i = 0; i < m; ++i ) {
    for ( idx_t j = 0; j < n; ++j )
      file << get( i, j ) << " ";
    file << std::endl;
  }
  file.close();
}

template <typename T>
inline void Array2D<T>::loadASCII ( const char *FileName ) {
  std::ifstream file( FileName );
  for ( idx_t i = 0; i < m; ++i )
    for ( idx_t j = 0; j < n; ++j ) {
      T value;
      file >> value;
      put( value, i, j );
    }
  file.close();
}

#endif
