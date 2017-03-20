#ifndef _ARRAY2D_HPP_
#define _ARRAY2D_HPP_

// To-do:
// [1] Doxygen
// [2] Sub-class HexArray2D (Array2D with hexagonal connectivity)

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

  void labelHexCC(T const value, const bool sign, idx_t const idx, Array2D<int> &label) const;
  void labelHexCCX(idx_t const idx, T const val, std::vector<idx_t> &changed);
  void labelHexCCY(idx_t const idx, T const val, std::vector<idx_t> &changed);

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
  inline std::vector<idx_t> connectedComponent(idx_t const idx, double const size);
  std::vector<idx_t> connectedComponent(idx_t const idx, double const size, T const changeval);
  std::vector<idx_t> connectedComponentPM(idx_t const idx, double const size);
  std::vector<idx_t> boundaryOfCC(const std::vector<idx_t> cc);
  idx_t diff1Nbors(T const val, idx_t idx);
  idx_t diff2Nbors(T const val, idx_t idx);
  inline T maxval() const;
  inline T minval() const;
  inline T minabsval() const;
  inline void fillWithValue(T const value);
  inline T fourNborMin(idx_t const idx) const;
  inline T fourNborMax(idx_t const idx) const;
  void dump() const;
  void dump(std::string const &filename) const;
  void dumpbooleq(T const val) const;
  void dumpboolgeq(T const val) const;
  void dumpboolleq(T const val) const;
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

  inline idx_t hexN2(idx_t const idx, int const nbor) const;
  inline idx_t hexN(idx_t const idx, int const nbor) const;
  inline idx_t hexNNP2(idx_t const idx, int const nbor) const;
  inline idx_t hexNNP(idx_t const idx, int const nbor) const;
  Array2D<int> labelHex(T const value, const bool sign) const;
  bool labelHexFromIdxX(idx_t const idx);
  bool labelHexFromIdxY(idx_t const idx);

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
  inline T dotProd ( const Array2D<T>& Vec ) const;
  inline T getNormSqr() const;
  inline T getNorm() const;

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
inline std::vector<idx_t> Array2D<T>::connectedComponent(idx_t const idx, double const size) 
{ // calls connectedComponent with changeval==-1
  return(connectedComponent(idx,size,T(-1)));
}

template <typename T>
std::vector<idx_t> Array2D<T>::connectedComponent(idx_t const idx, double const size, T const changeval) 
{ // returns a std::vector containing the index of all pixels in the connected component identified by q==q[idx] containing idx
  // Change values to changeval when touched to indicate they've been
  // checked; then change all the values back to the original
  // value at the end.
  // Assumes that changeval is not a value used by the array for any
  // other purpose.
  // size is an estimate of the proportion of the grid covered
  // by the connected component
  T value = get(idx);
  idx_t curr = 0;
  std::vector<idx_t> v;
  v.reserve(static_cast<idx_t>(ceil(size)*N));
  v.push_back(idx);
  put(-1,idx); 
  
  while(curr < v.size())
  {
    if(get(xp(v[curr])) == value)
    { v.push_back(xp(v[curr])); put(changeval, xp(v[curr]));}
    if(get(xm(v[curr])) == value)
    { v.push_back(xm(v[curr])); put(changeval, xm(v[curr]));}
    if(get(yp(v[curr])) == value)
    { v.push_back(yp(v[curr])); put(changeval, yp(v[curr]));}
    if(get(ym(v[curr])) == value)
    { v.push_back(ym(v[curr])); put(changeval, ym(v[curr]));}
    curr++;
  }
  for(idx_t ii=0; ii<v.size(); ++ii)
    put(value,v[ii]);

  return(v);
}

template <typename T>
std::vector<idx_t> Array2D<T>::connectedComponentPM(idx_t const idx, double const size) 
{ // returns a std::vector containing the index of all pixels in the connected component identified by q== \pm q[idx] containing idx
  // Change values to 0 when touched to indicate they've been
  // checked; then change all the values back to the original
  // value at the end.
  // Assumes that 0 is not a value used by the array for any
  // other purpose.
  // size is an estimate of the proportion of the grid covered
  // by the connected component
  T value = abs(get(idx));
  idx_t curr = 0;
  std::vector<idx_t> v;  // indices
  std::vector<bool> s; // sign (+:true -:false)
  v.reserve(static_cast<idx_t>(ceil(size)*N));
  s.reserve(static_cast<idx_t>(ceil(size)*N));
  v.push_back(idx);
  s.push_back(get(idx) > 0.0f);
  put(0,idx); 
  
  while(curr < v.size())
  {
    if(abs(getxp(v[curr])) == value)
    { 
      v.push_back(xp(v[curr])); 
      s.push_back(getxp(v[curr])>0.0f);
      put(0,xp(v[curr]));
    }
    if(abs(getxm(v[curr])) == value)
    { 
      v.push_back(xm(v[curr])); 
      s.push_back(getxm(v[curr])>0.0f);
      put(0,xm(v[curr]));
    }
    if(abs(getyp(v[curr])) == value)
    { 
      v.push_back(yp(v[curr])); 
      s.push_back(getyp(v[curr])>0.0f);
      put(0,yp(v[curr]));
    }
    if(abs(getym(v[curr])) == value)
    { 
      v.push_back(ym(v[curr])); 
      s.push_back(getym(v[curr])>0.0f);
      put(0,ym(v[curr]));
    }
    curr++;
  }
  for(idx_t ii=0; ii<v.size(); ++ii)
    if(s[ii])
      put(value,v[ii]);
    else
      put(-value,v[ii]);

  return(v);
}

template <typename T>
std::vector<idx_t> Array2D<T>::boundaryOfCC(const std::vector<idx_t> cc)
{
  T value = get(cc[0]);
  std::vector<idx_t> v;
  v.reserve(cc.size());
  for(idx_t ii=0;ii<cc.size();++ii)
    if(get(xp(cc[ii]))!=value)
      v.push_back(cc[ii]);
    else if(get(xm(cc[ii]))!=value)
      v.push_back(cc[ii]);
    else if(get(yp(cc[ii]))!=value)
      v.push_back(cc[ii]);
    else if(get(ym(cc[ii]))!=value)
      v.push_back(cc[ii]);
  return(v);
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
void Array2D<T>::dump() const
{
  std::cout.setf(std::ios::fixed,std::ios::floatfield);
  std::cout.setf(std::ios::showpos);
  std::cout.precision(5);
  for(idx_t ii=0; ii < m; ++ii)
  {
    for(idx_t jj=0; jj < n; ++jj)
      std::cout << get(ii,jj) << " ";
    std::cout << std::endl;
  }
  std::cout.unsetf(std::ios::fixed);
  std::cout.unsetf(std::ios::floatfield);
  std::cout.unsetf(std::ios::showpos);
}

template <typename T>
void Array2D<T>::dump(std::string const &filename) const
{
  std::ofstream of;
  of.open(filename.c_str());
  
  of.setf(std::ios::fixed,std::ios::floatfield);
  of.setf(std::ios::showpos);
  of.precision(5);
  for(idx_t ii=0; ii < m; ++ii)
  {
    for(idx_t jj=0; jj < n; ++jj)
      of << get(ii,jj) << " ";
    of << std::endl;
  }
  of.unsetf(std::ios::fixed);
  of.unsetf(std::ios::floatfield);
  of.unsetf(std::ios::showpos);
  of.close();
  of.clear();
}

template <typename T>
void Array2D<T>::dumpbooleq(T const val) const
{ // outputs 1 
  for(idx_t ii=0; ii < m; ++ii)
  {
    for(idx_t jj=0; jj < n; ++jj)
      std::cout << (get(ii,jj)==val);
    std::cout << std::endl;
  }
}

template <typename T>
void Array2D<T>::dumpboolgeq(T const val) const
{
  for(idx_t ii=0; ii < m; ++ii)
  {
    for(idx_t jj=0; jj < n; ++jj)
      std::cout << (get(ii,jj)>=val);
    std::cout << std::endl;
  }
}

template <typename T>
void Array2D<T>::dumpboolleq(T const val) const
{
  for(idx_t ii=0; ii < m; ++ii)
  {
    for(idx_t jj=0; jj < n; ++jj)
      std::cout << (get(ii,jj)<=val);
    std::cout << std::endl;
  }
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
idx_t Array2D<T>::diff1Nbors(T const val, idx_t idx)
{ // count the number of 4-neighbors of idx different than val
  return (val != get(xp(idx))) + (val != get(xm(idx))) 
    + (val != get(yp(idx))) + (val != get(ym(idx)));
}

template <typename T>
idx_t Array2D<T>::diff2Nbors(T const val, idx_t idx)
{ // count the number of 2nd-nearest-neighbors of idx different than val
  return  (val != get(xp(yp(idx)))) + (val != get(xm(yp(idx)))) 
    + (val != get(xp(ym(idx)))) + (val != get(xm(ym(idx))));
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
inline idx_t Array2D<T>::hexN2(idx_t const idx, int const nbor) const
{
  assert((nbor>=1) && (nbor<=6));
  switch(nbor)
  {
  case 1:
    return(xp(idx));
    break;
  case 2:
    return(xp(yp(idx)));
    break;
  case 3:
    return(yp(idx));
    break;
  case 4:
    return(xm(idx));
    break;
  case 5:
    return(xm(ym(idx)));
    break;
  case 6:
    return(ym(idx));
    break;
  default:
    std::cout << "hexN2 only accepts values between 1 and 6 in the second argument." << std::endl;
    return(-1);
  }
}

template <typename T>
inline idx_t Array2D<T>::hexN(idx_t const idx, int const nbor) const
{
  assert((nbor>=1) && (nbor<=6));
  switch(nbor)
  {
  case 1:
    return(xp(idx));
    break;
  case 2:
    if(((idx%m)%2)==0) // even row
      return(xp(yp(idx)));
    else
      return(yp(idx));
    break;
  case 3:
    if(((idx%m)%2)==0) // even row
      return(yp(idx));
    else
      return(xm(yp(idx)));
    break;
  case 4:
    return(xm(idx));
    break;
  case 5:
    if(((idx%m)%2)==0) // even row
      return(ym(idx));
    else
      return(xm(ym(idx)));
    break;
  case 6:
    if(((idx%m)%2)==0) // even row
      return(xp(ym(idx)));
    else
      return(ym(idx));
    break;
  default:
    std::cout << "hexN only accepts values between 1 and 6 in the second argument." << std::endl;
    return(-1);
  }
}

template <typename T>
inline idx_t Array2D<T>::hexNNP2(idx_t const idx, int const nbor) const
{
  assert((nbor>=1) && (nbor<=6));
  // check boundary cases, if you would cross a periodic boundary, return idx instead
  if((idx/m)==0) // left boundary
  {
    if((nbor ==4)||(nbor == 5))
      return(idx);
  }
  else if((idx/m)==(m-1)) // right boundary
  {
    if((nbor == 1)||(nbor == 2))
      return(idx);
  }
  if((idx%m)==0) // bottom row
  {
    if((nbor == 5) || (nbor == 6))
      return(idx);
  }
  else if((idx%m)==(m-1)) // top row
  {
    if((nbor == 2) || (nbor == 3))
      return(idx);
  }
  // if boundary cases fall through, return usual neighbor
  return(hexN2(idx,nbor));
}

template <typename T>
inline idx_t Array2D<T>::hexNNP(idx_t const idx, int const nbor) const
{
  assert((nbor>=1) && (nbor<=6));
  // check boundary cases, if you would cross a periodic boundary, return idx instead
  if((idx/m)==0) // left boundary
  {
    if(nbor == 4)
      return(idx);
    else if(((idx%m)%2)==1) // odd row
      if((nbor == 3) || (nbor == 5))
        return(idx);
  }
  else if((idx/m)==(m-1)) // right boundary
  {
    if(nbor == 1)
      return(idx);
    else if(((idx%m)%2)==0) // even row
      if((nbor == 2) || (nbor == 6))
        return(idx);
  }
  if((idx%m)==0) // bottom row
  {
    if((nbor == 5) || (nbor == 6))
      return(idx);
  }
  else if((idx%m)==(m-1)) // top row
  { // m, n must be even for grid to match up - top row is odd
    if((nbor == 2) || (nbor == 3))
      return(idx);
  }
  // if boundary cases fall through, return usual neighbor
  return(hexN(idx,nbor));
}

template <typename T>
Array2D<int> Array2D<T>::labelHex(T const value, const bool sign) const
{ // return labeling of the connected components of data>=value (if sign==true)
  // and data<value (if sign==false), with hexagonal connectivity
  Array2D<int> label(m,n);
  label.fillWithValue(-1);

  int labelnum=0;
  for(idx_t ii=0;ii<N;++ii)
    if(label.get(ii) == -1)
      if(((sign) && (data[ii] >= value)) || ((!sign) && (data[ii] < value)))
      {
        label.put(labelnum,ii);
        labelnum++;
        labelHexCC(value,sign,ii,label);
      }
  return(label);
}

template <typename T>
void Array2D<T>::labelHexCC(T const value, const bool sign, idx_t const idx, Array2D<int> &label) const
{ // recursively add pieces of this connected component

  std::stack<idx_t> s;
  s.push(idx);
  idx_t stackval;
  while(!s.empty())
  {
    stackval = s.top();
    s.pop();
    for(int ii=1; ii<=6; ++ii)
    {
      if(label.get(hexN(stackval,ii)) == -1)
        if(((sign) && (data[hexN(stackval,ii)] >= value)) || ((!sign) && (data[hexN(stackval,ii)] < value)))
        {
          assert(label.get(stackval) != -1);
          label.put(label.get(stackval),hexN(stackval,ii));
          s.push(hexN(stackval,ii));
        }
    }
  }
}

template <typename T>
bool Array2D<T>::labelHexFromIdxY(idx_t const idx)
{ // Try to label from idx without crossing Y=0, and see if the component is "infinite" from top to bottom
  std::vector<idx_t> changed;
  T const val = data[idx];

  for(idx_t ii=0; ii<n; ++ii)
  {
    if(data[m*ii] == val)
    {
      data[m*ii] = 0;
      changed.push_back(m*ii);
      labelHexCCY(m*ii,val,changed);
    }
  }

  // check to see if a connection was made. Then return true/false
  bool connection = false;
  for(idx_t ii=0; ii<n; ++ii)
  {
    if((get(0,ii)==0) && (get(m-1,ii)==0))
    {
      connection = true;
      break;
    }
    if((get(0,ii)==0) && (get(m-1,(ii+1)%n)==0))
    {
      connection = true;
      break;
    }
  }

  for(idx_t ii=0; ii<changed.size(); ++ii)
    data[changed[ii]] = val;

  return(connection);
}

template <typename T>
bool Array2D<T>::labelHexFromIdxX(idx_t const idx)
{ // Try to label from idx without crossing X=0, and see if the component is "infinite" from left to right
  std::vector<idx_t> changed;
  T const val = data[idx];

  for(int ii=0; ii<m; ++ii)
  {
    if(data[ii] == val)
    {
      data[ii] = 0;
      changed.push_back(ii);
      labelHexCCX(ii,val,changed);
    }
  }

  bool connection = false;
  // check to see if a connection was made. Then return true/false
  for(idx_t ii=0; ii<m; ++ii)
  {
    if((get(ii,0)==0)&&(get(ii,n-1)==0))
    {
      connection = true;
      break;
    }
    if((ii%2) == 1) // odd row, check above and below, too
    {
      if((get(ii,0)==0) && (get((ii+1)%m,n-1)==0))
      {
        connection = true;
        break;
      }
      if((get(ii,0)==0) && (get((ii-1)%m,n-1)==0))
      {
        connection = true;
        break;
      }
    }
  }

  // change values back
  for(idx_t ii=0; ii<changed.size(); ++ii)
    data[changed[ii]] = val;
  return(connection);
}

template <typename T>
void Array2D<T>::labelHexCCY(idx_t const idx, T const val, std::vector<idx_t> &changed)
{ // recursively add pieces of this connected component
  std::stack<idx_t> s;
  s.push(idx);
  idx_t stackval;
  while(!s.empty())
  {
    stackval = s.top();
    s.pop();
    if((stackval % m) != 0) // check all 6 neighbors
      for(int ii=1; ii<=6; ++ii)
      {
        idx_t cstackval = hexN(stackval,ii);
        if(data[cstackval] == val)
        {
          data[cstackval] = 0;
          changed.push_back(cstackval);
          s.push(cstackval);
        }
      }
    else // (stackval%m) == 0: top row
      for(int ii=1; ii<=4; ++ii)
      {
        idx_t cstackval = hexN(stackval,ii);
        if(data[cstackval] == val)
        {
          data[cstackval] = 0;
          changed.push_back(cstackval);
          s.push(cstackval);
        }
      }
  }
}

template <typename T>
void Array2D<T>::labelHexCCX(idx_t const idx, T const val, std::vector<idx_t> &changed)
{ // recursively add pieces of this connected component
  std::stack<idx_t> s;
  s.push(idx);
  idx_t stackval;
  while(!s.empty())
  {
    stackval = s.top();
    s.pop();
    if((stackval / m) != 0) // check all 6 neighbors
      for(int ii=1; ii<=6; ++ii)
      {
        idx_t cstackval = hexN(stackval,ii);
        if(data[cstackval] == val)
        {
          data[cstackval] = 0;
          changed.push_back(cstackval);
          s.push(cstackval);
        }
    }
    else // (stackval/m) == 0: left column
    {
      int jj[3] = {1,2,6};
      for(int ii=0; ii<3; ++ii)
      {
        idx_t cstackval = hexN(stackval,jj[ii]);
        if(data[cstackval] == val)
        {
          data[cstackval] = 0;
          changed.push_back(cstackval);
          s.push(cstackval);
        }
      }
      if((stackval%m) == 0) // even row, also check 3 and 5
      {
        int jj[2] = {3,5};
        for(int ii=0; ii<2; ++ii)
        {
          idx_t cstackval = hexN(stackval,jj[ii]);
          if(data[cstackval] == val)
          {
            data[cstackval] = 0;
            changed.push_back(cstackval);
            s.push(cstackval);
          }
        }
      }
    }
  }
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
inline T Array2D<T>::dotProd ( const Array2D<T>& Vec ) const
{
  assert(getN() == Vec.getN());
  return std::inner_product(data.begin(), data.end(), Vec.returnDataArray(), T(0));
}

template <typename T>
inline T Array2D<T>::getNormSqr() const
{
  return std::accumulate(data.begin(), data.end(), T(0), [](T const &partSum, T const &t)->T { return partSum + t*t; });
}

template <typename T>
inline T Array2D<T>::getNorm() const
{
  return sqrt(getNormSqr());
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
