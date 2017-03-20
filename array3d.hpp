#ifndef _ARRAY3D_HPP_
#define _ARRAY3D_HPP_

#include<vector>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<stack>
#include "defs.h"
#include<cassert>

using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::stack;

template <typename T> 
class Array3D
{

protected:
  vector<T> data; // contains (m x n x p) = N elements of type T
  const size_t m;          // number of elements in the y-direction
  const size_t n;          // number of elements in the x-direction
  const size_t k;          // number of elements in the z-direction
  const size_t N;          // total number of elements
  const double dx;         // grid spacing in x
  const double dy;         // grid spacing in y
  const double dz;         // grid spacing in z
  explicit Array3D<T>(); // don't use empty constructor

  inline size_t sub2ind(const size_t ii, const size_t jj, const size_t kk) const;

public: 
  explicit Array3D<T>(const size_t nn);
  explicit Array3D<T>(const size_t mm, const size_t nn, const size_t kk);
  explicit Array3D<T>(const T* indata, const size_t mm, const size_t nn, const size_t kk);
  explicit Array3D<T>(const size_t mm, const size_t nn, const size_t kk, const double _dx, const double _dy, const double _dz);
  explicit Array3D<T>(const T* indata, const size_t mm, const size_t nn, const size_t kk, const double _dx, const double _dy, const double _dz);
  explicit Array3D<T>(const Array3D<T> &Input, const int CopyFlag);
  Array3D<T>(const Array3D<T> &Input); 
  inline T get(const size_t idx) const;
  inline T get(const size_t ii, const size_t jj, const size_t kk) const;
  inline size_t getm() const;
  inline size_t getn() const;
  inline size_t getk() const;
  inline size_t getN() const;
  inline double getdx() const;
  inline double getdy() const;
  inline double getdz() const;
  inline double lenx() const;
  inline double leny() const;
  inline double lenz() const;
  inline size_t size() const;
  inline void put(const T value, const size_t idx);
  inline void put(const T value, const size_t ii, const size_t jj, const size_t kk);
  inline void putadd(const T value, const size_t idx);
  inline void putadd(const T value, const size_t ii, const size_t jj, const size_t kk);
  void copyData(const Array3D<T> &input);

  inline T maxval() const;
  inline T minval() const;
  inline T minabsval() const;
  inline T maxabsval() const;
  inline T sum() const;
  inline void fillWithValue(const T value);
  inline T sixNborMin(const size_t idx) const;
  inline T sixNborMax(const size_t idx) const;
  T* dataAddress();
  vector<T> returnData() const;
  T* returnDataArray() const;
  Array3D<T> duplicateArray3D() const;
  inline double getX(const size_t idx) const;
  inline double getY(const size_t idx) const;
  inline double getZ(const size_t idx) const;
  inline size_t getXidx(const size_t idx) const;
  inline size_t getYidx(const size_t idx) const;
  inline size_t getZidx(const size_t idx) const;
  inline size_t xp(const size_t idx) const;
  inline size_t xm(const size_t idx) const;
  inline size_t yp(const size_t idx) const;
  inline size_t ym(const size_t idx) const;
  inline size_t zp(const size_t idx) const;
  inline size_t zm(const size_t idx) const;
  inline T getxp(const size_t idx) const;
  inline T getxm(const size_t idx) const;
  inline T getyp(const size_t idx) const;
  inline T getym(const size_t idx) const;
  inline T getzp(const size_t idx) const;
  inline T getzm(const size_t idx) const;

  size_t countVal(const T value) const;

  // algebraic operations
  inline T& operator[] ( size_t I );
  inline const T& operator[] ( size_t I ) const;
  inline Array3D<T>& operator= ( const Array3D<T>& Vec );
  inline Array3D<T>& operator*= ( const T Value );
  inline Array3D<T>& operator/= ( const T Value );
  inline Array3D<T>& operator+= ( const Array3D<T> &Vec );
  inline Array3D<T>& operator-= ( const Array3D<T> &Vec );
  inline void setAll( const T val );
  inline void addMultiple ( const Array3D<T>& Vec, T Factor );
};

template <typename T> 
Array3D<T>::Array3D(const size_t nn) :
  m(nn),
  n(nn),
  k(nn),
  N(m*n*k),
  dx(1.0f/static_cast<double>(n)),
  dy(1.0f/static_cast<double>(m)),
  dz(1.0f/static_cast<double>(k))
{
  assert(m>0);
  data.resize(N);
}

template <typename T> 
Array3D<T>::Array3D(const size_t mm, const size_t nn, const size_t kk) :
  m(mm),
  n(nn),
  k(kk),
  N(m*n*k),
  dx(1.0f/static_cast<double>(n)),
  dy(1.0f/static_cast<double>(m)),
  dz(1.0f/static_cast<double>(k))
{
  assert(m>0);
  assert(n>0);
  assert(k>0);
  data.resize(N);
}

template <typename T>
Array3D<T>::Array3D(const T* indata, const size_t mm, const size_t nn, const size_t kk) :
  m(mm),
  n(nn),
  k(kk),
  N(mm*nn*kk),
  dx(1.0f/static_cast<double>(n)),
  dy(1.0f/static_cast<double>(m)),
  dz(1.0f/static_cast<double>(k))
{ // assume indata contains at least (exactly) N entries
  assert(m>0);
  assert(n>0);
  assert(k>0);
  data.resize(N);
  for(size_t ii=0; ii<N; ++ii) // deep copy
    data[ii] = indata[ii];
}

template <typename T> 
Array3D<T>::Array3D(const size_t mm, const size_t nn, const size_t kk, const double _dx, const double _dy, const double _dz) :
  m(mm),
  n(nn),
  k(kk),
  N(mm*nn*kk),
  dx(_dx),
  dy(_dy),
  dz(_dz)
{
  assert(m>0);
  assert(n>0);
  assert(k>0);
  assert(dx>0.f);
  assert(dy>0.f);
  assert(dz>0.f);
  data.resize(N);
}

template <typename T>
Array3D<T>::Array3D(const T* indata, const size_t mm, const size_t nn, const size_t kk, const double _dx, const double _dy, const double _dz) :
  m(mm),
  n(nn),
  k(kk),
  N(mm*nn*kk),
  dx(_dx),
  dy(_dy),
  dz(_dz)
{ // assume indata contains at least (exactly) N entries
  assert(m>0);
  assert(n>0);
  assert(k>0);
  assert(dx>0.f);
  assert(dy>0.f);
  assert(dz>0.f);
  data.resize(N);
  for(size_t ii=0; ii<N; ++ii) // deep copy
    data[ii] = indata[ii];
}

template <typename T>
Array3D<T>::Array3D( const Array3D<T> &Input, const int CopyFlag) :
  m(Input.m),
  n(Input.n),
  k(Input.k),
  N(Input.N),
  dx(Input.dx),
  dy(Input.dy),
  dz(Input.dz)
{
  switch ( CopyFlag ) {
  case 0:  // deep copy
    data = Input.data;
    break;
  default: // structure copy
    data.resize(N);
    break;
  }
}

template <typename T>
Array3D<T>::Array3D( const Array3D<T> &Input) :
  m(Input.m),
  n(Input.n),
  k(Input.k),
  N(Input.N),
  dx(Input.dx),
  dy(Input.dy),
  dz(Input.dz)
{
  data = Input.data;
}

template <typename T>
inline size_t Array3D<T>::sub2ind(const size_t ii, const size_t jj, const size_t kk) const
{
  assert((ii+m*jj+m*n*kk) < N);
  return(ii+m*(jj+n*kk));
}

template <typename T>
inline T Array3D<T>::get(const size_t idx) const
{
  return(data[idx]);
}

template <typename T>
inline T Array3D<T>::get(const size_t ii, const size_t jj, const size_t kk) const
{
  return(data[sub2ind(ii,jj,kk)]);
}

template <typename T>
inline size_t Array3D<T>::getm() const
{
  return(m);
}

template <typename T>
inline size_t Array3D<T>::getn() const
{
  return(n);
}

template <typename T>
inline size_t Array3D<T>::getk() const
{
  return(k);
}

template <typename T>
inline size_t Array3D<T>::getN() const
{
  return(N);
}

template <typename T>
inline double Array3D<T>::getdx() const
{
  return(dx);
}

template <typename T>
inline double Array3D<T>::getdy() const
{
  return(dy);
}

template <typename T>
inline double Array3D<T>::getdz() const
{
  return(dz);
}

template <typename T>
inline double Array3D<T>::lenx() const
{
  return(static_cast<double>(n) * dx);
}

template <typename T>
inline double Array3D<T>::leny() const
{
  return(static_cast<double>(m) * dy);
}

template <typename T>
inline double Array3D<T>::lenz() const
{
  return(static_cast<double>(k) * dz);
}

template <typename T>
inline size_t Array3D<T>::size() const
{
  return(N);
}

template <typename T>
inline void Array3D<T>::put(const T value, const size_t idx)
{
  data[idx] = value;
}

template <typename T>
inline void Array3D<T>::put(const T value, const size_t ii, const size_t jj, const size_t kk)
{
  data[sub2ind(ii,jj,kk)] = value;
}

template <typename T>
inline void Array3D<T>::putadd(const T value, const size_t idx)
{
  data[idx] += value;
}

template <typename T>
inline void Array3D<T>::putadd(const T value, const size_t ii, const size_t jj, const size_t kk)
{
  data[sub2ind(ii,jj,kk)] += value;
}

template <typename T>
void Array3D<T>::copyData(const Array3D<T> &input)
{ // copy data after checking for agreement of sizes
  if((input.getm() != m) || (input.getn() != n) || (input.getk() != k))
    cout << "Data does not agree in size. Skipping call to copyData" << endl;
  else
    data = input.returnData();
}

template <typename T>
inline T Array3D<T>::maxval() const
{
  assert(N>0);
  T mval = data[0];
  for(size_t ii=1;ii<N;++ii)
    mval = mymax(mval,data[ii]);
  return(mval);
}

template <typename T>
inline T Array3D<T>::minval() const
{
  assert(N>0);
  T mval = data[0];
  for(size_t ii=1;ii<N;++ii)
    mval = mymin(mval,data[ii]);
  return(mval);
}

template <typename T>
inline T Array3D<T>::minabsval() const
{
  assert(N>0);
  T mval = data[0]*mysign(data[0]);
  for(size_t ii=1;ii<N;++ii)
    mval = mymin(mval,data[ii]*mysign(data[ii]));
  return(mval);
}

template <typename T>
inline T Array3D<T>::maxabsval() const
{
  assert(N>0);
  T mval = data[0]*mysign(data[0]);
  for(size_t ii=1;ii<N;++ii)
    mval = mymax(mval,data[ii]*mysign(data[ii]));
  return(mval);
}

template <typename T>
inline T Array3D<T>::sum() const
{
  T val = 0;
  for( int ii=0; ii<N; ++ii )
    val += data[ii];
  return val;
}

template <typename T>
inline void Array3D<T>::fillWithValue(const T value)
{
  for(size_t ii=0; ii<N; ++ii)
    data[ii] = value;
}

template <typename T>
inline T Array3D<T>::sixNborMin(const size_t idx) const
{
  return(mymin(mymin(mymin(data[xp(idx)],data[xm(idx)]),mymin(data[yp(idx)],data[ym(idx)])),mymin(data[zp(idx)],data[zm(idx)])));
}

template <typename T>
inline T Array3D<T>::sixNborMax(const size_t idx) const
{
  return(mymax(mymax(mymax(data[xp(idx)],data[xm(idx)]),mymax(data[yp(idx)],data[ym(idx)])),mymax(data[zp(idx)],data[zm(idx)])));
}

template <typename T>
T* Array3D<T>::dataAddress()
{ // returns address to first data element
  return(&(data.front()));
}

template <typename T>
vector<T> Array3D<T>::returnData() const
{ // returns a copy of the data vector
  return(data);
}

template <typename T>
T* Array3D<T>::returnDataArray() const
{ // returns a pointer to a new (independent) array with data copied in.
  T* newArray = new T[getN()];
  for(size_t ii=0;ii<N;++ii)
    newArray[ii] = get(ii);
  return(newArray);
}

template <typename T>
Array3D<T> Array3D<T>::duplicateArray3D() const
{ // deep copy of Array3D data
  Array3D<T> v(m,n,k);
  v.data = this->data;
  return(v);
}

template <typename T>
inline double Array3D<T>::getX(const size_t idx) const
{
  return(static_cast<double>(getXidx(idx)) / static_cast<double>(n));
}

template <typename T>
inline double Array3D<T>::getY(const size_t idx) const
{
  return(static_cast<double>(getYidx(idx)) / static_cast<double>(m));
}

template <typename T>
inline double Array3D<T>::getZ(const size_t idx) const
{
  return(static_cast<double>(getZidx(idx)) / static_cast<double>(k));
}

template <typename T>
inline size_t Array3D<T>::getXidx(const size_t idx) const
{
  return((idx%(m*n))/m);
}

template <typename T>
inline size_t Array3D<T>::getYidx(const size_t idx) const
{
  return((idx%(m*n))%m);
}

template <typename T>
inline size_t Array3D<T>::getZidx(const size_t idx) const
{
  return(idx/(m*n));
}

template <typename T>
inline size_t Array3D<T>::xp(const size_t idx) const
{
  return( (idx%(m*n))<(n*m-m) ? (idx +m) : (idx+m-m*n) );
}

template <typename T>
inline size_t Array3D<T>::xm(const size_t idx) const
{
  return( (idx%(m*n))>(m-1) ? (idx -m) : (idx-m+m*n) );
}

template <typename T>
inline size_t Array3D<T>::yp(const size_t idx) const
{
  return( ((idx%(m*n))%m)==(m-1) ? (idx-m+1) : (idx+1) );
}

template <typename T>
inline size_t Array3D<T>::ym(const size_t idx) const
{
  return( ((idx%(m*n))%m)==(0) ? (idx+m-1) : (idx-1) );
}

template <typename T>
inline size_t Array3D<T>::zp(const size_t idx) const
{
  return( (idx < (m*n*(k-1))) ? (idx+m*n) : (idx+m*n-m*n*k) );
}

template <typename T>
inline size_t Array3D<T>::zm(const size_t idx) const
{
  return( (idx > (m*n-1)) ? (idx-m*n) : (idx-m*n+m*n*k) );
}

template <typename T>
inline T Array3D<T>::getxp(const size_t idx) const
{
  return(data[xp(idx)]);
}

template <typename T>
inline T Array3D<T>::getxm(const size_t idx) const
{
  return(data[xm(idx)]);
}

template <typename T>
inline T Array3D<T>::getyp(const size_t idx) const
{
  return(data[yp(idx)]);
}

template <typename T>
inline T Array3D<T>::getym(const size_t idx) const
{
  return(data[ym(idx)]);
}

template <typename T>
inline T Array3D<T>::getzp(const size_t idx) const
{
  return(data[zp(idx)]);
}

template <typename T>
inline T Array3D<T>::getzm(const size_t idx) const
{
  return(data[zm(idx)]);
}

template <typename T>
size_t Array3D<T>::countVal(const T value) const
{
  size_t count = 0;
  for(size_t ii=0; ii<N; ++ii)
    if(data[ii] == value)
      count++;
  return(count);
}

template <typename T>
inline T& Array3D<T>::operator[] ( size_t I )
{
  return data[I];
}

template <typename T>
inline const T& Array3D<T>::operator[] ( size_t I ) const
{
  return(data[I]);
}

template <typename T>
inline Array3D<T>& Array3D<T>::operator= ( const Array3D<T>& Vec )
{
  assert(m == Vec.getm());
  assert(n == Vec.getn());
  assert(k == Vec.getk());
  data = Vec.data;
  return(*this);
}

template <typename T>
inline Array3D<T>& Array3D<T>::operator*= ( const T Value )
{
  for ( size_t i = 0; i < N; ++i )
    data[i] *= Value;
  return(*this);
}

template <typename T>
inline Array3D<T>& Array3D<T>::operator/= ( const T Value )
{
  for ( size_t i = 0; i < N; ++i )
    data[i] /= Value;
  return(*this);
}

template <typename T>
inline Array3D<T>& Array3D<T>::operator+= ( const Array3D<T> &Vec )
{
  for ( size_t i = 0; i < N; ++i )
    data[i] += Vec.get( i );
  return(*this);
}

template <typename T>
inline Array3D<T>& Array3D<T>::operator-= ( const Array3D<T> &Vec )
{
  for ( size_t i = 0; i < N; ++i )
    data[i] -= Vec.get( i );
  return(*this);
}

template <typename T>
inline void Array3D<T>::setAll( const T val )
{
  for ( size_t i = 0; i < N; ++i )
    data[i] = val;
}

template <typename T>
inline void Array3D<T>::addMultiple ( const Array3D<T>& Vec, T Factor )
{
  for ( size_t i = 0; i < N; ++i )
    data[i] += Vec.get( i ) * Factor;
}

#endif
