#ifndef _ARRAY3D_HPP_
#define _ARRAY3D_HPP_

#include<vector>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include "defs.h"
#include<cassert>

using idx_t = std::size_t;

template <typename T> 
class Array3D
{
protected:
  std::vector<T> data; // contains (m x n x p) = N elements of type T
  const idx_t m;          // number of elements in the y-direction
  const idx_t n;          // number of elements in the x-direction
  const idx_t k;          // number of elements in the z-direction
  const idx_t N;          // total number of elements
  const double dx;         // grid spacing in x
  const double dy;         // grid spacing in y
  const double dz;         // grid spacing in z
  explicit Array3D<T>(); // don't use empty constructor

  idx_t sub2ind(const idx_t ii, const idx_t jj, const idx_t kk) const;

public: 
  explicit Array3D<T>(const idx_t nn);
  explicit Array3D<T>(const idx_t mm, const idx_t nn, const idx_t kk);
  explicit Array3D<T>(const T* indata, const idx_t mm, const idx_t nn, const idx_t kk);
  explicit Array3D<T>(const idx_t mm, const idx_t nn, const idx_t kk, const double _dx, const double _dy, const double _dz);
  explicit Array3D<T>(const T* indata, const idx_t mm, const idx_t nn, const idx_t kk, const double _dx, const double _dy, const double _dz);
  explicit Array3D<T>(const Array3D<T> &Input, const int CopyFlag);
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
  T sum() const;
  void fillWithValue(const T value);
  T sixNborMin(const idx_t idx) const;
  T sixNborMax(const idx_t idx) const;
  std::vector<T> returnData() const;
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
  void addMultiple ( const Array3D<T>& Vec, T Factor );
};

template <typename T> 
Array3D<T>::Array3D(const idx_t nn) :
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
Array3D<T>::Array3D(const idx_t mm, const idx_t nn, const idx_t kk) :
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
Array3D<T>::Array3D(const T* indata, const idx_t mm, const idx_t nn, const idx_t kk) :
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
  for(idx_t ii=0; ii<N; ++ii) // deep copy
    data[ii] = indata[ii];
}

template <typename T> 
Array3D<T>::Array3D(const idx_t mm, const idx_t nn, const idx_t kk, const double _dx, const double _dy, const double _dz) :
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
Array3D<T>::Array3D(const T* indata, const idx_t mm, const idx_t nn, const idx_t kk, const double _dx, const double _dy, const double _dz) :
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
  for(idx_t ii=0; ii<N; ++ii) // deep copy
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
idx_t Array3D<T>::sub2ind(const idx_t ii, const idx_t jj, const idx_t kk) const
{
  assert((ii+m*jj+m*n*kk) < N);
  return(ii+m*(jj+n*kk));
}

template <typename T>
T Array3D<T>::get(const idx_t idx) const
{
  return(data[idx]);
}

template <typename T>
T Array3D<T>::get(const idx_t ii, const idx_t jj, const idx_t kk) const
{
  return(data[sub2ind(ii,jj,kk)]);
}

template <typename T>
idx_t Array3D<T>::getm() const
{
  return(m);
}

template <typename T>
idx_t Array3D<T>::getn() const
{
  return(n);
}

template <typename T>
idx_t Array3D<T>::getk() const
{
  return(k);
}

template <typename T>
idx_t Array3D<T>::getN() const
{
  return(N);
}

template <typename T>
double Array3D<T>::getdx() const
{
  return(dx);
}

template <typename T>
double Array3D<T>::getdy() const
{
  return(dy);
}

template <typename T>
double Array3D<T>::getdz() const
{
  return(dz);
}

template <typename T>
double Array3D<T>::lenx() const
{
  return(static_cast<double>(n) * dx);
}

template <typename T>
double Array3D<T>::leny() const
{
  return(static_cast<double>(m) * dy);
}

template <typename T>
double Array3D<T>::lenz() const
{
  return(static_cast<double>(k) * dz);
}

template <typename T>
idx_t Array3D<T>::size() const
{
  return(N);
}

template <typename T>
void Array3D<T>::put(const T value, const idx_t idx)
{
  data[idx] = value;
}

template <typename T>
void Array3D<T>::put(const T value, const idx_t ii, const idx_t jj, const idx_t kk)
{
  data[sub2ind(ii,jj,kk)] = value;
}

template <typename T>
void Array3D<T>::putadd(const T value, const idx_t idx)
{
  data[idx] += value;
}

template <typename T>
void Array3D<T>::putadd(const T value, const idx_t ii, const idx_t jj, const idx_t kk)
{
  data[sub2ind(ii,jj,kk)] += value;
}

template <typename T>
T Array3D<T>::maxval() const
{
  assert(N>0);
  T mval = data[0];
  for(idx_t ii=1;ii<N;++ii)
    mval = mymax(mval,data[ii]);
  return(mval);
}

template <typename T>
T Array3D<T>::minval() const
{
  assert(N>0);
  T mval = data[0];
  for(idx_t ii=1;ii<N;++ii)
    mval = mymin(mval,data[ii]);
  return(mval);
}

template <typename T>
T Array3D<T>::minabsval() const
{
  assert(N>0);
  T mval = data[0]*mysign(data[0]);
  for(idx_t ii=1;ii<N;++ii)
    mval = mymin(mval,data[ii]*mysign(data[ii]));
  return(mval);
}

template <typename T>
T Array3D<T>::maxabsval() const
{
  assert(N>0);
  T mval = data[0]*mysign(data[0]);
  for(idx_t ii=1;ii<N;++ii)
    mval = mymax(mval,data[ii]*mysign(data[ii]));
  return(mval);
}

template <typename T>
T Array3D<T>::sum() const
{
  T val = 0;
  for( int ii=0; ii<N; ++ii )
    val += data[ii];
  return val;
}

template <typename T>
void Array3D<T>::fillWithValue(const T value)
{
  for(idx_t ii=0; ii<N; ++ii)
    data[ii] = value;
}

template <typename T>
T Array3D<T>::sixNborMin(const idx_t idx) const
{
  return(mymin(mymin(mymin(data[xp(idx)],data[xm(idx)]),mymin(data[yp(idx)],data[ym(idx)])),mymin(data[zp(idx)],data[zm(idx)])));
}

template <typename T>
T Array3D<T>::sixNborMax(const idx_t idx) const
{
  return(mymax(mymax(mymax(data[xp(idx)],data[xm(idx)]),mymax(data[yp(idx)],data[ym(idx)])),mymax(data[zp(idx)],data[zm(idx)])));
}

template <typename T>
std::vector<T> Array3D<T>::returnData() const
{ // returns a copy of the data vector
  return(data);
}

template <typename T>
double Array3D<T>::getX(const idx_t idx) const
{
  return(static_cast<double>(getXidx(idx)) / static_cast<double>(n));
}

template <typename T>
double Array3D<T>::getY(const idx_t idx) const
{
  return(static_cast<double>(getYidx(idx)) / static_cast<double>(m));
}

template <typename T>
double Array3D<T>::getZ(const idx_t idx) const
{
  return(static_cast<double>(getZidx(idx)) / static_cast<double>(k));
}

template <typename T>
idx_t Array3D<T>::getXidx(const idx_t idx) const
{
  return((idx%(m*n))/m);
}

template <typename T>
idx_t Array3D<T>::getYidx(const idx_t idx) const
{
  return((idx%(m*n))%m);
}

template <typename T>
idx_t Array3D<T>::getZidx(const idx_t idx) const
{
  return(idx/(m*n));
}

template <typename T>
idx_t Array3D<T>::xp(const idx_t idx) const
{
  return( (idx%(m*n))<(n*m-m) ? (idx +m) : (idx+m-m*n) );
}

template <typename T>
idx_t Array3D<T>::xm(const idx_t idx) const
{
  return( (idx%(m*n))>(m-1) ? (idx -m) : (idx-m+m*n) );
}

template <typename T>
idx_t Array3D<T>::yp(const idx_t idx) const
{
  return( ((idx%(m*n))%m)==(m-1) ? (idx-m+1) : (idx+1) );
}

template <typename T>
idx_t Array3D<T>::ym(const idx_t idx) const
{
  return( ((idx%(m*n))%m)==(0) ? (idx+m-1) : (idx-1) );
}

template <typename T>
idx_t Array3D<T>::zp(const idx_t idx) const
{
  return( (idx < (m*n*(k-1))) ? (idx+m*n) : (idx+m*n-m*n*k) );
}

template <typename T>
idx_t Array3D<T>::zm(const idx_t idx) const
{
  return( (idx > (m*n-1)) ? (idx-m*n) : (idx-m*n+m*n*k) );
}

template <typename T>
T Array3D<T>::getxp(const idx_t idx) const
{
  return(data[xp(idx)]);
}

template <typename T>
T Array3D<T>::getxm(const idx_t idx) const
{
  return(data[xm(idx)]);
}

template <typename T>
T Array3D<T>::getyp(const idx_t idx) const
{
  return(data[yp(idx)]);
}

template <typename T>
T Array3D<T>::getym(const idx_t idx) const
{
  return(data[ym(idx)]);
}

template <typename T>
T Array3D<T>::getzp(const idx_t idx) const
{
  return(data[zp(idx)]);
}

template <typename T>
T Array3D<T>::getzm(const idx_t idx) const
{
  return(data[zm(idx)]);
}

template <typename T>
idx_t Array3D<T>::countVal(const T value) const
{
  idx_t count = 0;
  for(idx_t ii=0; ii<N; ++ii)
    if(data[ii] == value)
      count++;
  return(count);
}

template <typename T>
T& Array3D<T>::operator[] ( idx_t I )
{
  return data[I];
}

template <typename T>
const T& Array3D<T>::operator[] ( idx_t I ) const
{
  return(data[I]);
}

template <typename T>
Array3D<T>& Array3D<T>::operator= ( const Array3D<T>& Vec )
{
  assert(m == Vec.getm());
  assert(n == Vec.getn());
  assert(k == Vec.getk());
  data = Vec.data;
  return(*this);
}

template <typename T>
Array3D<T>& Array3D<T>::operator*= ( const T Value )
{
  for ( idx_t i = 0; i < N; ++i )
    data[i] *= Value;
  return(*this);
}

template <typename T>
Array3D<T>& Array3D<T>::operator/= ( const T Value )
{
  for ( idx_t i = 0; i < N; ++i )
    data[i] /= Value;
  return(*this);
}

template <typename T>
Array3D<T>& Array3D<T>::operator+= ( const Array3D<T> &Vec )
{
  for ( idx_t i = 0; i < N; ++i )
    data[i] += Vec.get( i );
  return(*this);
}

template <typename T>
Array3D<T>& Array3D<T>::operator-= ( const Array3D<T> &Vec )
{
  for ( idx_t i = 0; i < N; ++i )
    data[i] -= Vec.get( i );
  return(*this);
}

template <typename T>
void Array3D<T>::addMultiple ( const Array3D<T>& Vec, T Factor )
{
  for ( idx_t i = 0; i < N; ++i )
    data[i] += Vec.get( i ) * Factor;
}

#endif
