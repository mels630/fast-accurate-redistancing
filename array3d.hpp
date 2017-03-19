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
  inline bool onBndry(const size_t idx) const;
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
  inline vector<size_t> connectedComponent(const size_t idx, const double size);
  vector<size_t> connectedComponent(const size_t idx, const double size, const int changeval);
  vector<size_t> connectedComponentPM(const size_t idx, const double size);
  vector<int> connectedComponent(const int idx, const double size, const T compareval, bool(*compareFn)(const T, const T));
  inline static bool compareG(const T currval, const T compareval);
  inline static bool compareGeq(const T currval, const T compareval);
  inline static bool compareEq(const T currval, const T compareval);
  inline static bool compareLeq(const T currval, const T compareval);
  inline static bool compareL(const T currval, const T compareval);

  vector<size_t> boundaryOfCC(const vector<size_t> cc);
  size_t diff1Nbors(T val, size_t idx);
  size_t diff2Nbors(T val, size_t idx);
  inline T maxval() const;
  inline T minval() const;
  inline T minabsval() const;
  inline T maxabsval() const;
  inline T sum() const;
  inline void fillWithValue(const T value);
  inline T sixNborMin(const size_t idx) const;
  inline T sixNborMax(const size_t idx) const;
  void dump() const;
  void dump(string filename) const;
  void dumpbooleq(const T val) const;
  void dumpboolgeq(const T val) const;
  void dumpboolleq(const T val) const;
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
  inline T dotProd ( const Array3D<T>& Vec ) const;
  inline T getNormSqr() const;
  inline T getNorm() const;  
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
inline bool Array3D<T>::onBndry(const size_t idx) const
{ // assumes input Array3D is a sign array, e.g. \pm 1, 0 are only values contained
  return((abs(getxp(idx)-get(idx)) + abs(getxm(idx)-get(idx)) + abs(getyp(idx)-get(idx)) + abs(getym(idx)-get(idx)) + abs(getzp(idx)-get(idx)) + abs(getzm(idx)-get(idx))) > 0);
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
inline vector<size_t> Array3D<T>::connectedComponent(const size_t idx, const double size) 
{ // calls connectedComponent with changeval==-1
  return(connectedComponent(idx,size,-1));
}

template <typename T>
vector<size_t> Array3D<T>::connectedComponent(const size_t idx, const double size, const int changeval) 
{ // returns a vector containing the index of all pixels in the connected component identified by q==q[idx] containing idx
  // Change values to changeval when touched to indicate they've been
  // checked; then change all the values back to the original
  // value at the end.
  // Assumes that changeval is not a value used by the array for any
  // other purpose.
  // size is an estimate of the proportion of the grid covered
  // by the connected component
  T value = get(idx);
  size_t curr = 0;
  vector<size_t> v;
  v.reserve(static_cast<size_t>(ceil(size)*N));
  v.push_back(idx);
  put(-1,idx); 
  
  while(curr < v.size())
  {
    if(get(xp(v[curr])) == value)
    { v.push_back(xp(v[curr])); put(-1,xp(v[curr]));}
    if(get(xm(v[curr])) == value)
    { v.push_back(xm(v[curr])); put(-1,xm(v[curr]));}
    if(get(yp(v[curr])) == value)
    { v.push_back(yp(v[curr])); put(-1,yp(v[curr]));}
    if(get(ym(v[curr])) == value)
    { v.push_back(ym(v[curr])); put(-1,ym(v[curr]));}
    if(get(zp(v[curr])) == value)
    { v.push_back(zp(v[curr])); put(-1,zp(v[curr]));}
    if(get(zm(v[curr])) == value)
    { v.push_back(zm(v[curr])); put(-1,zm(v[curr]));}
    curr++;
  }
  for(size_t ii=0; ii<v.size(); ++ii)
    put(value,v[ii]);

  return(v);
}

template <typename T>
vector<size_t> Array3D<T>::connectedComponentPM(const size_t idx, const double size) 
{ // returns a vector containing the index of all pixels in the connected component identified by q== \pm q[idx] containing idx
  // Change values to 0 when touched to indicate they've been
  // checked; then change all the values back to the original
  // value at the end.
  // Assumes that 0 is not a value used by the array for any
  // other purpose.
  // size is an estimate of the proportion of the grid covered
  // by the connected component
  T value = abs(get(idx));
  size_t curr = 0;
  vector<size_t> v;  // indices
  vector<bool> s; // sign (+:true -:false)
  v.reserve(static_cast<size_t>(ceil(size)*N));
  s.reserve(static_cast<size_t>(ceil(size)*N));
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
    if(abs(getzp(v[curr])) == value)
    { 
      v.push_back(zp(v[curr])); 
      s.push_back(getzp(v[curr])>0.0f);
      put(0,zp(v[curr]));
    }
    if(abs(getzm(v[curr])) == value)
    { 
      v.push_back(zm(v[curr])); 
      s.push_back(getzm(v[curr])>0.0f);
      put(0,zm(v[curr]));
    }
    curr++;
  }
  for(size_t ii=0; ii<v.size(); ++ii)
    if(s[ii])
      put(value,v[ii]);
    else
      put(-value,v[ii]);

  return(v);
}

template <typename T>
vector<int> Array3D<T>::connectedComponent(const int idx, const double size, const T compareval, bool(*compareFn)(const T, const T))  // this could be rewritten with Heap2
{ // compareType: 1 >, 2 >=, 3 ==, 4 <=, 5 <
  // labels connected component of data[idx] for which data[*] satisfies comparison with compareval

  // auxilary array to indicate which grid cells have been visited.
  Array3D<bool> touched(m,n,k);
  touched.fillWithValue(false);
  
  size_t curr = 0;
  vector<int> v;
  v.reserve(static_cast<int>(ceil(size)*N));

  if((*compareFn)(get(idx),compareval))
  {
    v.push_back(idx);
    touched.put(true,idx);
  }
  
  while(curr < v.size())
  {
    if(((*compareFn)(getxp(v[curr]),compareval)) && (!touched.getxp(v[curr])))
    { v.push_back(xp(v[curr])); touched.put(true,xp(v[curr]));}
    if(((*compareFn)(getxm(v[curr]),compareval)) && (!touched.getxm(v[curr])))
    { v.push_back(xm(v[curr])); touched.put(true,xm(v[curr]));}
    if(((*compareFn)(getyp(v[curr]),compareval)) && (!touched.getyp(v[curr])))
    { v.push_back(yp(v[curr])); touched.put(true,yp(v[curr]));}
    if(((*compareFn)(getym(v[curr]),compareval)) && (!touched.getym(v[curr])))
    { v.push_back(ym(v[curr])); touched.put(true,ym(v[curr]));}
    if(((*compareFn)(getzp(v[curr]),compareval)) && (!touched.getzp(v[curr])))
    { v.push_back(zp(v[curr])); touched.put(true,zp(v[curr]));}
    if(((*compareFn)(getzm(v[curr]),compareval)) && (!touched.getzm(v[curr])))
    { v.push_back(zm(v[curr])); touched.put(true,zm(v[curr]));}
    curr++;
  }

  return(v);
}

template <typename T>
inline bool Array3D<T>::compareG(const T currval, const T compareval) 
{
  return(currval > compareval);
}

template <typename T>
inline bool Array3D<T>::compareGeq(const T currval, const T compareval) 
{
  return(currval >= compareval);
}

template <typename T>
inline bool Array3D<T>::compareEq(const T currval, const T compareval)
{
  return(currval == compareval);
}

template <typename T>
inline bool Array3D<T>::compareLeq(const T currval, const T compareval) 
{
  return(currval <= compareval);
}

template <typename T>
inline bool Array3D<T>::compareL(const T currval, const T compareval) 
{
  return(currval < compareval);
}

template <typename T>
vector<size_t> Array3D<T>::boundaryOfCC(const vector<size_t> cc)
{
  T value = get(cc[0]);
  vector<size_t> v;
  v.reserve(cc.size());
  for(size_t ii=0;ii<cc.size();++ii)
    if(get(xp(cc[ii]))!=value)
      v.push_back(cc[ii]);
    else if(get(xm(cc[ii]))!=value)
      v.push_back(cc[ii]);
    else if(get(yp(cc[ii]))!=value)
      v.push_back(cc[ii]);
    else if(get(ym(cc[ii]))!=value)
      v.push_back(cc[ii]);
    else if(get(zp(cc[ii]))!=value)
      v.push_back(cc[ii]);
    else if(get(zm(cc[ii]))!=value)
      v.push_back(cc[ii]);
  return(v);
}

template <typename T>
size_t Array3D<T>::diff1Nbors(T val, size_t idx)
{ // count the number of 6-neighbors of idx different than val
  size_t diff = 0;
  diff += (val != get(xp(idx))) + (val != get(xm(idx))) 
        + (val != get(yp(idx))) + (val != get(ym(idx)))
        + (val != get(zp(idx))) + (val != get(zm(idx)));
  return(diff);
}

template <typename T>
size_t Array3D<T>::diff2Nbors(T val, size_t idx)
{ // count the number of 2nd-nearest-neighbors of idx different than val
  size_t diff = 0;
  diff += (val != get(xp(yp(idx)))) + (val != get(xm(yp(idx)))) 
        + (val != get(xp(ym(idx)))) + (val != get(xm(ym(idx))))
        + (val != get(zp(xp(idx)))) + (val != get(zm(xp(idx))))
        + (val != get(zp(xm(idx)))) + (val != get(zm(xm(idx))))
        + (val != get(zp(yp(idx)))) + (val != get(zm(yp(idx))))
        + (val != get(zp(ym(idx)))) + (val != get(zm(ym(idx))));
  return(diff);
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
void Array3D<T>::dump() const
{
  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.setf(std::ios::showpos);
  cout.precision(5);
  for(size_t kk=0; kk < k; ++kk)
  {
    for(size_t ii=0; ii < m; ++ii)
    {
      for(size_t jj=0; jj < n; ++jj)
        cout << get(ii,jj,kk) << " ";
      cout << endl;
    }
    cout << endl;
  }
  cout.unsetf(std::ios::fixed);
  cout.unsetf(std::ios::floatfield);
  cout.unsetf(std::ios::showpos);
}

template <typename T>
void Array3D<T>::dump(string filename) const
{
  ofstream of;
  of.open(filename.c_str());
  
  of.setf(std::ios::fixed,std::ios::floatfield);
  of.setf(std::ios::showpos);
  of.precision(5);
  for(size_t kk=0; kk < k; ++kk)
    for(size_t ii=0; ii < m; ++ii)
    {
      for(size_t jj=0; jj < n; ++jj)
        of << get(ii,jj,kk) << " ";
      of << endl;
    }
  of.unsetf(std::ios::fixed);
  of.unsetf(std::ios::floatfield);
  of.unsetf(std::ios::showpos);
  of.close();
  of.clear();
}

template <typename T>
void Array3D<T>::dumpbooleq(const T val) const
{ // outputs 1 
  for(size_t kk=0; kk < k; ++kk)
  {
    for(size_t ii=0; ii < m; ++ii)
    {
      for(size_t jj=0; jj < n; ++jj)
        cout << (get(ii,jj,kk)==val);
      cout << endl;
    }
    cout << endl;
  }
}

template <typename T>
void Array3D<T>::dumpboolgeq(const T val) const
{
  for(size_t kk=0; kk < k; ++kk)
  {
    for(size_t ii=0; ii < m; ++ii)
    {
      for(size_t jj=0; jj < n; ++jj)
        cout << (get(ii,jj,kk)>=val);
      cout << endl;
    }
    cout << endl;
  }
}

template <typename T>
void Array3D<T>::dumpboolleq(const T val) const
{
  for(size_t kk=0; kk < k; ++kk)
  {
    for(size_t ii=0; ii < m; ++ii)
    {
      for(size_t jj=0; jj < n; ++jj)
        cout << (get(ii,jj,kk)<=val);
      cout << endl;
    }
    cout << endl;
  }
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

template <typename T>
inline T Array3D<T>::dotProd ( const Array3D<T>& Vec ) const
{
  T val = 0;
  for ( size_t i = 0 ; i < N ; ++i )
    val += data[i] * Vec.data[i];
  return(val);
}

template <typename T>
inline T Array3D<T>::getNormSqr() const
{
  T val = 0;
  for ( size_t i = 0; i < N; ++i )
    val += data[i] * data[i];
  return(val);
}

template <typename T>
inline T Array3D<T>::getNorm() const
{
  return(sqrt(getNormSqr()));
}


#endif
