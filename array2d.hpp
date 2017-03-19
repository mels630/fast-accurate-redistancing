#ifndef _ARRAY2D_HPP_
#define _ARRAY2D_HPP_

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
class Array2D
{

protected:
  vector<T> data; // contains (m x n) = N elements of type T
  const int m;          // number of elements in the y-direction
  const int n;          // number of elements in the x-direction
  const int N;          // total number of elements
  const double dx;      // grid spacing in x
  const double dy;      // grid spacing in y
  explicit Array2D<T>(); // don't use empty constructor

  inline int sub2ind(const int ii, const int jj) const;
  void transposeData();
  void flipData(const int ii, const int jj);

  void labelHexCC(const T value, const bool sign, const int idx, Array2D<int> &label) const;
  void labelHexCCX(const int idx, const int val, vector<int> &changed);
  void labelHexCCY(const int idx, const int val, vector<int> &changed);

public: 
  explicit Array2D<T>(const int nn);
  explicit Array2D<T>(const int mm, const int nn);
  explicit Array2D<T>(const T* indata, const int mm, const int nn);
  explicit Array2D<T>(const int mm, const int nn, const double _dx, const double _dy);
  explicit Array2D<T>(const T* indata, const int mm, const int nn, const double _dx, const double _dy);
  explicit Array2D<T>(const Array2D<T> &Input, const int CopyFlag);
  Array2D<T>(const Array2D<T> &Input); 
  inline T get(const int idx) const;
  inline T get(const int ii, const int jj) const;
  inline bool onBndry(const int idx) const;
  inline int getm() const;
  inline int getn() const;
  inline int getN() const;
  inline double getdx() const;
  inline double getdy() const;
  inline double lenx() const;
  inline double leny() const;
  inline int size() const;
  inline void put(const T value, const int idx);
  inline void put(const T value, const int ii, const int jj);
  inline void putadd(const T value, const int idx);
  inline void putadd(const T value, const int ii, const int jj);
  void copyData(const Array2D<T> &input);
  inline vector<int> connectedComponent(const int idx, const double size);
  vector<int> connectedComponent(const int idx, const double size, const int changeval);
  vector<int> connectedComponentPM(const int idx, const double size);
  vector<int> boundaryOfCC(const vector<int> cc);
  int diff1Nbors(int val, int idx);
  int diff2Nbors(int val, int idx);
  inline T maxval() const;
  inline T minval() const;
  inline T minabsval() const;
  inline void fillWithValue(const T value);
  inline T fourNborMin(const int idx) const;
  inline T fourNborMax(const int idx) const;
  void dump() const;
  void dump(string filename) const;
  void dumpbooleq(const T val) const;
  void dumpboolgeq(const T val) const;
  void dumpboolleq(const T val) const;
  T* dataAddress();
  vector<T> returnData() const;
  T* returnDataArray() const;
  Array2D<T> duplicateArray2D() const;
  inline double getX(const int idx) const;
  inline double getY(const int idx) const;
  inline int getXidx(const int idx) const;
  inline int getYidx(const int idx) const;
  inline int xp(const int idx) const;
  inline int xm(const int idx) const;
  inline int yp(const int idx) const;
  inline int ym(const int idx) const;
  inline T getxp(const int idx) const;
  inline T getxm(const int idx) const;
  inline T getyp(const int idx) const;
  inline T getym(const int idx) const;

  inline int hexN2(const int idx, const int nbor) const;
  inline int hexN(const int idx, const int nbor) const;
  inline int hexNNP2(const int idx, const int nbor) const;
  inline int hexNNP(const int idx, const int nbor) const;
  Array2D<int> labelHex(const T value, const bool sign) const;
  bool labelHexFromIdxX(const int idx);
  bool labelHexFromIdxY(const int idx);

  int countVal(const T value) const;
  
  // algebraic operations
  inline T& operator[] ( int I );
  inline const T& operator[] ( int I ) const;
  inline Array2D<T>& operator= ( const Array2D<T>& Vec );
  inline Array2D<T>& operator*= ( const T Value );
  inline Array2D<T>& operator/= ( const T Value );
  inline Array2D<T>& operator+= ( const Array2D<T> &Vec );
  inline Array2D<T>& operator-= ( const Array2D<T> &Vec );
  inline void setAll( const T val );
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

template <typename T>
inline int Array2D<T>::sub2ind(const int ii, const int jj) const
{
  assert((ii+m*jj) < N);
  return(ii+m*jj);
}

template <typename T> 
Array2D<T>::Array2D(const int nn) :
  m(nn),
  n(nn),
  N(nn*nn),
  dx(1.0f/static_cast<double>(n)),
  dy(1.0f/static_cast<double>(m))
{
  if(N <= 0)
  {
    cout << "N = " << N << ", should be positive. Aborting ..." << endl;
    abort();
  }
  data.resize(N);
}

template <typename T> 
Array2D<T>::Array2D(const int mm, const int nn) :
  m(mm),
  n(nn),
  N(mm*nn),
  dx(1.0f/static_cast<double>(n)),
  dy(1.0f/static_cast<double>(m))
{
  if(N <= 0)
  {
    cout << "N = " << N << ", should be positive. Aborting ..." << endl;
    abort();
  }
  data.resize(N);
}

template <typename T>
Array2D<T>::Array2D(const T* indata, const int mm, const int nn) :
  m(mm),
  n(nn),
  N(mm*nn),
  dx(1.0f/static_cast<double>(n)),
  dy(1.0f/static_cast<double>(m))
{
  if(N <= 0)
  {
    cout << "N = " << N << ", should be positive. Aborting ..." << endl;
    abort();
  }
  data.resize(N);
  for(int ii=0; ii<N; ++ii) // deep copy
    data[ii] = indata[ii];
}

template <typename T> 
Array2D<T>::Array2D(const int mm, const int nn, const double _dx, const double _dy) :
  m(mm),
  n(nn),
  N(mm*nn),
  dx(_dx),
  dy(_dy)
{
  if(N <= 0)
  {
    cout << "N = " << N << ", should be positive. Aborting ..." << endl;
    abort();
  }
  data.resize(N);
}

template <typename T>
Array2D<T>::Array2D(const T* indata, const int mm, const int nn, const double _dx, const double _dy) :
  m(mm),
  n(nn),
  N(mm*nn),
  dx(_dx),
  dy(_dy)
{
  if(N <= 0)
  {
    cout << "N = " << N << ", should be positive. Aborting ..." << endl;
    abort();
  }
  data.resize(N);
  for(int ii=0; ii<N; ++ii) // deep copy
    data[ii] = indata[ii];
}

template <typename T>
Array2D<T>::Array2D( const Array2D<T> &Input, const int CopyFlag) :
  m(Input.m),
  n(Input.n),
  N(Input.N),
  dx(Input.dx),
  dy(Input.dy)
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
Array2D<T>::Array2D( const Array2D<T> &Input) :
  m(Input.m),
  n(Input.n),
  N(Input.N),
  dx(Input.dx),
  dy(Input.dy)
{
  data = Input.data;
}

template <typename T>
inline int Array2D<T>::xp(const int idx) const
{
  return((idx)<(N-m) ? (idx +m) : (idx+m-N));
}

template <typename T>
inline int Array2D<T>::xm(const int idx) const
{
  return((idx)>(m-1) ? (idx -m) : (idx-m+N));
}

template <typename T>
inline int Array2D<T>::yp(const int idx) const
{
  return( (idx%m)==(m-1) ? (idx-m+1) : (idx+1));
}

template <typename T>
inline int Array2D<T>::ym(const int idx) const
{
  return( (idx%m)==(0) ? (idx+m-1) : (idx-1));
}

template <typename T>
inline T Array2D<T>::getxp(const int idx) const
{
  return(data[xp(idx)]);
}

template <typename T>
inline T Array2D<T>::getxm(const int idx) const
{
  return(data[xm(idx)]);
}

template <typename T>
inline T Array2D<T>::getyp(const int idx) const
{
  return(data[yp(idx)]);
}

template <typename T>
inline T Array2D<T>::getym(const int idx) const
{
  return(data[ym(idx)]);
}

template <typename T>
inline int Array2D<T>::getm() const
{
  return(m);
}

template <typename T>
inline int Array2D<T>::getn() const
{
  return(n);
}

template <typename T>
inline int Array2D<T>::getN() const
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
inline int Array2D<T>::size() const
{
  return(N);
}

template <typename T>
inline T Array2D<T>::get(const int idx) const
{
  return(data[idx]);
}

template <typename T>
inline T Array2D<T>::get(const int ii, const int jj) const
{
  return(data[sub2ind(ii,jj)]);
}

template <typename T>
inline bool Array2D<T>::onBndry(const int idx) const
{ // assumes input Array2D is a sign array, e.g. \pm 1, 0 are only values contained
  return((abs(getxp(idx)-get(idx)) + abs(getxm(idx)-get(idx)) + abs(getyp(idx)-get(idx)) + abs(getym(idx)-get(idx))) > 0);
}

template <typename T>
inline void Array2D<T>::put(const T value, const int idx)
{
  data[idx] = value;
}

template <typename T>
inline void Array2D<T>::put(const T value, const int ii, const int jj)
{
  data[sub2ind(ii,jj)] = value;
}

template <typename T>
inline void Array2D<T>::putadd(const T value, const int idx)
{
  data[idx] += value;
}

template <typename T>
inline void Array2D<T>::putadd(const T value, const int ii, const int jj)
{
  data[sub2ind(ii,jj)] += value;
}

template <typename T>
void Array2D<T>::copyData(const Array2D<T> &input)
{ // copy data after checking for agreement of sizes
  if((input.getm() != m) || (input.getn() != n))
    cout << "Data does not agree in size. Skipping call to copyData" << endl;
  else
    data = input.returnData();
}

template <typename T>
inline vector<int> Array2D<T>::connectedComponent(const int idx, const double size) 
{ // calls connectedComponent with changeval==-1
  return(connectedComponent(idx,size,-1));
}


template <typename T>
vector<int> Array2D<T>::connectedComponent(const int idx, const double size, const int changeval) 
{ // returns a vector containing the index of all pixels in the connected component identified by q==q[idx] containing idx
  // Change values to changeval when touched to indicate they've been
  // checked; then change all the values back to the original
  // value at the end.
  // Assumes that changeval is not a value used by the array for any
  // other purpose.
  // size is an estimate of the proportion of the grid covered
  // by the connected component
  int value = get(idx);
  int curr = 0;
  vector<int> v;
  v.reserve(static_cast<int>(ceil(size)*N));
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
    curr++;
  }
  for(int ii=0; ii<v.size(); ++ii)
    put(value,v[ii]);

  return(v);
}

template <typename T>
vector<int> Array2D<T>::connectedComponentPM(const int idx, const double size) 
{ // returns a vector containing the index of all pixels in the connected component identified by q== \pm q[idx] containing idx
  // Change values to 0 when touched to indicate they've been
  // checked; then change all the values back to the original
  // value at the end.
  // Assumes that 0 is not a value used by the array for any
  // other purpose.
  // size is an estimate of the proportion of the grid covered
  // by the connected component
  int value = abs(get(idx));
  int curr = 0;
  vector<int> v;  // indices
  vector<bool> s; // sign (+:true -:false)
  v.reserve(static_cast<int>(ceil(size)*N));
  s.reserve(static_cast<int>(ceil(size)*N));
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
  for(int ii=0; ii<v.size(); ++ii)
    if(s[ii])
      put(value,v[ii]);
    else
      put(-value,v[ii]);

  return(v);
}

template <typename T>
vector<int> Array2D<T>::boundaryOfCC(const vector<int> cc)
{
  int value = get(cc[0]);
  vector<int> v;
  v.reserve(cc.size());
  for(int ii=0;ii<cc.size();++ii)
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
  T mval = data[0];
  for(int ii=1;ii<N;++ii)
    mval = mymax(mval,data[ii]);
  return(mval);
}

template <typename T>
inline T Array2D<T>::minval() const
{
  assert(N>0);
  T mval = data[0];
  for(int ii=1;ii<N;++ii)
    mval = mymin(mval,data[ii]);
  return(mval);
}

template <typename T>
inline T Array2D<T>::minabsval() const
{
  if(N == 0)
    return(0);
  T mval = data[0]*mysign(data[0]);
  for(int ii=1;ii<N;++ii)
    mval = mymin(mval,data[ii]*mysign(data[ii]));
  return(mval);
}

template <typename T>
inline void Array2D<T>::fillWithValue(const T value)
{
  for(int ii=0; ii<N; ++ii)
    data[ii] = value;
}

template <typename T>
inline T Array2D<T>::fourNborMin(const int idx) const
{
  return(mymin(mymin(data[xp(idx)],data[xm(idx)]),mymin(data[yp(idx)],data[ym(idx)])));
}

template <typename T>
inline T Array2D<T>::fourNborMax(const int idx) const
{
  return(mymax(mymax(data[xp(idx)],data[xm(idx)]),mymax(data[yp(idx)],data[ym(idx)])));
}

template <typename T>
void Array2D<T>::dump() const
{
  cout.setf(std::ios::fixed,std::ios::floatfield);
  cout.setf(std::ios::showpos);
  cout.precision(5);
  for(int ii=0; ii < m; ++ii)
  {
    for(int jj=0; jj < n; ++jj)
      cout << get(ii,jj) << " ";
    cout << endl;
  }
  cout.unsetf(std::ios::fixed);
  cout.unsetf(std::ios::floatfield);
  cout.unsetf(std::ios::showpos);
}

template <typename T>
void Array2D<T>::dump(string filename) const
{
  ofstream of;
  of.open(filename.c_str());
  
  of.setf(std::ios::fixed,std::ios::floatfield);
  of.setf(std::ios::showpos);
  of.precision(5);
  for(int ii=0; ii < m; ++ii)
  {
    for(int jj=0; jj < n; ++jj)
      of << get(ii,jj) << " ";
    of << endl;
  }
  of.unsetf(std::ios::fixed);
  of.unsetf(std::ios::floatfield);
  of.unsetf(std::ios::showpos);
  of.close();
  of.clear();
}

template <typename T>
void Array2D<T>::dumpbooleq(const T val) const
{ // outputs 1 
  for(int ii=0; ii < m; ++ii)
  {
    for(int jj=0; jj < n; ++jj)
      cout << (get(ii,jj)==val);
    cout << endl;
  }
}

template <typename T>
void Array2D<T>::dumpboolgeq(const T val) const
{
  for(int ii=0; ii < m; ++ii)
  {
    for(int jj=0; jj < n; ++jj)
      cout << (get(ii,jj)>=val);
    cout << endl;
  }
}

template <typename T>
void Array2D<T>::dumpboolleq(const T val) const
{
  for(int ii=0; ii < m; ++ii)
  {
    for(int jj=0; jj < n; ++jj)
      cout << (get(ii,jj)<=val);
    cout << endl;
  }
}

template <typename T>
T* Array2D<T>::dataAddress()
{ // returns address to first data element
  return(&(data.front()));
}

template <typename T>
vector<T> Array2D<T>::returnData() const
{ // returns a copy of the data vector
  return(data);
}

template <typename T>
T* Array2D<T>::returnDataArray() const
{ // returns a pointer to a new (independent) array with data copied in.
  T* newArray = new T[getN()];
  for(int ii=0;ii<N;++ii)
    newArray[ii] = get(ii);
  return(newArray);
}


template <typename T>
Array2D<T> Array2D<T>::duplicateArray2D() const
{ // deep copy of Array2D data
  Array2D<T> v(m,n,dx,dy);
  v.data = this->data;
  return(v);
}

template <typename T>
int Array2D<T>::diff1Nbors(int val, int idx)
{ // count the number of 4-neighbors of idx different than val
  int diff = 0;
  diff += (val != get(xp(idx))) + (val != get(xm(idx))) 
           + (val != get(yp(idx))) + (val != get(ym(idx)));
  return(diff);
}

template <typename T>
int Array2D<T>::diff2Nbors(int val, int idx)
{ // count the number of 2nd-nearest-neighbors of idx different than val
  int diff = 0;
  diff += (val != get(xp(yp(idx)))) + (val != get(xm(yp(idx)))) 
          + (val != get(xp(ym(idx)))) + (val != get(xm(ym(idx))));
  return(diff);
}

template <typename T>
void Array2D<T>::transposeData()
{
  for(int ii=0;ii<n;++ii)
    for(int jj=0;jj<n;++jj)
      flipData(ii,jj);
}

template <typename T>
void Array2D<T>::flipData(const int ii, const int jj)
{
  T tmp = data[ii+n*jj];
  data[ii+n*jj] = data[jj+n*ii];
  data[jj+n*ii] = tmp;
}

template <typename T>
inline double Array2D<T>::getX(const int idx) const
{
  return(static_cast<double>(getXidx(idx)) * dx);
}

template <typename T>
inline double Array2D<T>::getY(const int idx) const
{
  return(static_cast<double>(getYidx(idx)) * dy);
}

template <typename T>
inline int Array2D<T>::getXidx(const int idx) const
{
  return(idx/m);
}

template <typename T>
inline int Array2D<T>::getYidx(const int idx) const
{
  return(idx%m);
}

template <typename T>
inline int Array2D<T>::hexN2(const int idx, const int nbor) const
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
    cout << "hexN2 only accepts values between 1 and 6 in the second argument." << endl;
    return(-1);
  }
}

template <typename T>
inline int Array2D<T>::hexN(const int idx, const int nbor) const
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
    cout << "hexN only accepts values between 1 and 6 in the second argument." << endl;
    return(-1);
  }
}

template <typename T>
inline int Array2D<T>::hexNNP2(const int idx, const int nbor) const
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
inline int Array2D<T>::hexNNP(const int idx, const int nbor) const
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
Array2D<int> Array2D<T>::labelHex(const T value, const bool sign) const
{ // return labeling of the connected components of data>=value (if sign==true)
  // and data<value (if sign==false), with hexagonal connectivity
  Array2D<int> label(m,n);
  label.fillWithValue(-1);

  int labelnum=0;
  for(int ii=0;ii<N;++ii)
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
void Array2D<T>::labelHexCC(const T value, const bool sign, const int idx, Array2D<int> &label) const
{ // recursively add pieces of this connected component

  stack<int> s;
  s.push(idx);
  int stackval;
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
bool Array2D<T>::labelHexFromIdxY(const int idx)
{ // Try to label from idx without crossing Y=0, and see if the component is "infinite" from top to bottom
  vector<int> changed;
  const int val = data[idx];

  for(int ii=0; ii<n; ++ii)
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
  for(int ii=0; ii<n; ++ii)
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

  //if(connection)
  //{
  //  cout << "Connection in Y found for " << val << endl;
  //  dumpbool(0);
  //}

  for(int ii=0; ii<changed.size(); ++ii)
    data[changed[ii]] = val;

  return(connection);
}

template <typename T>
bool Array2D<T>::labelHexFromIdxX(const int idx)
{ // Try to label from idx without crossing X=0, and see if the component is "infinite" from left to right
  vector<int> changed;
  const int val = data[idx];

  //if(idx==0)
  //{
  //  int zv = 0;
  //  for(int ii=0;ii<N;++ii)
  //    zv += (data[ii]==0);
  //  cout << "Prelooking for connection at 0, there are " << zv << " zero values in data." << endl;
  //}

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
  for(int ii=0; ii<m; ++ii)
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

  //if(connection && (idx==0))
  //{
  //  cout << "Connection in X found for " << val << endl;
  //  dumpbool(0);
  //}

  // change values back
  for(int ii=0; ii<changed.size(); ++ii)
    data[changed[ii]] = val;
  return(connection);
}

template <typename T>
void Array2D<T>::labelHexCCY(const int idx, const int val, vector<int> &changed)
{ // recursively add pieces of this connected component
  stack<int> s;
  s.push(idx);
  int stackval;
  while(!s.empty())
  {
    stackval = s.top();
    s.pop();
    if((stackval % m) != 0) // check all 6 neighbors
      for(int ii=1; ii<=6; ++ii)
      {
        int cstackval = hexN(stackval,ii);
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
        int cstackval = hexN(stackval,ii);
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
void Array2D<T>::labelHexCCX(const int idx, const int val, vector<int> &changed)
{ // recursively add pieces of this connected component
  stack<int> s;
  s.push(idx);
  int stackval;
  while(!s.empty())
  {
    stackval = s.top();
    s.pop();
    if((stackval / m) != 0) // check all 6 neighbors
      for(int ii=1; ii<=6; ++ii)
      {
        int cstackval = hexN(stackval,ii);
        if(data[cstackval] == val)
        {
          //if(cstackval == (N-1))
          //  cout << "HEY1" << endl;
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
        int cstackval = hexN(stackval,jj[ii]);
        if(data[cstackval] == val)
        {
          //if(cstackval == (N-1))
          //  cout << "HEY2" << endl;
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
          int cstackval = hexN(stackval,jj[ii]);
          if(data[cstackval] == val)
          {
            //if(cstackval == (N-1))
            //  cout << "HEY3, stackval = " << stackval << ", jj[ii] = " << jj[ii] << ", ii = " << ii << endl;
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
int Array2D<T>::countVal(const T value) const
{
  int count = 0;
  for(int ii=0; ii<N; ++ii)
    if(data[ii] == value)
      count++;
  return(count);
}

template <typename T>
inline T& Array2D<T>::operator[] ( int I )
{
  return data[I];
}

template <typename T>
inline const T& Array2D<T>::operator[] ( int I ) const
{
  return data[I];
}

template <typename T>
inline Array2D<T>& Array2D<T>::operator= ( const Array2D<T>& Vec )
{
  assert(m == Vec.getm());
  assert(n == Vec.getn());
  data = Vec.data;
  return *this;
}

template <typename T>
inline Array2D<T>& Array2D<T>::operator*= ( const T Value )
{
  for ( int i = 0; i < N; ++i )
    data[i] *= Value;
  return *this;
}

template <typename T>
inline Array2D<T>& Array2D<T>::operator/= ( const T Value )
{
  for ( int i = 0; i < N; ++i )
    data[i] /= Value;
  return *this;
}

template <typename T>
inline Array2D<T>& Array2D<T>::operator+= ( const Array2D<T> &Vec )
{
  for ( int i = 0; i < N; ++i )
    data[i] += Vec.get( i );
  return *this;
}

template <typename T>
inline Array2D<T>& Array2D<T>::operator-= ( const Array2D<T> &Vec )
{
  for ( int i = 0; i < N; ++i )
    data[i] -= Vec.get( i );
  return *this;
}

template <typename T>
inline void Array2D<T>::setAll( const T val )
{
  for ( int i = 0; i < N; ++i )
    data[i] = val;
}

template <typename T>
inline void Array2D<T>::addMultiple ( const Array2D<T>& Vec, T Factor )
{
  for ( int i = 0; i < N; ++i )
    data[i] += Vec.get( i ) * Factor;
}

template <typename T>
inline T Array2D<T>::dotProd ( const Array2D<T>& Vec ) const
{
  T val = 0;
  for ( int i = 0 ; i < N ; ++i )
    val += data[i] * Vec.data[i];
  return val;
}

template <typename T>
inline T Array2D<T>::getNormSqr() const
{
  T val = 0;
  for ( int i = 0; i < N; ++i )
    val += data[i] * data[i];
  return val;
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
  for ( int i = 0; i < m; ++i ) {
    for ( int j = 0; j < n; ++j )
      file << get( i, j ) << " ";
    file << endl;
  }
  file.close();
}

template <typename T>
inline void Array2D<T>::loadASCII ( const char *FileName ) {
  std::ifstream file( FileName );
  for ( int i = 0; i < m; ++i )
    for ( int j = 0; j < n; ++j ) {
      T value;
      file >> value;
      put( value, i, j );
    }
  file.close();
}

#endif
