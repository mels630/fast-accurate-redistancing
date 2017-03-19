#ifndef _REDISTC_H_
#define _REDISTC_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array2d.hpp"
#include "idarray2d.hpp"
#include "heap.hpp"
#include "toolbox.hpp"
#include "defs.h"
#include "assert.h"

//#include "mex.h"

class Redist
{ 
private:
  const int width;
  const int flag;
  const int m;
  const int n;
  const int N;
  const double dx;
  const double dy;
  const double thres;
  const bool cpflag;
  Heap h;
  Array2D<bool> state;
  vector<int> bnd;
  IDArray2D u0;
  Array2D<double> u;
  Array2D<double> cpx, cpy;

  void fastMarchingRedist();
  void directionalOptimization();
  void updateAndAddNeighborsToHeap(const int idx);
  void updateAndAddNeighborsToHeapDO(const int idx);
  double estimateUpdate(const int idx);
  void setInterfaceValues();
  void setInterfaceValuesDO();
  void thresholdAwayFromInterface();
  void applyResult(const struct helt &h);
  struct helt performDO(const int idx);
  struct helt performDOSurf(const int idx);
  struct helt performDO(const int idx, const double cpxguess, const double cpyguess);
  void lineSearch(const int idx, double (&grad)[2], double (&cpguess)[2]);
  bool findOppSign(const int idx, double (&guess)[2]);
  void bisect(const int idx, double (&guess)[2]);
  double bisect(double (&result)[2], double (&xm)[3], double (&xp)[3]);
  explicit Redist(); // no empty constructor
  bool bracket(const int idx, double (&cpguess)[2], double (&xm)[3], double (&xp)[3]);
  void findNborDirections(const double (&x0)[2], const double (&xc)[2], double (&xl)[2], double (&xr)[2], const double delta);
  double search1D(const int idx, double (&x)[2]);
  void secondOrderIterations();
  bool diffSign(int idx);
  
public:
  explicit Redist(const Array2D<double> &_u, const int _width, const int _flag);
  explicit Redist(const Array2D<double> &_u, const double _dx, const double _dy, const int _width, const int _flag);
  void redistance();
  void dump_u(double *v);
  inline Array2D<double> dump_u();

  void dump_cp(double *cpxArr, double *cpyArr);
  
};

inline Array2D<double> Redist::dump_u()
{ // assumes return argument already has correct size
  return(u);
}

#endif
