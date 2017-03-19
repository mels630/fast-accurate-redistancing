#ifndef _REDIST_UNSIGNED_H_
#define _REDIST_UNSIGNED_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array2d.hpp"
#include "idarray2d.hpp"
#include "heap.hpp"
#include "toolbox.hpp"
#include "defs.h"
#include "assert.h"

#include "mex.h"

class RedistU
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
  void thresholdAwayFromInterface();
  void applyResult(const struct helt &h);
  explicit RedistU(); // no empty constructor
  void findNborDirections(const double (&x0)[2], const double (&xc)[2], double (&xl)[2], double (&xr)[2], const double delta);
  struct helt performDOSurf_unsigned(const int idx);
  struct helt performDO_unsigned(const int idx, const double cpxguess, const double cpyguess);
  double search1D_unsigned(const int idx, double (&x)[2]);
  bool bracket_unsigned(const int idx, double (&cpguess)[2], double (&x1)[3], double (&x2)[3]);
  double bisect_unsigned(double (&result)[2], double (&x1)[3], double (&x2)[3]);
  inline double dirDerivative(double *grad, double *dir);
  void setInterfaceValuesDO_unsigned(const double cutoff);
  void lineSearch_unsigned(const int idx, double (&grad)[2], double (&cpguess)[2]);
  
public:
  explicit RedistU(const Array2D<double> &_u, const int _width, const int _flag);
  explicit RedistU(const Array2D<double> &_u, const double _dx, const double _dy, const int _width, const int _flag);
  void redistance();
  void dump_u(double *v);
  inline Array2D<double> dump_u();

  void dump_cp(double *cpxArr, double *cpyArr);
  
};

inline Array2D<double> RedistU::dump_u()
{ // assumes return argument already has correct size
  return(u);
}

inline double RedistU::dirDerivative(double *grad, double *dir)
{ // assume grad and dir each hold at least 2 doubles
  return(grad[0]*dir[0]+grad[1]*dir[1]);
}

#endif
