#ifndef _REDIST3C_H_
#define _REDIST3C_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array3d.hpp"
#include "idarray3d.hpp"
#include "heap.hpp"
//#include "toolbox.hpp"
#include "toolbox3d.hpp"
#include "defs.h"
#include "assert.h"

// Later, it would be even better to have a SINGLE redist class
// which would determine which redist is called simply based on 
// initializer (Array2D or Array3D)

class Redist3
{ 
private:
  const int width;
  const int flag;
  const int m, n, k, N;
  const double dx, dy, dz;
  const double thres;
  const bool cpflag;
  Heap h;
  Array3D<bool> state;
  vector<int> bnd;

  void fastMarchingRedist();
  void directionalOptimization();
  void updateAndAddNeighborsToHeap(const int idx);
  void updateAndAddNeighborsToHeapDO(const int idx);
  double estimateUpdate(const int idx);
  void sort(double &a1, double &a2, double &a3);
  void setInterfaceValues();
  void setInterfaceValuesDO();
  void thresholdAwayFromInterface(vector<int> &bndry);
  void applyResult(const struct helt &h);
  struct helt performDO(const int idx);
  struct helt performDO(const int idx, const double cpxguess, const double cpyguess, const double cpzguess);
  void lineSearch(const int idx, double (&grad)[3], double (&cpguess)[3]);
  bool findOppSign(const int idx, double (&guess)[3]);
  void bisect(const int idx, double (&guess)[3]);
  double bisect(double (&result)[3], double (&xm)[4], double (&xp)[4]);
  explicit Redist3(); // no empty constructor
  bool bracket(const int idx, double (&cpguess)[3], double (&xm)[4], double (&xp)[4]);
  double findNborDirections(const double (&x0)[3], double (&xx)[5][3], const double delta, double (&rvec)[3], double (&v1)[3], double (&v2)[3]);
  double search1D(const int idx, double (&x)[3]);
  void secondOrderIterations();
  bool diffSign(int idx);
  
public:
  explicit Redist3(const Array3D<double> &_u, const int _width, const int _flag);
  void redistance();
  void dump_u(double *v);
  void dump_cp(double *cpx_d, double *cpy_d, double *cpz_d);

  IDArray3D u0;
  Array3D<double> u;
  Array3D<double> cpx, cpy, cpz;
};

#endif
