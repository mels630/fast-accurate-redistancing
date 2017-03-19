#ifndef _REINIT_H_
#define _REINIT_H_

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

class Reinit
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
  Heap h;
  Array2D<bool> state;
  vector<int> bnd;
  //Array2D<double> tpx;
  //Array2D<double> tpy;
  IDArray2D cpx;
  IDArray2D cpy;
  IDArray2D dd;

  explicit Reinit(); // no empty constructor
  void tempInterfaceValues();
  void applyResult(const struct helt &h);
  void updateAndAddNeighborsToHeapDO(const int idx);
  void thresholdAwayFromInterface();
  double distToCP(const int idx, const double cpxidx, const double cpyidx) const;
  void findNborDirections(const struct point x0, const struct point xc, struct point &xl, struct point &xr, const double delta) const;
  double search1D(const int idx, struct point &x) const;
  struct helt performDO(const int idx, const struct point cpguess);
  bool bracket(const struct point x0, const struct point d, const double initT, double (&t)[2], double (&htout)[2],double (&hderout)[2]) const;
  bool bracketOld(const struct point x0, const struct point d, const double initT, double (&t)[2], double (&htout)[2],double (&hderout)[2]) const;
  void computeHtAndDer(const struct point x, const struct point d, const double t, double &ht, double &htder) const;
  void computeHtAndDerOld(const struct point x, const struct point d, const double t, double &ht, double &htder) const;
  double bisect(const struct point &x0, const struct point &d, const double (&t)[2], const double (&ht)[2], const double (&htder)[2], struct point &result) const;

  /*
  void fastMarchingReinit();
  Array2D<double>* directionalOptimization();
  void updateAndAddNeighborsToHeap(const int idx);
  double estimateUpdate(const int idx);
  void setInterfaceValues();
  void setInterfaceValuesDO();
  struct helt performDO(const int idx);
  void lineSearch(const int idx, double (&grad)[2], double (&cpguess)[2]);
  bool findOppSign(const int idx, double (&guess)[2]);
  void bisect(const int idx, double (&guess)[2]);
  double bisect(double (&result)[2], double (&xm)[3], double (&xp)[3]);
  bool bracket(const int idx, double (&cpguess)[2], double (&xm)[3], double (&xp)[3]);
  void secondOrderIterations();
  bool diffSign(int idx);
  */  
public:
  explicit Reinit(const Array2D<double> &_tpx, const Array2D<double> &_tpy, const int _width, const int _flag);
  void reinit();
  void dump_cp(double *dcpx, double *dcpy) const;

};

#endif
