#ifndef _REDIST3C_H_
#define _REDIST3C_H_

#include "array3d.hpp"
#include "idarray3d.hpp"
#include "heap.hpp"
#include "defs.h"

#include <utility>
#include <array>

/// Auxilary data-carrying POD struct for heap operations
struct Aux3
{
public:
  Aux3(idx_t const _i, Point3 const &_aux);
  explicit Aux3(idx_t const _i);
  Aux3() = default;                         ///< Default empty constructor
  Aux3(Aux3 const &a) = default;            ///< Default copy constructor
  Aux3& operator=(Aux3 const &a) = default; ///< Default assignment operator
  idx_t i;                                  ///< Index data
  Point3 aux;                               ///< Closest point location
  bool operator<(Aux3 const &other) const { return i < other.i; }
  bool operator>(Aux3 const &other) const { return other < *this; }
};

using Helt3 = std::pair<double, Aux3>;
using PrPtD = std::pair<Point3, double>;

class Redist3
{ 
private:
  idx_t const width;      ///< Number of pixels to redistance away from the interface
  int const flag;         ///< Flag for interpolation order
  idx_t const m;          ///< Number of rows in each array
  idx_t const n;          ///< Number of columns in each array
  idx_t const k;          ///< Number of slices in each array
  idx_t const N;          ///< Number of elements in each array
  double const dx;        ///< x grid spacing
  double const dy;        ///< y grid spacing
  double const dz;        ///< z grid spacing
  double const thres;     ///< Maximum distance to compute to
  bool const cpflag;      ///< Flag indicating whether to compute closest point information
  Heap<double, Aux3> h;   ///< Heap for maintaining data <compare performance to std::priority_queue>
  Array3D<bool> state;    ///< State of each element (redistancing is complete here)
  std::vector<idx_t> bnd; ///< List of indices on the boundary of the interface
  IDArray3D u0;           ///< Interpolation array containing initial values
  Array3D<double> u;      ///< Updated array with signed distances
  Array3D<double> cpx;    ///< x-value of closest point on interface
  Array3D<double> cpy;    ///< y-value of closest point on interface
  Array3D<double> cpz;    ///< z-value of closest point on interface
  bool bWarn;             ///< Flag for warning user about non-existence of secondOrderIterations
  
  void fastMarchingRedist();
  void directionalOptimization();
  void updateAndAddNeighborsToHeap(idx_t const idx);
  void updateAndAddNeighborsToHeapDO(idx_t const idx);
  double estimateUpdate(idx_t const idx);
  void setInterfaceValues();
  void setInterfaceValuesDO();
  void thresholdAwayFromInterface();
  void applyResult(Helt3 const &h);
  Helt3 performDO(idx_t const idx);
  Helt3 performDOSurf(idx_t const idx);
  Helt3 performDO(idx_t const idx, Point3 const &cpguess);
  Point3 lineSearch(idx_t const idx, Point3 &grad);
  bool findOppSign(idx_t const idx, Point3 &guess);
  void bisect(idx_t const idx, Point3 &guess);
  double bisect(Point3 &result, std::pair<Point3, double> const &xm, std::pair<Point3, double> const &xp);
  Redist3() = delete; // no empty constructor
  bool bracket(idx_t const idx, Point3 const &cpguess, std::pair<Point3, double> &xm, std::pair<Point3, double> &xp);
  double findNborDirections(Point3 const &x0, std::array<PrPtD, 5> &xx, double const delta, Point3 &rvec, Point3 &v1, Point3 &v2);
  double search1D(idx_t const idx, Point3 &x);
  void secondOrderIterations();
  bool diffSign(idx_t const idx);
  
public:
  Redist3(Array3D<double> const &_u, idx_t const _width, int const _flag);
  void redistance();
  void dump_u(double *v);
  Array3D<double> const& dump_u();
  void dump_cp(double *cpx_d, double *cpy_d, double *cpz_d);
};

#endif
