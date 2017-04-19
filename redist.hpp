#ifndef _REDISTC_H_
#define _REDISTC_H_

#include "array2d.hpp"
#include "idarray2d.hpp"
#include "heap.hpp"
#include "defs.h"

#include <utility>
#include <array>

/// Auxilary data-carrying POD struct for heap operations
struct Aux
{
public:
  Aux(idx_t const _i, Point const &_aux);
  explicit Aux(idx_t const _i);
  Aux() = default;                        ///< Default empty constructor
  Aux(Aux const &a) = default;            ///< Default copy constructor
  Aux& operator=(Aux const &a) = default; ///< Default assignment operator
  idx_t i;                                ///< Index data
  Point aux;                              ///< Closest point location

  /// Comparison operator
  /// \param[in] other : Aux object to compare to
  /// \return            True if *this is lexicographically less than other
  bool operator<(Aux const &other) const { return i < other.i; } 

  /// Comparison operator
  /// \param[in] other : Aux object to compare to
  /// \return            True if *this is lexicographically greater than other
  bool operator>(Aux const &other) const { return other < *this; }
};

using Helt = std::pair<double, Aux>;

class Redist
{ 
private:
  idx_t const width;      ///< Number of pixels to redistance away from the interface
  int const flag;         ///< Flag for interpolation order
  idx_t const m;          ///< Number of rows in each array
  idx_t const n;          ///< Number of columns in each array
  idx_t const N;          ///< Number of elements in each array
  double const dx;        ///< x grid spacing
  double const dy;        ///< y grid spacing
  double const thres;     ///< Maximum distance to compute to
  bool const cpflag;      ///< Flag indicating whether to compute closest point information
  Heap<double, Aux> h;    ///< Heap for maintaining data <compare performance to std::priority_queue>
  Array2D<bool> state;    ///< State of each element (redistancing is complete here)
  std::vector<idx_t> bnd; ///< List of indices on the boundary of the interface
  IDArray2D u0;           ///< Interpolation array containing initial values
  Array2D<double> u;      ///< Updated array with signed distances
  Array2D<double> cpx;    ///< x-value of closest point on interface
  Array2D<double> cpy;    ///< y-value of closest point on interface

  void fastMarchingRedist();
  void directionalOptimization();
  void updateAndAddNeighborsToHeap(idx_t const idx);
  void updateAndAddNeighborsToHeapDO(idx_t const idx);
  double estimateUpdate(idx_t const idx);
  void setInterfaceValues();
  void setInterfaceValuesDO();
  void thresholdAwayFromInterface();
  void applyResult(Helt const &h);
  Helt performDO(idx_t const idx);
  Helt performDOSurf(idx_t const idx);
  Helt performDO(idx_t const idx, Point const &cpguess);
  Point lineSearch(idx_t const idx, Point &grad);
  bool findOppSign(idx_t const idx, Point &guess);
  void bisect(idx_t const idx, Point &guess);
  double bisect(Point &result, std::pair<Point, double> const &xm, std::pair<Point, double> const &xp);
  Redist() = delete; // no empty constructor
  bool bracket(idx_t const idx, Point const &cpguess, std::pair<Point, double> &xm, std::pair<Point, double> &xp);
  void findNborDirections(Point const &x0, Point const &xc, Point &xl, Point &xr, double const delta);
  double search1D(idx_t const idx, Point &x);
  void secondOrderIterations();
  bool diffSign(idx_t const idx);
  
public:
  Redist(Array2D<double> const &_u, idx_t const _width, int const _flag);
  Redist(Array2D<double> const &_u, double const _dx, double const _dy, idx_t const _width, int const _flag);
  void redistance();
  void dump_u(double *v);
  Array2D<double> const& dump_u();
  
  void dump_cp(double *cpxArr, double *cpyArr);
  
};

#endif
