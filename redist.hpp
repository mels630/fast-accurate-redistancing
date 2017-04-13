#ifndef _REDISTC_H_
#define _REDISTC_H_

#include "array2d.hpp"
#include "idarray2d.hpp"
#include "heap.hpp"

#include <utility>
#include <array>

using Point = std::array<double, 2>;

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
  std::array<double, 2> lineSearch(idx_t const idx, Point &grad);
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
  Redist(const Array2D<double> &_u, int const _width, int const _flag);
  Redist(const Array2D<double> &_u, double const _dx, double const _dy, int const _width, int const _flag);
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
