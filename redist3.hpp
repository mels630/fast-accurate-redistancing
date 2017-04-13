#ifndef _REDIST3C_H_
#define _REDIST3C_H_

#include "array3d.hpp"
#include "idarray3d.hpp"
#include "heap.hpp"

#include <utility>
#include <array>

using Point3 = std::array<double, 3>;

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
};

using Helt3 = std::pair<double, Aux3>;

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

  void fastMarchingRedist();
  void directionalOptimization();
  void updateAndAddNeighborsToHeap(idx_t const idx);
  void updateAndAddNeighborsToHeapDO(idx_t const idx);
  double estimateUpdate(idx_t const idx);
  void setInterfaceValues();
  void setInterfaceValuesDO();
  void thresholdAwayFromInterface();
  void applyResult(Helt3 const &h);
  struct helt performDO(idx_t const idx);
  struct helt performDO(idx_t const idx, Point3 const &cpguess)
  void lineSearch(idx_t const idx, double (&grad)[3], double (&cpguess)[3]);
  bool findOppSign(idx_t const idx, double (&guess)[3]);
  void bisect(idx_t const idx, double (&guess)[3]);
  double bisect(double (&result)[3], double (&xm)[4], double (&xp)[4]);
  explicit Redist3(); // no empty constructor
  bool bracket(idx_t const idx, double (&cpguess)[3], double (&xm)[4], double (&xp)[4]);
  double findNborDirections(const double (&x0)[3], double (&xx)[5][3], const double delta, double (&rvec)[3], double (&v1)[3], double (&v2)[3]);
  double search1D(idx_t const idx, double (&x)[3]);
  void secondOrderIterations();
  bool diffSign(int idx);
  
public:
  explicit Redist3(const Array3D<double> &_u, idx_t const _width, int const _flag);
  void redistance();
  void dump_u(double *v);
  void dump_cp(double *cpx_d, double *cpy_d, double *cpz_d);
};

#endif
