#include <algorithm>
#include <numeric>
#include <functional>

#include "defs.h"
#include "toolbox3d.hpp"

/// \todo Move toolbox3d code into a namespace

/// Create a sphere centered at (xc,yc,zc) with radius r on a square grid of n elements in each direction discretizing [0,1)^3
/// \param[in] n  : Number of elements in each dimension
/// \param[in] r  : Radius of the sphere to create
/// \param[in] xc : X-center of the sphere
/// \param[in] yc : Y-center of the sphere
/// \param[in] zc : Z-center of the sphere
/// \return         Signed distance to the sphere
Array3D<double> makeSphere(idx_t const n, double const r, double const xc, double const yc, double const zc)
{
  Array3D<double> u(n);
  for(idx_t kk=0; kk<n; ++kk) {
    double const z = static_cast<double>(kk) / static_cast<double>(n);
    for(idx_t ii=0; ii<n; ++ii) {
      double const x = static_cast<double>(ii) / static_cast<double>(n);
      for(idx_t jj=0; jj<n; ++jj) {
          double y = static_cast<double>(jj) / static_cast<double>(n);
          u.put(r-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)),jj, ii, kk);
      }
    }
  }
  return(u);
}

/// Create a sphere centered at (1/2-1/n,1/2-1/n,1/2-1/n) with radius r on a square grid of n elements in each direction discretizing [0,1)^3
/// \param[in] n  : Number of elements in each dimension
/// \param[in] r  : Radius of the sphere to create
/// \return         Signed distance to the sphere
Array3D<double> makeSphere(idx_t const n, double const r)
{
  double const ctr = static_cast<double>(0.5 - 1.0f/static_cast<double>(n));
  return(makeSphere(n,r,ctr,ctr,ctr));
}

/// Compute the L2 difference between arrays (treated as vectors)
/// \param[in] u : First array to compare
/// \param[in] v : Second array to compare
/// \return        L2-difference, or -1. if shapes don't match
double l2err(Array3D<double> const &u, Array3D<double> const &v)
{
  if((u.getn() != v.getn()) || (u.getm() != v.getm()) || (u.getk() != v.getk()))
    return(-1.0f);
  double const err = std::inner_product(u.returnData().begin(), u.returnData().end(), v.returnData().begin(), 0.,
                                        std::plus<double>(), [](double const &d1, double const &d2)->double{return (d1-d2)*(d1-d2); });
  return std::sqrt(err / static_cast<double>(u.getN()));
}

/// Normalize the vector
/// \param[in,out] vec : Vector to normalize, in-place
/// \return              Length of the vector prior to normalization
double normalize(Point3 &vec)
{ // assume vec = <vec[0],vec[1],vec[2]>
  double const normsq = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
  double const nrm = std::sqrt(normsq);
  if(normsq < 1e-14)
    if(std::abs(vec[0]) > std::abs(vec[1])) {
      vec[1] = 0.0f;
      if(std::abs(vec[0]) > std::abs(vec[2])) {
        vec[0] = mysign(vec[0]);
        vec[2] = 0.0f;
      } else {
        vec[2] = mysign(vec[2]);
        vec[0] = 0.0f;
      }
    } else {
      vec[0] = 0.0f;
      if(std::abs(vec[1]) > std::abs(vec[2])) {
        vec[1] = mysign(vec[1]);
        vec[2] = 0.0f;
      } else {
        vec[2] = mysign(vec[2]);
        vec[1] = 0.0f;
      }
    }
  else {
    vec[0] /= nrm;
    vec[1] /= nrm;
    vec[2] /= nrm;
  }
  return nrm;
}

/// Return a convex combination of two points (theta x1 + (1-theta) x2)
/// \param[in] x1    : First point
/// \param[in] x2    : Second point
/// \param[in] theta : Combination parameter
Point3 ccomb(Point3 const &x1, Point3 const &x2, double const theta)
{ // forms convex combination xr = theta*x1 + (1-theta)*x2
  return ccomb(x1, x2, theta, 1., 1., 1.);
}

/// Return a convex combination of two points (theta x1 + (1-theta) x2)
/// \param[in] x1    : First point
/// \param[in] x2    : Second point
/// \param[in] theta : Combination parameter
/// \param[in] xlen  : Length of space in x (for periodicity)
/// \param[in] ylen  : Length of space in y (for periodicity)
/// \param[in] zlen  : Length of space in z (for periodicity)
Point3 ccomb(Point3 const &x1, Point3 const &x2, double const theta, double const xlen, double const ylen, double const zlen)
{ // forms convex combination xr = theta*x1 + (1-theta)*x2
  assert((theta >= 0.) && (theta <= 1.) && (!std::isnan(theta)));
  Point3 const diff({pdl(x2[0],x1[0],xlen), pdl(x2[1],x1[1],ylen), pdl(x2[2],x1[2],zlen)});
  return Point3({x1[0] + (1.-theta)*diff[0],
        x1[1] + (1.-theta)*diff[1],
        x1[2] + (1.-theta)*diff[2]});
}

/// Compute the distance between two points
/// \param[in] x1 : First point
/// \param[in] x2 : Second point
/// \return         L2 distance between x1 and x2
double dist(Point3 const &x1, Point3 const &x2)
{
  Point3 const diff({pd(x1[0],x2[0]), pd(x1[1],x2[1]), pd(x1[2],x2[2])});
  return std::sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
}

/// Compute the periodic distance between two points
/// \param[in] x1   : First point
/// \param[in] x2   : Second point
/// \param[in] xlen : Length of space in x
/// \param[in] ylen : Length of space in y
/// \param[in] zlen : Length of space in z
/// \return         L2 distance between x1 and x2
double dist(Point3 const &x1, Point3 const &x2, double const xlen, double const ylen, double const zlen)
{
  Point3 const diff({pdl(x2[0],x1[0],xlen), pdl(x2[1],x1[1],ylen), pdl(x2[2],x1[2],zlen)});
  return std::sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
}

/// Generates vectors out1 and out2 which form an ONB for R^3 along with the vector in
/// \param[in]  in   : Input vector (assumed to be unit length
/// \param[out] out1 : Output vector ON to in and out2
/// \param[out] out2 : Output vector ON to in and out1
void orthoVecs(Point3 const &in, Point3 &out1, Point3 &out2)
{
  out1 = generateRandomVec(in);
  out2 = mycross(in,out1); /* out2 should be O.N. to in */
  normalize(out2);
  out1 = mycross(in,out2); /* out1 should be O.N to in, out2 */
}

/// Create a random vector which does not point in nearly the same direction as in
/// \param[in] in : Input vector, assumed unit length
/// \return         Vector in reasonably different direction than in
Point3 generateRandomVec(Point3 const &in)
{ /* Assume "in" is 3x1 vector with unit length. Generates a unit length vector
     "out" with |dot(in,out)| < 0.5 (not pointing in the same direction */
  unsigned int const VERBOSITY = 0;
  unsigned int const MAXCOUNT = 200;
  unsigned int count = 0;
  Point3 out;
  do {
    for (double &d : out)
      d = ((2.0*(double)rand()) / (double)RAND_MAX) - 1.0; /* U(-1,1) */
    normalize(out);
  } while((std::abs(std::inner_product(in.begin(), in.end(), out.begin(), 0.)) >= 0.5) && (count++ < MAXCOUNT));
    
  if((count == MAXCOUNT) && VERBOSITY)
    printf("  count equalled MAXCOUNT!\n");
  
  return out;
}

/// Returns the cross product A x B
/// \param[in] A : first vector
/// \param[in] B : second vector
/// \return        A x B
Point3 mycross(Point3 const &A, Point3 const &B)
{ /* returns C = A x B */
  return Point3({A[1]*B[2]-A[2]*B[1],
        A[2]*B[0]-A[0]*B[2],
        A[0]*B[1]-A[1]*B[0]});
}

/*
int minabs(double const *val, idx_t const amt)
{ // returns the index of the entry in val with smallest absolute value. Val
  //  assumed to contain amt entries valued < 1e-64 
  int ind, ii;
  double minval;

  // initialize 
  minval = 1e64;
  ind = -1;
  // loop over entries of val 
  for(ii=0;ii<amt;ii++)
    if(std::abs(val[ii]) < minval)
    {
      minval = std::abs(val[ii]);
      ind = ii;
    }
  return(ind);
}
*/
