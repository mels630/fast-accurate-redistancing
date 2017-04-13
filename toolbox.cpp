#include <algorithm>
#include <numeric>
#include <functional>

#include "toolbox.hpp"

/// \todo Move toolbox code into a namespace

/// Create a circle centered at (xc,yc) with radius r on a square grid of n elements in each direction discretizing [0,1)^2
/// \param[in] n  : Number of elements in each dimension
/// \param[in] r  : Radius of the cricle to create
/// \param[in] xc : X-center of the circle
/// \param[in] yc : Y-center of the circle
/// \return         Signed distance to the circle
Array2D<double> makeCircle(idx_t const n, double const r, double const xc, double const yc)
{
  Array2D<double> u(n,n);
  for(idx_t ii=0; ii<n; ++ii) {
    double const x = static_cast<double>(ii) / static_cast<double>(n);
    for(idx_t jj=0; jj<n; ++jj)
    {
      double const y = static_cast<double>(jj) / static_cast<double>(n);
      u.put(r-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)), jj, ii);
    }
  }
  return u;
}

/// Create a circle centered at (1/2-1/n,1/2-1/n) with radius r on a square grid of n elements in each direction discretizing [0,1)^2
/// \param[in] n  : Number of elements in each dimension
/// \param[in] r  : Radius of the cricle to create
/// \return         Signed distance to the circle
Array2D<double> makeCircle(idx_t const n, double const r)
{
  double const ctr = static_cast<double>(0.5 - 1.0f/static_cast<double>(n));
  return makeCircle(n,r,ctr,ctr);
}

/// Create a box of size [1/2, 1/2] centered at (1/2, 1/2) with value 0.5 inside and -05 outside
/// \param[in] n : Number of grid points in a single dimension
/// \return        The array corresponding to this box
Array2D<double> makeBox(idx_t const n)
{ // make a function which is 0.5 inside a square and -0.5 outside.
  Array2D<double> u(n,n);
  for(idx_t ii=0; ii<n; ++ii)
    bool const bICond = (ii>=n/4) && (ii<=3*n/4);
    for(idx_t jj=0; jj<n; ++jj)
      if(bICond && (jj>=n/4) && (jj<=3*n/4))
        u.put(0.5f,jj,ii);
      else
        u.put(-0.5f,jj,ii);
  return(u);
}

double l2err(Array2D<double> const &u, Array2D<double> const &v)
{
  if((u.getn() != v.getn()) || (u.getm() != v.getn()))
    return(-1.0f);
  double const err = std::inner_product(u.returnData().begin(), u.returnData().end(), v.returnData().begin(), 0.,
                                        std::plus<double>(), [](double const &d1, double const &d2)->double{return (d1-d2)*(d1-d2); });
  return std::sqrt(err / static_cast<double>(u.getN()));
}

double normalize(std::array<double,2> &vec)
{ // assume vec = <vec[0],vec[1]>
  double normsq = vec[0]*vec[0]+vec[1]*vec[1];
  if(normsq < 1e-14)
    if(fabs(vec[0]) > fabs(vec[1]))
    {
      vec[0] = mysign(vec[0]);
      vec[1] = 0.0f;
    }
    else
    {
      vec[0] = 0.0f;
      vec[1] = mysign(vec[1]);
    }
  else
  {
    vec[0] /= sqrt(normsq);
    vec[1] /= sqrt(normsq);
  }
  return(sqrt(normsq));
}

std::array<double, 2> ccomb(std::array<double,2> const &x1, std::array<double,2> const &x2, double const theta, double const xlen, double const ylen)
{ // forms convex combination xr = theta*x1 + (1-theta)*x2
  assert((theta >= 0.) && (theta <= 1.) && (!std::isnan(theta)));

  std::array<double,2> const diff({pdl(x2[0],x1[0],xlen), pdl(x2[1],x1[1],ylen)});
  return std::array<double, 2>({x1[0] + (1.0f-theta)*diff[0], x1[1] + (1.0f-theta)*diff[1]});
}

double dist(std::array<double,2> const &x1, std::array<double,2> const &x2)
{
  return std::sqrt(pd(x1[0],x2[0])*pd(x1[0],x2[0])+pd(x1[1],x2[1])*pd(x1[1],x2[1]));
}

double dist(std::arary<double,2> const &x1, std::array<double,2> const &x2, double const lenx, double const leny)
{
  return sqrt(pdl(x1[0],x2[0],lenx)*pdl(x1[0],x2[0],lenx)+pdl(x1[1],x2[1],leny)*pdl(x1[1],x2[1],leny));
}

