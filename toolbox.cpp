#include <algorithm>
#include <numeric>
#include <functional>
#include "defs.h"

#include "toolbox.hpp"

Array2D<double> makeCircle(int const n, double const r, double const xc, double const yc)
{
  Array2D<double> u(n,n);
  for(int ii=0; ii<n; ++ii)
  {
    double x = static_cast<double>(ii) / static_cast<double>(n);
    for(int jj=0; jj<n; ++jj)
    {
      double y = static_cast<double>(jj) / static_cast<double>(n);
      u.put(r-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)),jj,ii);
    }
  }
  return(u);
}

Array2D<double> makeCircle(int const n, double const r)
{
  double const ctr = static_cast<double>(0.5 - 1.0f/static_cast<double>(n));
  return(makeCircle(n,r,ctr,ctr));
}

Array2D<double> makeBox(int const n)
{ // make a function which is 0.5 inside a square and -0.5 outside.
  Array2D<double> u(n,n);
  for(int ii=0; ii<n; ++ii)
    for(int jj=0; jj<n; ++jj)
      if((ii>=n/4) && (ii<=3*n/4) && (jj>=n/4) && (jj<=3*n/4))
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

double normalize(double (&grad)[2])
{ // assume grad = <grad[0],grad[1]>
  double normsq = grad[0]*grad[0]+grad[1]*grad[1];
  if(normsq < 1e-14)
    if(fabs(grad[0]) > fabs(grad[1]))
    {
      grad[0] = mysign(grad[0]);
      grad[1] = 0.0f;
    }
    else
    {
      grad[0] = 0.0f;
      grad[1] = mysign(grad[1]);
    }
  else
  {
    grad[0] /= sqrt(normsq);
    grad[1] /= sqrt(normsq);
  }
  return(sqrt(normsq));
}

void ccomb(double const * const x1, double const * const x2, double* const xr, double theta, double const xlen, double const ylen)
{ // forms convex combination xr = theta*x1 + (1-theta)*x2
  double diff[2];
  if(theta < 0.0f || theta > 1.0f || std::isnan(theta))
    abort();
  if(std::isnan(theta))
    theta = 0.5f;

  diff[0] = pdl(x2[0],x1[0],xlen);
  diff[1] = pdl(x2[1],x1[1],ylen);
  xr[0] = x1[0] + (1.0f-theta)*diff[0];
  xr[1] = x1[1] + (1.0f-theta)*diff[1];
}

double dist(double const * const x1, double const * const x2)
{
  return sqrt(pd(x1[0],x2[0])*pd(x1[0],x2[0])+pd(x1[1],x2[1])*pd(x1[1],x2[1]));
}

double dist(double const * const x1, double const * const x2, double const lenx, double const leny)
{
  return sqrt(pdl(x1[0],x2[0],lenx)*pdl(x1[0],x2[0],lenx)+pdl(x1[1],x2[1],leny)*pdl(x1[1],x2[1],leny));
}

