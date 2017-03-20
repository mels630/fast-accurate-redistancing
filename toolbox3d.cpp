#include <algorithm>
#include <functional>
#include "defs.h"

#include "toolbox3d.hpp"

double l2err(Array3D<double> const &u, Array3D<double> const &v)
{
  if((u.getn() != v.getn()) || (u.getm() != v.getm()) || (u.getk() != v.getk()))
    return(-1.0f);
  double const err = std::inner_product(u.returnData().begin(), u.returnData().end(), v.returnData().begin(), 0.,
                                        std::plus<double>(), [](double const &d1, double const &d2)->double{return (d1-d2)*(d1-d2); });
  return std::sqrt(err / static_cast<double>(u.getN()));
}

void ccomb(const double (&x1)[3], const double (&x2)[3], double (&xr)[3], double theta)
{ // forms convex combination xr = theta*x1 + (1-theta)*x2
  double diff[3];
  if(theta < 0.0f || theta > 1.0f || std::isnan(theta))
    abort();
  if(std::isnan(theta))
    theta = 0.5f;

  diff[0] = pd(x2[0],x1[0]);
  diff[1] = pd(x2[1],x1[1]);
  diff[2] = pd(x2[2],x1[2]);
  xr[0] = x1[0] + (1.0f-theta)*diff[0];
  xr[1] = x1[1] + (1.0f-theta)*diff[1];
  xr[2] = x1[2] + (1.0f-theta)*diff[2];
}

void ccomb(const double (&x1)[4], const double (&x2)[4], double (&xr)[3], double theta)
{ // forms convex combination xr = theta*x1 + (1-theta)*x2. Ignores third value in x1/x2
  double diff[3];
  if(theta < 0.0f || theta > 1.0f || std::isnan(theta))
    abort();
  if(std::isnan(theta))
    theta = 0.5f;

  diff[0] = pd(x2[0],x1[0]);
  diff[1] = pd(x2[1],x1[1]);
  diff[2] = pd(x2[2],x1[2]);
  xr[0] = x1[0] + (1.0f-theta)*diff[0];
  xr[1] = x1[1] + (1.0f-theta)*diff[1];
  xr[2] = x1[2] + (1.0f-theta)*diff[2];
  //printf(" theta = %.3e, xr = (%f,%f).\n",theta,xr[0],xr[1]); fflush(stdout);
}

void orthoVecs(const double (&in)[3], double (&out1)[3], double (&out2)[3])
{ /* generates vectors out1 and out2 which form an ONB for R^3 along with 
     the vector in. Memory is assumed to be preallocated. in is assumed to be unit length;
     resulting in out1 / out2 which are also unit length */
  generateRandomVec(out1,in);
  mycross(in,out1,out2); /* out2 should be O.N. to in */
  normalize(out2);
  mycross(in,out2,out1); /* out1 should be O.N to in, out2 */
}

void generateRandomVec(double (&out)[3], const double (&in)[3])
{ /* Assume "in" is 3x1 vector with unit length. Generates a unit length vector
     "out" with |dot(in,out)| < 0.5 (not pointing in the same direction */
  int ii, count;
  const int VERBOSITY = 0;
  const int MAXCOUNT = 20;
  double dot;
  dot = 1.0;
  count = 0;
  while((fabs(dot) >= 0.5) && (count < MAXCOUNT))
  {
    for(ii=0;ii<3;ii++)
      out[ii] = ((2.0*(double)rand()) / (double)RAND_MAX) - 1.0; /* U(-1,1) */
    normalize(out);
    dot = 0.0;
    for(ii=0;ii<3;ii++)
      dot += in[ii]*out[ii];
    count++;
  }
  if((count == MAXCOUNT) && VERBOSITY)
    printf("  count equalled MAXCOUNT!\n");
}

void mycross(const double (&A)[3], const double (&B)[3], double (&C)[3])
{ /* returns C = A x B */
  C[0] = A[1]*B[2]-A[2]*B[1];
  C[1] = A[2]*B[0]-A[0]*B[2];
  C[2] = A[0]*B[1]-A[1]*B[0];
}

Array3D<double> makeSphere(const int n, const double r, const double xc, const double yc, const double zc)
{
  Array3D<double> u(n);
  for(int kk=0; kk<n; ++kk)
  {
    double z = static_cast<double>(kk) / static_cast<double>(n);
    for(int ii=0; ii<n; ++ii)
    {
      double x = static_cast<double>(ii) / static_cast<double>(n);
      for(int jj=0; jj<n; ++jj)
      {
        double y = static_cast<double>(jj) / static_cast<double>(n);
        u.put(r-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)),jj,ii,kk);
      }
    }
  }
  return(u);
}

Array3D<double> makeSphere(const int n, const double r)
{
  const double ctr = static_cast<double>(0.5 - 1.0f/static_cast<double>(n));
  return(makeSphere(n,r,ctr,ctr,ctr));
}

double normalize(double (&vec)[3])
{ // assume vec = <vec[0],vec[1],vec[2]>
  double normsq = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
  const double nrm = sqrt(normsq);
  if(normsq < 1e-14)
    if(fabs(vec[0]) > fabs(vec[1]))
    {
      vec[1] = 0.0f;
      if(fabs(vec[0]) > fabs(vec[2]))
      {
        vec[0] = mysign(vec[0]);
        vec[2] = 0.0f;
      }
      else
      {
        vec[2] = mysign(vec[2]);
        vec[0] = 0.0f;
      }
    }
    else
    {
      vec[0] = 0.0f;
      if(fabs(vec[1]) > fabs(vec[2]))
      {
        vec[1] = mysign(vec[1]);
        vec[2] = 0.0f;
      }
      else
      {
        vec[2] = mysign(vec[2]);
        vec[1] = 0.0f;
      }
    }
  else
  {
    vec[0] /= nrm;
    vec[1] /= nrm;
    vec[2] /= nrm;
  }
  return(nrm);
}

int minabs(const double *val, const int amt)
{ /* returns the index of the entry in val with smallest absolute value. Val
     assumed to contain amt entries valued < 1e-64 */
  int ind, ii;
  double minval;

  /* initialize */
  minval = 1e64;
  ind = -1;
  /* loop over entries of val */
  for(ii=0;ii<amt;ii++)
    if(fabs(val[ii]) < minval)
    {
      minval = fabs(val[ii]);
      ind = ii;
    }
  return(ind);
}

double dist3(const double *x1, const double *x2)
{ return(sqrt(pd(x1[0],x2[0])*pd(x1[0],x2[0])+pd(x1[1],x2[1])*pd(x1[1],x2[1])+pd(x1[2],x2[2])*pd(x1[2],x2[2]))); }

