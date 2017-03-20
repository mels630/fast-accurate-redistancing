#ifndef _TOOLBOX3D_HPP_
#define _TOOLBOX3D_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <climits>
#include "sys/stat.h"
#include "sys/resource.h"
#include "dirent.h"
#include "defs.h"
#include "array3d.hpp"
#include "assert.h"

using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::vector;

void fillSign(Array3D<double> &u, const double nullval, const double newval);
void fillSignSubfn(Array3D<double> &u, const double nullval, vector<int> &list, const int curr, const double newval, const double TOL);
double l2err(const Array3D<double> u, const Array3D<double> v);
void ccomb(const double (&x1)[3], const double (&x2)[3], double (&xr)[3], double theta);
void ccomb(const double (&x1)[4], const double (&x2)[4], double (&xr)[3], double theta);
inline double dist3(const double *x1, const double *x2);
void orthoVecs(const double (&in)[3], double (&out1)[3], double (&out2)[3]);
void generateRandomVec(double (&out)[3], const double (&in)[3]);
void mycross(const double (&A)[3], const double (&B)[3], double (&C)[3]);
Array3D<double> makeSphere(const int n, const double r, const double xc, const double yc, const double zc);
Array3D<double> makeSphere(const int n, const double r);
double normalize(double (&vec)[3]);
double normalize(struct point3 &vec);
inline double dotProd(const struct point3 &a, const struct point3 &b);
int minabs(const double *val, const int amt);
inline struct point3 crossProd(const struct point3 &a, const struct point3 &b);

inline double dist3(const double *x1, const double *x2)
{ return(sqrt(pd(x1[0],x2[0])*pd(x1[0],x2[0])+pd(x1[1],x2[1])*pd(x1[1],x2[1])+pd(x1[2],x2[2])*pd(x1[2],x2[2]))); }

inline double dotProd(const struct point3 &a, const struct point3 &b)
{
  return(a.x*b.x+a.y*b.y+a.z*b.z);
}

inline struct point3 crossProd(const struct point3 &a, const struct point3 &b)
{
  struct point3 lval;
  lval.x = a.y*b.z-a.z*b.y;
  lval.y = a.z*b.x-a.x*b.z;
  lval.z = a.x*b.y-a.y*b.x;
  return(lval);
}

#endif
