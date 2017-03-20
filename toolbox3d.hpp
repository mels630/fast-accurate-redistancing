#ifndef _TOOLBOX3D_HPP_
#define _TOOLBOX3D_HPP_

#include "array3d.hpp"

double l2err(Array3D<double> const &u, Array3D<double> const &v);
void ccomb(const double (&x1)[3], const double (&x2)[3], double (&xr)[3], double theta);
void ccomb(const double (&x1)[4], const double (&x2)[4], double (&xr)[3], double theta);
double dist3(const double *x1, const double *x2);
void orthoVecs(const double (&in)[3], double (&out1)[3], double (&out2)[3]);
void generateRandomVec(double (&out)[3], const double (&in)[3]);
void mycross(const double (&A)[3], const double (&B)[3], double (&C)[3]);
Array3D<double> makeSphere(const int n, const double r, const double xc, const double yc, const double zc);
Array3D<double> makeSphere(const int n, const double r);
double normalize(double (&vec)[3]);
int minabs(const double *val, const int amt);

#endif
