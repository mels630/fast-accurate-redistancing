#ifndef _TOOLBOX_HPP_
#define _TOOLBOX_HPP_

#include "array2d.hpp"

// Function definitions 
Array2D<double> makeCircle(int const n, double const r, double const xc, double const yc);
Array2D<double> makeCircle(int const n, double const r);
Array2D<double> makeBox(int const n);
double l2err(Array2D<double> const &u, Array2D<double> const &v);
double normalize(double (&grad)[2]);
double dist(double const * const x1, double const * const x2);
double dist(double const * const x1, double const * const x2, double const lenx, double const leny);
void ccomb(double const * const x1, double const * const x2, double* const xr, double theta, double const xlen, double const ylen);

#endif
