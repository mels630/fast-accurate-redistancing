#ifndef _TOOLBOX_HPP_
#define _TOOLBOX_HPP_

#include "array2d.hpp"
#include "defs.h"

// Function definitions 
Array2D<double> makeCircle(idx_t const n, double const r, double const xc, double const yc);
Array2D<double> makeCircle(idx_t const n, double const r);
Array2D<double> makeBox(idx_t const n);
double l2err(Array2D<double> const &u, Array2D<double> const &v);
double normalize(Point &vec);
double dist(Point const &x1, Point const &x2);
double dist(Point const &x1, Point const &x2, double const xlen, double const ylen);
Point ccomb(Point const &x1, Point const &x2, double const theta, double const xlen, double const ylen);
#endif
