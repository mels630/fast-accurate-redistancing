#ifndef _TOOLBOX_HPP_
#define _TOOLBOX_HPP_

#include "array2d.hpp"

// Function definitions 
Array2D<double> makeCircle(idx_t const n, double const r, double const xc, double const yc);
Array2D<double> makeCircle(idx_t const n, double const r);
Array2D<double> makeBox(idx_t const n);
double l2err(Array2D<double> const &u, Array2D<double> const &v);
double normalize(std::array<double,2> &vec);
double dist(std::array<double,2> const &x1, std::array<double,2> const &x2);
double dist(std::array<double,2> const &x1, std::array<double,2> const &x2, double const lenx, double const leny);
std::array<double, 2> ccomb(std::array<double, 2> const &x1, std::array<double,2> const &x2, double const theta, double const xlen, double const ylen);
#endif
