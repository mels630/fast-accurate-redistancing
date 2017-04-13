#ifndef _TOOLBOX3D_HPP_
#define _TOOLBOX3D_HPP_

#include "array3d.hpp"
#include "defs.h"

// Function definitions
Array3D<double> makeSphere(idx_t const n, double const r, double const xc, double const yc, double const zc);
Array3D<double> makeSphere(idx_t const n, double const r);
double l2err(Array3D<double> const &u, Array3D<double> const &v);
double normalize(Point3 &vec);
double dist3(Point3 const &x1, Point3 const &x2);
double dist3(Point3 const &x1, Point3 const &x2, double const xlen, double const ylen, double const zlen);
Point3 ccomb(Point3 const &x1, Point3 const &x2, double const theta);
Point3 ccomb(Point3 const &x1, Point3 const &x2, double const theta, double const xlen, double const ylen, double const zlen);
void orthoVecs(Point3 const &in, Point3 &out1, Point3 &out2);
void generateRandomVec(Point3 &out, Point3 const &in);
Point3 mycross(Point3 const &A, Point3 const &B);
//int minabs(double const *val, idx_t const amt);

#endif
