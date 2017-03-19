#ifndef _IDARRAY2D_HPP_
#define _IDARRAY2D_HPP_

#include "array2d.hpp"

class IDArray2D : public Array2D<double> {
private:
  explicit IDArray2D();                             // Empty constructor not permitted

  mutable std::vector<double*> qe; // The entries of the vector are pointers to arrays containing
                                   // the interpolation coeffections for that index of the Array2D<double> data
                                   // These are only filled in once interpolation is asked for around the index.
  const int flag;                  // Interpolation type: (1=bilinear, 2=biquadratic, 3=bicubic)
public:
  explicit IDArray2D(const int nn, const int _flag);                 // Constructor for nn x nn array
  explicit IDArray2D(const int mm, const int nn, const int _flag);   // Constructor for mm x nn array
  explicit IDArray2D(const int mm, const int nn, const double _dx, const double _dy, const int _flag);   // Constructor for mm x nn array
  explicit IDArray2D(const Array2D<double> &input, const int _flag); // Constructor to deep copy "input" and initialize qe
  ~IDArray2D();                                     // Destructor deallocates memory allocated in entries of qe
  double interpolate(const double xx, double yy) const; // returns interpolated value of data(xx,yy), where data is imagined as a discretization of [0,1)^2
  void interpolate(const double xx, const double yy, double(&result)[3]) const; // returns interpolated value of data(xx,yy) in result[0] and gradient of data(x,y) in result[1:2]
  void freeQeIdx(const int idx) const; // clear interpolation coefficients at idx
  void freeQeIdxAll() const; // clear all interpolation coefficients
};

#endif
