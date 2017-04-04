#ifndef _IDARRAY2D_HPP_
#define _IDARRAY2D_HPP_

#include "array2d.hpp"
#include <vector>

/// IDArray2D: class for interpolable Array2D<double>
class IDArray2D : public Array2D<double> {
private:
  ///\todo Convert qe to vector<unique_ptr<double>> // note, can we assign arrays this way?
  mutable std::vector<double*> qe; ///< The entries of the vector are pointers to arrays containing
                                   ///< the interpolation coeffections for that index of the Array2D<double> data
                                   ///< These are only filled in once interpolation is asked for around the index.
  int const flag;                  ///< Interpolation type: (1=bilinear, 2=biquadratic, 3=bicubic)
  void SetQe(idx_t const idx) const;
public:
  IDArray2D() = delete; // Empty constructor not permitted
  IDArray2D(int const nn, int const _flag);                 // Constructor for nn x nn array
  IDArray2D(int const mm, int const nn, int const _flag);   // Constructor for mm x nn array
  IDArray2D(int const mm, int const nn, double const _dx, double const _dy, int const _flag);   // Constructor for mm x nn array
  IDArray2D(const Array2D<double> &input, int const _flag); // Constructor to deep copy "input" and initialize qe
  ~IDArray2D();                                     // Destructor deallocates memory allocated in entries of qe
  double interpolate(double x, double y) const; // returns interpolated value of data(xx,yy), where data is imagined as a discretization of [0,1)^2
  void interpolate(double x, double y, double(&result)[3]) const; // returns interpolated value of data(xx,yy) in result[0] and gradient of data(x,y) in result[1:2]
};

#endif
