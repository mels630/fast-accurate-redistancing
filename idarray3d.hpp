#ifndef _IDARRAY3D_HPP_
#define _IDARRAY3D_HPP_

#include <vector>

#include "array3d.hpp"

class IDArray3D : public Array3D<double> {
private:
  explicit IDArray3D();                             // Empty constructor not permitted

  std::vector<double*> qe; // The entries of the vector are pointers to arrays containing
                      // the interpolation coeffections for that index of the Array3D<double> data
                      // These are only filled in once interpolation is asked for around the index.
  const int flag;     // Interpolation type: (2=triquadratic) <-- currently only option
public:
  explicit IDArray3D(const int nn, const int _flag);                 // Constructor for nn x nn array
  explicit IDArray3D(const int mm, const int nn, const int kk, const int _flag);   // Constructor for mm x nn x kk array
  explicit IDArray3D(const Array3D<double> &input, const int _flag); // Constructor to deep copy "input" and initialize qe
  ~IDArray3D();                                     // Destructor deallocates memory allocated in entries of qe
  double interpolate(double x, double y, double z); // returns interpolated value of data(x,y,z), where data is imagined as a discretization of [0,1)^3
};

#endif
