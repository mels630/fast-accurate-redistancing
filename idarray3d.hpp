#ifndef _IDARRAY3D_HPP_
#define _IDARRAY3D_HPP_

#include "array3d.hpp"
#include <vector>

/// IDArray3D: class for interpolable Array3D<double>
class IDArray3D : public Array3D<double> {
private:
  ///\todo Convert qe to vector<unique_ptr<double>> // note, can we assign arrays this way?

  mutable std::vector<double*> qe; ///< The entries of the vector are pointers to arrays containing
                                   ///< the interpolation coeffections for that index of the Array3D<double> data
                                   ///< These are only filled in once interpolation is asked for around the index.
  int const flag;          ///< Interpolation type: (2=triquadratic) <-- currently only option
  void SetQe(idx_t const idx) const;
public:
  IDArray3D() = delete; // Empty constructor not permitted
  IDArray3D(idx_t const nn, int const _flag);                 // Constructor for nn x nn array
  IDArray3D(idx_t const mm, idx_t const nn, idx_t const kk, int const _flag);   // Constructor for mm x nn x kk array
  IDArray3D(Array3D<double> const &input, int const _flag); // Constructor to deep copy "input" and initialize qe
  ~IDArray3D();                                     // Destructor deallocates memory allocated in entries of qe
  double interpolate(double x, double y, double z) const; // returns interpolated value of data(x,y,z), where data is imagined as a discretization of [0,1)^3
};

#endif
