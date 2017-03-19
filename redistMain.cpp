#include <iostream>
#include <fstream>
#include <string>
#include "redist.hpp"
#include "redist3.hpp"
#include "defs.h"
#include "array2d.hpp"
#include "array3d.hpp"
#include "toolbox.hpp"
#include <cassert>
using std::cout;
using std::cin;

int main(int argc, char **argv)
{
  // Need to obtain a Array2D<double> [u] for which the zero-level set should
  // be redistanced to, out to [width] pixels (thinking of domain as [0,1)^2)
  // with flag = 1 : fast marching method (Tsitsiklis's algorithm)
  //             2 : directional optimization with biquadratic interpolation
  //             3 : directional optimization with bicubic interpolation
  // (abort on errors)

  // Assume square initially ... (?)
  int d, n, width, flag;
  if(argc == 1) {
    cout << "Enter dimension: ";
    cin >> d;
    cout << "Enter n: ";
    cin >> n;
    cout << "Enter width: ";
    cin >> width;
    cout << "Enter flag: ";
    cin >> flag;
  } else {
    assert (argc > 4);
    d = std::atoi(argv[1]);
    n = std::atoi(argv[2]);
    width = std::atoi(argv[3]);
    flag = std::atoi(argv[4]);
  }
  
  std::cout << "Settings:" << std::endl;
  std::cout << "        d: " << d << std::endl;
  std::cout << "        n: " << n << std::endl;
  std::cout << "    width: " << width << std::endl;
  std::cout << "     flag: " << flag << std::endl;
      
  if(d == 2)
  {
    Array2D<double> u0 = makeCircle(n,.25,.5,.5);
    //Array2D<double> v = makeCircle(n,.35,.5,.5);
    //Array2D<double> u0 = makeBox(n);

    // Perform the redistancing
    Redist r(u0,width,flag);
    r.redistance();
    cout << endl << "Calculated L2-error is " << l2err(r.dump_u(),u0) << endl;  
  }
  else if(d == 3)
  {
    Array3D<double> u = makeSphere(n,.25,.5,.5,.5);

    // Perform the redistancing
    Redist3 r(u,width,flag);
    r.redistance();
    cout << endl << "Calculated L2-error is " << l2err(r.u,r.u0) << endl;  
  }
  else
    cout << "d must be 2 or 3." << endl;

  return(0);
}

