#include <iostream>
#include <fstream>
#include <string>
#include "redist.hpp"
#include "redist3.hpp"
#include "defs.h"
#include "array2d.hpp"
#include "array3d.hpp"
#include "toolbox.hpp"
#include "toolbox3d.hpp"
#include "redistMain.hpp"
#include <cassert>

void RedistDriver::Run(int const d, int const n, int const width, int const flag)
{
  std::vector<double> vd = Run(d, n, width, flag, false);
  assert(vd.empty());
}

std::vector<double> RedistDriver::Run(int const d, int const n, int const width, int const flag, bool const bReturnArray)
{
  std::cout << "Settings:" << std::endl;
  std::cout << "        d: " << d << std::endl;
  std::cout << "        n: " << n << std::endl;
  std::cout << "    width: " << width << std::endl;
  std::cout << "     flag: " << flag << std::endl;

  std::vector<double> vdu;
  if(d == 2)
  {
    Array2D<double> u0 = makeCircle(n,.25,.5,.5);

    // Perform the redistancing
    Redist r(u0,width,flag);
    r.redistance();
    std::cout << std::endl << "Calculated L2-error is " << l2err(r.dump_u(), u0) << std::endl;
    if (bReturnArray)
      vdu = r.dump_u().returnData();
  }
  else if(d == 3)
  {
    Array3D<double> u0 = makeSphere(n,.25,.5,.5,.5); 

    // Perform the redistancing
    Redist3 r(u0,width,flag);
    r.redistance();
    std::cout << std::endl << "Calculated L2-error is " << l2err(r.dump_u(), u0) << std::endl;

    if (bReturnArray)
      vdu = r.dump_u().returnData();
  }
  else
    std::cout << "d must be 2 or 3." << std::endl;
  return vdu;
}

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
    std::cout << "Enter dimension: ";
    std::cin >> d;
    std::cout << "Enter n: ";
    std::cin >> n;
    std::cout << "Enter width: ";
    std::cin >> width;
    std::cout << "Enter flag: ";
    std::cin >> flag;
  } else {
    assert (argc > 4);
    d = std::atoi(argv[1]);
    n = std::atoi(argv[2]);
    width = std::atoi(argv[3]);
    flag = std::atoi(argv[4]);
  }
  
  RedistDriver rd;
  rd.Run(d,n,width,flag);
  return(0);
}

