#ifndef __REDISTMAIN_HPP__
#define __REDISTMAIN_HPP__

#include <vector>

class RedistDriver
{
public:
  RedistDriver() = default;
  void Run(int const d, int const n, int const width, int const flag);
  std::vector<double> Run(int const d, int const n, int const width, int const flag, bool const bReturnArray);
};

#endif // __REDISTMAIN_HPP__
