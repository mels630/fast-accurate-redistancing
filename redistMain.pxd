# file: redistMain.pxd

from libcpp.vector cimport vector
from libcpp cimport bool as bool_t

cdef extern from "redistMain.hpp":
     cdef cppclass RedistDriver:
          RedistDriver() except+
          void Run(int d, int n, int width, int flag);
          vector[double] Run(int d, int n, int width, int flag, bool_t bUseArray);
          