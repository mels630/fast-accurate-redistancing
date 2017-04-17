# file: pyRedist.pyx

from cpython.mem cimport PyMem_Malloc

cimport redistMain
import numpy as np

cdef class pyRedist:
     cdef redistMain.RedistDriver c_rd
     def __cinit__(self):
         pass

     def Run(self, d, n, width, flag):
         self.c_rd.Run(d,n,width,flag)

     def RunRes(self, d, n, width, flag):
         res = self.c_rd.Run(d,n,width,flag,True)
         return np.reshape(np.asarray(res, dtype=float), (n,n) if d == 2 else (n,n,n))
