from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os

ext = Extension("pyRedist",
                ["pyRedist.pyx"],
                language="c++",
                include_dirs=["/home/melsey/Git/fast-accurate-redistancing/"],
                library_dirs=["/home/melsey/Git/fast-accurate-redistancing/"],
                libraries=["redistMain"],
                extra_compile_args=["-std=c++11"])

setup(
    ext_modules = cythonize([ext])
)
