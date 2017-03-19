% script to make mex-file for redistC
mex -c heap.cpp 
mex -c idarray2d.cpp 
mex -c idarray3d.cpp
mex -c toolbox.cpp
mex -c toolbox3d.cpp
mex -c redist.cpp
mex -c redist_unsigned.cpp
mex -c redist3.cpp
mex -c redistMex.cpp
mex -o redistMex redistMex.o redist.o redist3.o toolbox.o idarray2d.o idarray3d.o heap.o toolbox3d.o redist_unsigned.o -lm -cxx
