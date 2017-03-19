% script to make mex-file for redistC
mex -c heap.cpp 
mex -c idarray2d.cpp 
mex -c reinit.cpp
mex -c toolbox.cpp
mex -c cpReinitMex.cpp
mex -o cpReinitMex cpReinitMex.o reinit.o idarray2d.o heap.o toolbox.o -lm
