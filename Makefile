# Parameters to control Makefile operation
ERR = $(shell which icpc >/dev/null; echo $$?)
ifeq "$(ERR)" "0"
    CC = icpc
else
    CC = g++
endif

CFLAGS = -O2 -Wall -std=c++11 -fPIC -D__USE_STL__
#CFLAGS = -O2 -Wall -std=c++11 -fPIC
#CFLAGS= -pg -O2 -DDEBUG

# Entries to bring executable up to date

all: redistMain libredistMain.so

clean:
	rm *.o -f

redistMain: redistMain.o redist.o heap.o toolbox.o array2d.hpp defs.h idarray2d.o redist3.o array3d.hpp idarray3d.o toolbox3d.o 
	$(CC) $(CFLAGS) -o redistMain redistMain.o redist.o redist3.o heap.o toolbox.o idarray2d.o idarray3d.o -lm toolbox3d.o

libredistMain.so: redistMain.o redist.o heap.o toolbox.o array2d.hpp defs.h idarray2d.o redist3.o array3d.hpp idarray3d.o toolbox3d.o
	$(CC) $(CFLAGS) redistMain.o redist.o redist3.o heap.o toolbox.o idarray2d.o idarray3d.o -lm toolbox3d.o -shared -o libredistMain.so

redistMain.o: redistMain.cpp array2d.hpp array3d.hpp defs.h redist.hpp redist3.hpp
	$(CC) $(CFLAGS) -c redistMain.cpp

redist.o: redist.cpp redist.hpp array2d.hpp defs.h heap.hpp idarray2d.hpp toolbox.hpp
	$(CC) $(CFLAGS) -c redist.cpp

redist3.o: redist3.cpp redist3.hpp array3d.hpp defs.h heap.hpp idarray3d.hpp toolbox3d.hpp
	$(CC) $(CFLAGS) -c redist3.cpp

toolbox.o: toolbox.cpp toolbox.hpp array2d.hpp defs.h
	$(CC) $(CFLAGS) -c toolbox.cpp

toolbox3d.o: toolbox3d.cpp toolbox3d.hpp array3d.hpp defs.h
	$(CC) $(CFLAGS) -c toolbox3d.cpp

heap.o: heap.cpp defs.h 
	$(CC) $(CFLAGS) -c heap.cpp

idarray2d.o: idarray2d.cpp idarray2d.hpp array2d.hpp defs.h
	$(CC) $(CFLAGS) -c idarray2d.cpp

idarray3d.o: idarray3d.cpp idarray3d.hpp array3d.hpp defs.h
	$(CC) $(CFLAGS) -c idarray3d.cpp
