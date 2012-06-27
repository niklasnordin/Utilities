#!/usr/bin/make


EXE = solveDE

CPP =	g++
LD = 	ld
CPPFLAGS = -fPIC
SRC =	solve.C
OBJS =	solve.o

all: $(LIBNAME)

$(EXE): $(OBJS)
	$(CPP) -shared -o $(LIBNAME) $(OBJS)

.C.o:
	$(CPP) -c $(CPPFLAGS) $<

clean:
	rm -f $(OBJS) $(LIBNAME)

solve.C:
