#!/usr/bin/make


EXE = solveDE

LIBS = ../Matrix/libMatrix.a
CPP =	g++
LD = 	ld
CPPFLAGS = -O2
INC = -I../Matrix
SRC =	solve.C
OBJS =	solve.o

all: $(EXE)

$(EXE): $(OBJS)
	$(CPP) -o $(EXE) $(OBJS) $(LIBS)

.C.o:
	$(CPP) -c $(CPPFLAGS) $(INC) $<

clean:
	rm -f $(OBJS) $(LIBNAME)

solve.o:	solve.C
