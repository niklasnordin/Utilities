#!/usr/bin/make


EXE = solveDE

CPP =	g++
LD = 	ld
CPPFLAGS = -O2
SRC =	solve.C
OBJS =	solve.o

all: $(EXE)

$(EXE): $(OBJS)
	$(CPP) -o $(EXE) $(OBJS)

.C.o:
	$(CPP) -c $(CPPFLAGS) $<

clean:
	rm -f $(OBJS) $(LIBNAME)

solve.o:	solve.C
