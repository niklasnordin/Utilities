#!/usr/bin/make

LIBNAME = 	libMatrix.so

CPP =	g++
LD = 	ld
CPPFLAGS = 
SRC =	Matric.C Vector.C
OBJS =	Matrix.o Vector.o

all: $(LIBNAME)

$(LIBNAME): $(OBJS)
	$(CPP) -shared -o $(LIBNAME) $(OBJS)

.C.o:
	$(CPP) -c $<

clean:
	rm -f $(OBJS) $(LIBNAME)

Matrix.C:	Matrix.H Vector.H
Vector.C:	Vector.H
