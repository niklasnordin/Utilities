#!/usr/bin/make

LIBNAME = 	libMatrix.so
STATICLIB = 	libMatrix.a

CPP =	g++
LD = 	ld
CPPFLAGS = -fPIC
SRC =	Matric.C Vector.C
OBJS =	Matrix.o Vector.o

all: $(LIBNAME) $(STATICLIB)

$(STATICLIB): $(OBJS)
	ar -r $(STATICLIB) $(OBJS)

$(LIBNAME): $(OBJS)
	$(CPP) -shared -o $(LIBNAME) $(OBJS)

.C.o:
	$(CPP) -c $(CPPFLAGS) $<

clean:
	rm -f $(OBJS) $(LIBNAME)

Matrix.C:	Matrix.H Vector.H
Vector.C:	Vector.H
