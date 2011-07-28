#!/usr/bin/make

CC = 	gcc

EXE = hello

OBJ = main.o

EXE: $(OBJ)
	$(CC) -o $(EXE) $(OBJ)

ALL: $(EXE)

.m.o: 
	$(CC) -c -o $*.o $<

clean:
	rm $(EXE) *.o
