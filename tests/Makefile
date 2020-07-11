#
# Makefile
# haroldo, 2020-07-09 21:38
#
#`pkg-config --libs cbc cgl osi-clp clp osi coinutils`

CFLAGS=-Og -g -fsanitize=address,leak

all:c-interface-solver

c-interface-solver:c-interface-solver.o
	g++ c-interface-solver.o $(CFLAGS) -o c-interface-solver `pkg-config --libs cbc cgl osi-clp clp osi coinutils`

c-interface-solver.o:c-interface-solver.c
	gcc -c c-interface-solver.c $(CFLAGS) `pkg-config --cflags cbc cgl osi-clp clp osi coinutils` -DCBC_THREAD -o c-interface-solver.o





# vim:ft=make
#
