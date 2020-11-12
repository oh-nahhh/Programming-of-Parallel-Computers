F90C	= f90
CC	= g++
LIBS	=-lm
FFLAGS    = -O3 -xopenmp
CFLAGS    =  -O3 -fopenmp -Wno-write-strings

helloworld:	helloworld.c
	$(CC) $(CFLAGS) -o helloworld helloworld.c $(LIBS)

datasharing:	datasharing.c
	$(CC) $(CFLAGS) -o datasharing datasharing.c $(LIBS)

loop:	loop.c
	$(CC) $(CFLAGS) -o loop loop.c $(LIBS)

reduce:	reduce.c
	$(CC) $(CFLAGS) -o reduce reduce.c $(LIBS)

pi:	pi.c
	$(CC) $(CFLAGS) -o pi pi.c $(LIBS)

enumsort:	enumsort.c
	$(CC) $(CFLAGS) -o enumsort enumsort.c $(LIBS)

lu: lu.c
	$(CC) $(CFLAGS) -o lu lu.c $(LIBS)
main: main.c
	$(CC) $(CFLAGS) -o main main.c $(LIBS)

lulock: lulock.f90 walltime.o
	$(F90C) $(FFLAGS) -o lulock lulock.f90 walltime.o

walltime.o:	walltime.f90
	$(F90C) $(FFLAGS) -c walltime.f90


clean:
	rm -f *.o lu lulock enumsort pi reduce loop datasharing helloworld
