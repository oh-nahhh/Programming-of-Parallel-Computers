//Assignment 1
//Niklas Bergqvist



// mpitest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"      //Visual Studio stuff
#include "Include\mpi.h" //Visual Studio stuff

//Normal libraries
//#include <Include\mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

//For timer in Visual Studio
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64

//Need to create timeval manually in Visual Studios, since it can't find "sys/time.h" library
typedef struct timeval {
	long tv_sec;
	long tv_usec;
} timeval;

//Stolen code to get timeofday to work in Visual Studio
int gettimeofday(struct timeval *tp, struct timezone *tzp)
{
	// Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
	// This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
	// until 00:00:00 January 1, 1970
	static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

	SYSTEMTIME  system_time;
	FILETIME    file_time;
	uint64_t    time;

	GetSystemTime(&system_time);
	SystemTimeToFileTime(&system_time, &file_time);
	time = ((uint64_t)file_time.dwLowDateTime);
	time += ((uint64_t)file_time.dwHighDateTime) << 32;

	tp->tv_sec = (long)((time - EPOCH) / 10000000L);
	tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
	return 0;
}

double timer() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
	return seconds;
}

//Main function
int main(int argc, char *argv[]) {

	// Variables of world communication
	int rank, nproc;

	// Initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Nbr of intervals
	int N;
	N = 100000L;

	// Variables for Cartesian topology processor grid
	int coords[2], pos[2], reorder = 1, ndims = 2, dims[2] = { 0,0 }, periods[2] = { 0,0 };

	// 2D-grid
	MPI_Comm proc_grid;

	// MPI status variable
	MPI_Status status;

	// MPI datatype vector
	MPI_Datatype vec_type;
	MPI_Dims_create(nproc, ndims, dims);

	// Create Cartesian grid
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &proc_grid);

	// Set my coordinates in grid
	MPI_Cart_coords(proc_grid, rank, ndims, coords);

	// Set my rank in grid
	int my_rank;
	MPI_Cart_rank(proc_grid, coords, &my_rank);

	int nnx;
	int nny;

	//Divide the load on the processes
	nnx = (N) / dims[0];
	nny = (N) / dims[1];

	//std::cout << nnx;
	//std::cout << "\n";
	//std::cout << nny;

	double dx, sum, dy, globsum, t;

	dx = 1.0/nnx;
	dy = 1.0/nny;

	sum = 0.0;

	if (rank == 0) {
		t = timer();
	}

	//and integrate
	for (int i = 0; i < nnx; i++) {
		for (int j = 0; j < nny; j++) {
			double x = dx * (i - 0.5);
			double y = dy * (j - 0.5);
			sum += dx*dy*4.0/(1.0 + x*x + y*y);
		}
	}

	sum = sum / nproc;
	MPI_Reduce(&sum, &globsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		t = timer() - t;
		printf("Integral is approx. %.16f\n", globsum);
		printf("Time: %f s\n", t);
	}

	// Close MPI
	MPI_Finalize();

	return 0;
}
