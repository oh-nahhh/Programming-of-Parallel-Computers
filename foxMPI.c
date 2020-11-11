//Assignment - Fox's Algorithm MPI
//Niklas Bergqvist

//libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

//methods
void multiply(double *, double *, double *, int);

//Main
int main(int argc, char *argv[])
{

  int nproc;  //nbr of processors
  int order; //grid order

  //block sizes
  int n;
  int n_block;

  int my_rank;
  int col_rank;
  int row_rank;

  //local variables
  int i;
  int j;
  int k;

  int coords[2];

  int ndim=2;
  int dims[2]={0,0};
  int reorder=1;
  int periods[2] = {1,1};


 //local matrices
  double *A;
  double *B;
  double *C;
  double *blockA;
  double *blockB;
  double *blockC;
  double *tmpA;
  double *tmpB;


  int tmp_rank;
  int rank_source;
  int rank_dest;

  //start and stop times
  double start;
  double end;

  //MPI parameters
  MPI_Comm proc_grid;

  MPI_Request request, requesend;
  MPI_Status status2;

  MPI_Comm proc_row, proc_col;
  MPI_Status status;

  //initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc); //get the nbr of processors

  //create 2D cartesian topology
  MPI_Dims_create(nproc, ndim, dims);

  MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, reorder, &proc_grid);
  MPI_Comm_rank(proc_grid, &my_rank); //get my rank nbr

  n = atoi(argv[1]); //read in the problem size

  if (my_rank == 0)
  {
    //allocate local matrices
    A=(double *)malloc(n*n*sizeof(double));
    B=(double *)malloc(n*n*sizeof(double));
    C=(double *)malloc(n*n*sizeof(double));

    srand((unsigned) time(NULL));

    for (int i=0; i<n*n; i++)
    {
       //randomize matrices
      A[i]= ((double)rand()/(double)RAND_MAX);
      B[i]= ((double)rand()/(double)RAND_MAX);

      C[i]= 0;
    }

  }

  //initialize timer
  if (my_rank == 0)
  {
      start = MPI_Wtime();
  }

  order = sqrt(nproc); //we assume a square nbr of processes
  n_block = n/order;


  //determine the process coordinates in the cartesian topology
  MPI_Cart_coords(proc_grid, my_rank, ndim, coords);

  //create communicators for each row and column
  MPI_Comm_split(proc_grid,coords[1],coords[0],&proc_col);
  MPI_Comm_split(proc_grid,coords[0],coords[1],&proc_row);

  MPI_Comm_rank(proc_col, &col_rank);
  MPI_Comm_rank(proc_row, &row_rank);

  //allocate blocks
  blockA = (double *)calloc(n_block*n_block, sizeof(double));
  blockB = (double *)calloc(n_block*n_block, sizeof(double));
  blockC= (double *)calloc(n_block*n_block, sizeof(double));
  tmpA = (double *)calloc(n_block*n_block, sizeof(double));
  tmpB = (double *)calloc(n_block*n_block, sizeof(double));

  //create layout for the scatter and gathering process
  MPI_Datatype type, type2;
  MPI_Type_vector(n_block,n_block,n,MPI_DOUBLE,&type2);
  MPI_Type_create_resized(type2, 0, sizeof(double), &type);
  MPI_Type_commit(&type);

  int displays[nproc];
  int n_send[nproc];

 //determine the positions for the distribution of blocks in the scattering process
  for (int i=0; i<order; i++)
  {
      for (int j=0; j<order; j++)
      {
	  displays[i*order+j] = i*n*n_block+j*n_block;
	  n_send[i*order+j] = 1;
      }
  }

  //the scattering process divides the 2 matrices into blocks
  MPI_Scatterv(A, n_send, displays, type, blockA, n_block*n_block, MPI_DOUBLE, 0, proc_grid);
  MPI_Scatterv(B, n_send, displays, type, blockB, n_block*n_block, MPI_DOUBLE, 0, proc_grid);

  //copy the blocks
  memcpy(tmpB, blockB, (n_block*n_block)*sizeof(double));

  //Fox's Algorithm
  for (int k = 0; k < dims[0]; k++)
  {
    tmp_rank = (coords[0]+k) % dims[0];

    if (row_rank == tmp_rank)
    {
	memcpy(tmpA, blockA, (n_block*n_block)*sizeof(double));	//copy for broadcast
    }

    //broadcast per row
    MPI_Bcast(tmpA, (n_block*n_block), MPI_DOUBLE, tmp_rank, proc_row);

    //shift with respect to the row
    MPI_Cart_shift(proc_grid, 0, -1, &rank_source, &rank_dest);
    MPI_Isend(blockB, (n_block*n_block), MPI_DOUBLE, rank_dest, 111, proc_grid, &request);

    //multiply the matrices
    multiply(blockC,tmpA, tmpB, n_block);

    MPI_Irecv(tmpB, (n_block*n_block), MPI_DOUBLE, rank_source, 111, proc_grid, &requesend);

    MPI_Wait(&request, &status);
    MPI_Wait(&requesend, &status);
    memcpy(blockB, tmpB, (n_block*n_block)*sizeof(double));

  }

  //gather the blocks
  MPI_Gatherv(blockC, (n_block*n_block), MPI_DOUBLE, C, n_send, displays, type, 0, proc_grid);

  if (my_rank == 0)
  {
        //output time
      end = MPI_Wtime();
      printf("Total Runtime: %0.3f seconds.\n", end-start);
  }

  //Clear
  if (my_rank == 0)
  {
    free(A);
    free(B);
    free(C);
  }

  free(blockA);
  free(blockB);
  free(blockC);
  free(tmpA);
  free(tmpB);
  MPI_Type_free(&type);

  //end MPI
  MPI_Finalize();
  return 0;
}

//method for multiplication of two matrices
void multiply(double *c, double *a, double *b, int n)
{
	int i, j, k;
	for (i = 0; i < n ; i++)
	{
		for (j = 0; j < n ; j++)
		{
			double res = 0;
			for ( k = 0 ; k < n ; k++)
			{
				c[i*n+j]+= a[i*n+k]*b[k*n+j];
			}
		}
	}
}
