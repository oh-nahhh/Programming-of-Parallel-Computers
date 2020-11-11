//libraries
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define ROOT 0 //define the root process

//methods
void swap(double *v, int i, int j);
void split(int*, int*, double, double*, int);
void quickSort(double*, int, int);
double *quickSort_par(MPI_Comm comm, double *arr, int dim, int *arrDim);
void merge( double *arr, int dim, double *newDouble, int lftPivot, double *receiveData, int myDim, int newDim);
void merge2(double *arr, int dim, double *newDouble, int lftPivot, double *receiveData, int myDim, int newDim, int rgtPivot);

//Main function
int main (int argc, char *argv[]) {

  //parameters
  int i;

  int rank;
  int nproc;

  int length;
  int N;

  double *data;

  double T1;
  double T2;
  double end;

  double *arr;
  double  *local_array;

  double *lft;
  double *rght;


  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //Read the size of the vector
  int dim = atoi(argv[1]);

  //Nbr of elements per process
  N = dim/nproc;

  //Allocate array
  arr = (double*)malloc(N*sizeof(double));

  if (rank == ROOT)
  {
    data = (double *)malloc(dim*sizeof(double));
    //Generate random numbers
    srand(time(NULL));
    for (i = 0; i < dim; i++)
    {
      data[i] =  drand48();  //Pseudo random nbr
    }

    //Initialize timer
    T1 = MPI_Wtime();
  }

  //Scatter the data into smaller vectors
  MPI_Scatter(data, N, MPI_DOUBLE, arr, N, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

  //Quicksort in serial
  quickSort(arr, ROOT, N-1);

  //Quicksort in parallel
  local_array = quickSort_par(MPI_COMM_WORLD, arr, N, &length);

  // Allocate vectors
  int *getPlace = (int*)malloc(nproc*sizeof(int));
  int *procDim  = (int*)malloc(nproc*sizeof(int));

  //Gather all of the processors
  MPI_Gather(&length, 1, MPI_INT, procDim, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  //Fill the placement vector for gatherv
  getPlace[0] = 0;
  for (i = 1; i < nproc; i++)
  {
    getPlace[i] = getPlace[i-1] + procDim[i-1];
  }

  //Gather all subvectors
  MPI_Gatherv(local_array, length, MPI_DOUBLE, data, procDim, getPlace, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

  if (rank == ROOT)
  {
    T2 = MPI_Wtime();

    printf("%.12f\n", T2-T1);

  //clear MPI
  free(local_array);
  free(data);
  free(procDim);
  free(getPlace);
  }

  MPI_Finalize();
  return 0;
}

//Serial quicksort algorithm
void quickSort(double *v, int left, int right) {
    int i;
    int last;

    if(left>=right)
        return;
    swap(v,left,(left+right)/2);
    last = left;
    for(i=left+1;i<=right;i++)
        if(v[i]<v[left])
            swap(v,++last,i);
    swap(v,left,last);
    quickSort(v,left,last-1);
    quickSort(v,last+1,right);
}

//Parallel quicksort algorithm
double *quickSort_par (MPI_Comm comm, double *arr, int dim, int *arrDim) {
  int i;

  int my_rank;
  int new_nproc;

  int myDim;
  int newDim;

  int exchange;
  int lftPivot = 0;
  int rgtPivot = 0;

  double pivot;
  double *receiveData;
  double *newDouble;

  MPI_Request request;
  MPI_Status status;

  MPI_Comm_rank(comm, &my_rank); //get my rank
  MPI_Comm_size(comm, &new_nproc); //get the new nproc

  //Case if nproc > 1
  if (new_nproc > 1) {

    if (my_rank == ROOT) {
      //choose pivot to be the mean value of the subarrays
      pivot = arr[dim/2];
    }

    //Broadcast the pivot
    MPI_Bcast(&pivot, 1, MPI_DOUBLE, ROOT, comm);

    //Split data according to the pivot
    split(&lftPivot, &rgtPivot, pivot, arr, dim);

    //Exchange halves
    if (my_rank < new_nproc/2) {

      exchange = my_rank+(new_nproc/2);

      //Send array elements higher than the pivot from left processor
      MPI_Isend(arr+lftPivot, rgtPivot, MPI_DOUBLE, exchange, 111, comm, &request);

      //Recive the lower elements from right processor
      MPI_Probe(exchange, 222, comm, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &newDim);

      //Allocater recive vector
      receiveData = (double *)malloc(newDim*sizeof(double));

      MPI_Recv(receiveData, newDim, MPI_DOUBLE, exchange, 222, comm, &status);
      MPI_Wait(&request, &status);

      //Get the new dimensions for the array
      dim = lftPivot + newDim;
      newDouble = (double *) malloc(dim*sizeof(double));

      //Merge two arrays of sorted data
      merge(arr, dim, newDouble, lftPivot, receiveData, myDim, newDim);
    }

    //Case if rank >= new nproc/2
    if (my_rank >= new_nproc/2) {

      //Exchange halves
      exchange = my_rank-(new_nproc/2);

      MPI_Isend(arr, lftPivot, MPI_DOUBLE, exchange, 222, comm, &request);

      MPI_Probe(exchange, 111, comm, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &newDim);

      receiveData = (double *)malloc(newDim*sizeof(double));
      MPI_Recv(receiveData, newDim, MPI_DOUBLE, exchange, 111, comm, &status);
      MPI_Wait(&request, &status);

      dim = newDim + rgtPivot;
      newDouble = (double *) malloc(dim*sizeof(double));

      //Merge two arrays of sorted data
      merge2(arr, dim, newDouble, lftPivot, receiveData, myDim, newDim, rgtPivot) ;
    }

    //clear
    free(arr);
    free(receiveData);

    //Create communicators
    MPI_Comm my_comm;

    //Divide into subgroups
    int divide = my_rank/(new_nproc/2);
    MPI_Comm_split(comm, divide, ROOT, &my_comm);


    //Recursively return sorted array by parallel Quicksort
    *arrDim = dim;
    return quickSort_par(my_comm, newDouble, dim, arrDim);
  } else {
     *arrDim = dim;
    return arr;
  }
}

//Method for splitting
void split (int *left, int *right, double mid, double *arr, int dim) {
  int i;
  int j;

  if (dim < 4)
  {
    for (i = 0; i < dim; i++)
    {
      if (arr[i] <= mid)
      {
        (*left)++;
      }
      else
      {
        (*right)++;
      }
    }
  }
  else
  {
  j = dim/2;

    if (arr[j] <= mid)
    {
      for (i = j; i < dim; i++)
      {
        if (arr[i] > mid)
        {
          break;
        }
      }

    *left = i;
    *right = dim-i;
    }
    else if (arr[j] > mid)
    {
      for(i = j; i >= 0; i--)
      {
        if(arr[i] <= mid)
        {
          break;
        }
      }

      *left = i+1;
      *right = dim-i-1;
    }
  }
}

//Method for swaping elements
void swap(double *v, int i, int j)
{
    double t;
    t = v[i];
    v[i] = v[j];
    v[j] = t;
}

//Methods for merging two sorted arrays
void merge(double *arr, int dim, double *newDouble, int lftPivot, double *receiveData, int myDim,int newDim)
{
      int m = 0;
      int n = 0;
      int i;

      for (i = 0; i < dim; i++)
      {
        if (m < lftPivot && n < newDim)
        {
          if (arr[m] <= receiveData[n])
          {
             newDouble[i] = arr[m++];
          }
          else
          {
             newDouble[i] = receiveData[n++];
          }
        }
        else if (m < lftPivot)
        {
       	   newDouble[i] = arr[m++];
        }
        else if (n < myDim)
        {
   	   newDouble[i] = receiveData[n++];
        }
      }
}

void merge2(double *arr, int dim, double *newDouble, int lftPivot, double *receiveData, int myDim,int newDim, int rgtPivot)
{
      int m = 0;
      int n = 0;
      int i;
      double *right = arr + lftPivot;

      for(i = 0; i < newDim+rgtPivot; i++)
      {
        if (m < newDim && n < rgtPivot)
        {
          if (receiveData[m] <= right[n])
          {
            newDouble[i] = receiveData[m++];
          }
          else
          {
            newDouble[i] = right[n++];
          }
        }
        else if (m < newDim)
        {
          newDouble[i] = receiveData[m++];
        }
        else if (n < rgtPivot)
        {
          newDouble[i] = right[n++];
        }
      }
}
