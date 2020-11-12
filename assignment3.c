//Niklas Bergqvist
//Assigment 3 - OpenMP

//Libraries
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

//For visulaizing grid
#include <iostream>
using namespace std;


//Definitions
#define PI 3.14159265358979323846
double X0 = 0.0;
double XL = 1.0;
double Y0 = 0.0;
double YL = 1.0;

//Functions
double **allocateArray (int , int);
double **setBoundaries(double **, int , int);
double **setRHS(double **, double , double , int ,int);
double **poissonSolver(double ** (*iterate)(double **, double **, double, double, int, int),
                        double **, double **, double , double , int , int);
double **serialSOR(double **, double **, double , double , int , int);
double **redBlackSOR(double **, double **, double , double , int , int);
double **allocateArray(int, int);
double **visAllocateArray(int, int);
double f_rhs(double , double );
double residual(double **, double **, double , double , int , int);
void writeFileOutput (double **, int , int , char *);
void writeGridData(double , double , double , double ,
                    double , double ,
                    int , int , char *);
int getArgs(int , char **, int *, int *);
//Main function
int main(int argc, char *argv[]) {
    int  i, j;
    double **u, **rhs, time1, time2, **(*method)(double **, double **, double, double, int, int);

    //Read in grid size
     int sizeX = atoi(argv[1]), sizeY = atoi(argv[2]);

    //Calculate dx, dy
    double dx = (1.0)/(sizeX-1), dy = (1.0)/(sizeY-1);

/*
//Visualize 2D grid
    double ** myArray = visAllocateArray(sizeX, sizeY);

  for (int i=0; i<sizeY; ++i) {
    for (int j=0; j<sizeY; ++j)
      cout << myArray[i][j] << ' ';
    cout << endl;
  }//*/

    // Allocate arrays and initialize to zero
    u = allocateArray(sizeX, sizeY);
    rhs = allocateArray(sizeX, sizeY);

    //Set right-hand side

    rhs = setRHS(rhs, dx, dy, sizeX, sizeY);


    //Call serial or parallel SOR
    //method = &redBlackSOR;
    method = &serialSOR;

    //Take time for solver
    time1 = omp_get_wtime();

    //Call solver
    u = poissonSolver(method, u, rhs, dx, dy, sizeX, sizeY);

    //Set boundary conditions
    u = setBoundaries(u, sizeX, sizeY);

/*
  double ** myArray1 = u;

  for (int i=0; i<sizeY; ++i) {
    for (int j=0; j<sizeY; ++j)
      cout << myArray1[i][j] << ' ';
    cout << endl;
  }*/

    time2 = omp_get_wtime();

    writeGridData(X0, XL, Y0, YL, dx, dy, sizeX, sizeY, "output/gridData.txt");
    writeFileOutput(u, sizeX, sizeY, "output/output.txt");

    //Print execution time
    printf("Elapsed Time: %f seconds.\n", time2-time1);

    //Clear memory
    free(u);
    free(rhs);

    return 0;
}

//Serial implementation of SOR
double **serialSOR(double **u, double **f_rhs, double dx, double dy, int sizeX, int sizeY){
    int i,j;
    double w =1.9; //1.6 1.7 1.8 1.9 2.0;


    //Perform SOR iteration on the inner points
    for(i=1;i<sizeX-1;i++) {
        for(j=1;j<sizeY-1;j++) {
            u[i][j] = (1.0-w)*u[i][j]+(1.0/(2.0/(dx*dx)+2.0/(dy*dy)))*w*(((u[i-1][j]
                             + u[i+1][j])/(dy*dy)) +
                             (u[i][j-1] + u[i][j+1])/(dx*dx) -
                             f_rhs[i][j]);
        }

    }
 u = setBoundaries(u, sizeX, sizeY);
    return u;
}

//Parallelized SOR with red black coloring method
double **redBlackSOR(double **u, double **f_rhs, double dx, double dy, int sizeX, int sizeY){
    int i,j;
    double w = 1.9;

    //Run red sweep odd to odd, even to even for the inner points
 #pragma omp parallel for  private(i,j)
    for(i=1;i<sizeX-1;i++) {
        for(j=(i%2)+1;j<sizeY-1;j+=2) {

            u[i][j] = w*((dy*dy)*(u[i-1][j] + u[i+1][j]) +
                             (dx*dx)*(u[i][j-1] + u[i][j+1]) -
                             (dx*dx)*(dy*dy)*f_rhs[i][j])/(2.0*((dx*dx) + (dy*dy))) +
                      (1.0-w)*u[i][j];
        }
    }

    //Run the black sweep odd to even, even to odd for the inner points
   #pragma omp parallel for private(i,j)
    for(i=1;i<sizeX-1;i++) {
        for(j=((i+1)%2) + 1;j<sizeY-1;j+=2) {

            u[i][j] = w*((dy*dy)*(u[i-1][j] + u[i+1][j]) +
                             (dx*dx)*(u[i][j-1] + u[i][j+1]) -
                             (dx*dx)*(dy*dy)*f_rhs[i][j])/(2.0*((dx*dx) + (dy*dy))) +
                      (1.0-w)*u[i][j];
        }
    }
    return u;
}

//Initialize conditions on the boundary layer
double **setBoundaries(double **array, int sizeX, int sizeY){
    int i,j;

    //Left Boundary
    for(i=0;i<sizeX;i++) {
        array[i][0] = array[i][1];
    }
    //Right Boundary
    for(i=0;i<sizeX;i++) {
        array[i][sizeY-1] = array[i][sizeY-2];
    }
    //Top Boundary
    for(j=0;j<sizeY;j++) {
        array[0][j] = array[1][j];
    }
    //Bottom Boundary
    for(j=0;j<sizeY;j++) {
        array[sizeX-1][j] = array[sizeX-2][j];
    }

    return array;
}
//Functions for setting the right-hand side

double f_rhs(double x, double y){
    return sin(2*PI*x);
}
double **setRHS(double **array, double dx, double dy, int sizeX, int sizeY) {
    int i, j;

    for(i=0;i<sizeX;i++) {
        for(j=0;j<sizeY;j++) {
            array[i][j] = f_rhs(i*dx, j*dy);
        }
    }
    return array;
}
//Function for allocation
double **allocateArray (int row, int col){
  int i, j;
  double **array;

  //Allocate array storage:
  array = (double **)malloc(row * sizeof(double *));
  array[0] = (double *)malloc(row * col * sizeof(double));

  //initialize
  for (i=1;i<row;i++) {
    array[i] = array[i-1] + col;
  }
  //set array to zero
  for (i=0;i<row;i++) {
    for (j=0;j<col;j++) {
      array[i][j] = 0;
    }
  }
  return array;
}

//Poisson solver function
double **poissonSolver(double ** (*iterate)(double **, double **, double, double, int, int),
                        double **array, double **rhs, double dx, double dy, int sizeX, int sizeY){
    int t, count, interval, maxIter;
    double res, TOL;

    //Set the tolerence
    TOL = sqrt(sizeX*sizeY)*dx*dy;
    maxIter = 100*sizeX*sizeY;

    //Interval for printing the residual
    interval = 10;

    //Count the number of iterations
    count = 0;

    //Iterate until required tolerance is reached
    do{
        //Do the residual work
        res = residual(array, rhs, dx, dy, sizeX, sizeY);

        for (t=0;t<interval;t++) {
            array = (*iterate)(array, rhs, dx, dy, sizeX, sizeY);

            count += 1;
        }
        //Print current residual and iteration
        printf("iteration %d: residual = %f\t \n", count, res);

    } while (res > TOL && count <= maxIter);
    //Print the number number of iterations needed for convergence
    printf("Iterations required for convergence: %d\n", count);

    //Otherwise inform that convergence was not reached
    if(count>=maxIter) {
        printf("Convergence not reached after max nbr of iterations.\n");
    }
    return array;
}

//Function for calculating the residual
double residual(double **u, double **rhs, double dx, double dy, int sizeX, int sizeY) {
    int i, j;
    double Ax, res;
    res = 0;

    //Calculate "A*x - rho" and return it's norm
    #pragma omp parallel for private(i)
    for (i=1; i<sizeX-1; i++) {
        for (j=1; j<sizeY-1; j++) {
            Ax = ( 1.0/(dx*dx)*(u[i-1][j]-2.0*u[i][j]+u[i+1][j]) +
                     1.0/(dy*dy)*(u[i][j-1]-2.0*u[i][j]+u[i][j+1]) -
                     rhs[i][j] );
            res += Ax*Ax;
        }
    }
    return sqrt(res);
}

//function for visualizing the 2D grid
double **visAllocateArray (int row, int col){
  int i, j;
  double **array;

  array = (double **)malloc(row * sizeof(double *));
  array[0] = (double *)malloc(row * col * sizeof(double));

  for (i=1;i<col;i++) {
    array[i] = array[i-1] + row;
  }
for (i=0;i<row;i++) {
    for (j=0;j<col;j++) {
      array[i][j] = 1;
    }
  }

   //Left Boundary
    for(i=0;i<col;i++) {
        array[i][0] = 0;
    }
    //Right Boundary
    for(i=0;i<row;i++) {
        array[i][row-1] = 0;
    }
    //Top Boundary
    for(j=0;j<col;j++) {
        array[0][j] = 0;
    }
    //Bottom Boundary
    for(j=0;j<col;j++) {
        array[row-1][j] = 0;
    }

 //array[7][3] = 5;


  return array;
}


//Export grid
void writeGridData(double X0, double XL, double Y0, double YL,
                    double dx, double dy,
                    int row, int col, char *filename_GridData)
{
    FILE *outputFile;

    outputFile = fopen(filename_GridData, "w");

    // Write Data
    fprintf(outputFile, "# X - Y domain: Xstart Xend Ystart Yend : \ndomain %f %f %f %f\n", X0, XL, Y0, YL);
    fprintf(outputFile, "# step size dx, dy: \nstepsize %f %f\n", dx, dy);
    fprintf(outputFile, "# number of grid points in X and Y dimensions: \npoints %d %d\n", row, col);


    fclose(outputFile);
}

//export solution
void writeFileOutput (double **outputArray, int row, int col, char *filename_output)
{
    int i, j;
    FILE *outputFile;


    outputFile = fopen(filename_output, "w");

    for(i=0;i<row;i++) {
        for(j=0;j<col;j++) {
            fprintf(outputFile, "%f ", outputArray[i][j]);
        }
        fprintf(outputFile, "\n");
    }


    fclose(outputFile);
}
