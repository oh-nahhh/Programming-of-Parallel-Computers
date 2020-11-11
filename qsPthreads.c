//Assignment Quick Sort Pthreads
//Niklas Bergqvist

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

//Functions
void swap(double *, double *, double *);
void *sort_work(void *);
void quicksort_serial(double *, int, int );


//Constants
int threads;
volatile int count = 0; //count the nbr of active threads

pthread_mutex_t lock;

//Arguments
typedef struct {
    int start, end, n;
    double *a;
} sort_arg;

//Initialize timer function
int timer(void){
  struct timeval tv;
  gettimeofday(&tv, (struct timezone*)0);
  return (tv.tv_sec*1000000+tv.tv_usec);
}
//Main function
int main(int argc, char *argv[]){
    //Read in Length of array and nbr of threads
    int n = atoi(argv[1]), p = atoi(argv[2]);
    threads = p;

    //For use of nanosleep
    struct timespec ts;
    ts.tv_sec = 0;
    ts.tv_nsec = 1;

    //Allocate array
    double *a = (double *)malloc(n*sizeof(double));

    //Generate random numbers
    srand(time(NULL));

    //fill the array with random numbers
    int i;
    for (i = 0; i < n; i++)    {
	//a[i] = rand() % 10 + 1;
	a[i] = drand48(); //Pseudo randomizer
    }
    int time1, time2;

    time1 = timer();

    sort_arg *arg = (sort_arg *)malloc(sizeof(sort_arg));
    arg->a = a, arg->start = 0, arg->end = n-1, arg->n = 0;

    //Stop the timer and print the time result
    pthread_t t;
    count++;
    pthread_create(&t, NULL, sort_work, (void *)arg);

  while (count > 0){
	nanosleep(&ts, NULL);
    }

    time2 = timer();
    printf("Elapsed time: %f s \n",(time2-time1)/1000000.0);

    //Check if the list got sorted correctly
    for (i = 0; i < n-1; ++i){
	if (a[i] > a[i+1]){
	    printf("List did not get sorted\n\n");
	    break;
	}
    }
    pthread_mutex_destroy (&lock);
    pthread_exit(NULL);
}

//Function for swaping elements
void swap(double *a, double *b, double *tmp){
    *tmp = *b;
    *b = *a;
    *a = *tmp;
}


//Quicksort algorithm
void quicksort_serial(double *a, int p, int r){
    //Find the pivot element
    if (p < r){
	double pivot = a[r], tmp;
	int i = p - 1, j = p;

	while (j < r){
	    if (a[j] <= pivot){
		i++;
		swap(&a[i], &a[j], &tmp);
	    }
	    j++;
	}
	i++;
        //Move the pivot element to its final postition
	swap(&a[r], &a[i], &tmp);
        //Call the function recursively
	quicksort_serial(a, p, i-1);
	quicksort_serial(a, i+1, r);
    }
}

//Parallelized quicksort function
void *sort_work(void *arg){
    sort_arg *sarg = (sort_arg *)arg;
    double *a = sarg->a;
    int left = sarg->start, right = sarg->end, n = sarg->n;

    if (left < right){
	double pivot = a[right];
	int i = left - 1, j = left;
	double tmp;

	while (j < right){
	    if (a[j] <= pivot){
		i++;
		swap(&a[i], &a[j], &tmp);
	    }

	    j++;
	}
	i++;
	swap(&a[right], &a[i], &tmp);

        //if there are threads available use parallel sort
	if (n < threads){
            //Split the array int oleft and right and allocate two threads for each side of the list
	    sort_arg *left_arg = (sort_arg *)malloc(sizeof(sort_arg));
            //Initialize the boundaries for the threads
	    left_arg->a = a, left_arg->start = left, left_arg->end = i-1, left_arg->n = n+1;

	    sort_arg *right_arg = (sort_arg *)malloc(sizeof(sort_arg));
	    right_arg->a = a, right_arg->start = i+1, right_arg->end = right, right_arg->n = n+1;

            //Create the new thread
	    pthread_t newthread;
	    pthread_mutex_lock(&lock);
	    count++;

	    pthread_mutex_unlock(&lock);

	    //Start the threads
	    pthread_create(&newthread, NULL, sort_work, (void *)left_arg);
	    sort_work((void *)right_arg);
	}
        //Otherwise use serial sort
	else{
	    pthread_mutex_unlock(&lock);
	    quicksort_serial(a, left, i-1);
	    quicksort_serial(a, i+1, right);
	}
    }
    //clear threads
    free(arg);

    pthread_mutex_lock(&lock);
    count--;
    pthread_mutex_unlock(&lock);
    pthread_exit(NULL);
}
