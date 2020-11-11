//Final Project - Game of Life Pthreads
//Niklas Bergqvist

//Libraries
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

//Constants
#define WIDTH	25 //Define the dimension of the board
#define HEIGHT	25
#define TIME	500000 //For the use of usleep

unsigned int width;
unsigned int height;

//Nbr of threads
unsigned int nthreads;

int *t_ids;

//Synchronization barrier.
pthread_barrier_t barr;

//Arrays for allocation
int **array1;
int **array2;

//Pointers for current and next state.
int **curptr, **nextptr, **temp;

//Flag for time exection test
unsigned int bflag;

//Input filename
char *filename;

struct timeval t_start, t_end;

//Functions
void initialize_board (int **);
void read_file (int **, char *);
void copy_region (int **, int **);
int alive_neighbors (int **, int , int );
void play (int **, int **, int , int );
void *do_game_of_life(void *);
void print (int **);
int arg_check (int , char *argv[]);
void print_help (void);

//Main program
int main (int argc, char *argv[]){
	int i;
	//Default values
	bflag=0;
	filename=NULL;
	nthreads=1;

	//Check arguments
	if(arg_check(argc, argv)!=0)
		return -1;

	// Allocate the arrays according to the dimensions
	array1 = (int **)malloc(width*sizeof(int *));

	for(i = 0; i < width; i++){
		array1[i] = (int *)malloc(height*sizeof(int));
	}

	array2 = (int **)malloc(width*sizeof(int *));

        for(i = 0; i < width; i++){
                array2[i] = (int *)malloc(height*sizeof(int));

        }

	//Pointers for the arrays
	curptr=array1;
	nextptr=array2;

	//Initialize the board
	initialize_board(curptr);

	//Read the starting position from file
	read_file (curptr, filename);

	//Copy cells
	copy_region(curptr, nextptr);

        //initialize nthreads
	pthread_t thr[nthreads];

	//Allocate memory for the thread
	t_ids = malloc(nthreads * sizeof(int));

	// Barrier initialization
    	pthread_barrier_init(&barr, NULL, nthreads);


	//Start the timer if test flag is on
	if(bflag){
		gettimeofday(&t_start, NULL);
	}

	// Create the threads
	for(i = 0; i < nthreads; i++){
		t_ids[i]=i;
        	pthread_create(&thr[i], NULL, &do_game_of_life, (void *)&t_ids[i]);

    	}

	//Wait for all the threads to finish
	for(i = 0; i < nthreads; i++){
      		  pthread_join(thr[i], NULL);

    	}

	//Stop the timer and print execution time if test mode
	if(bflag){
		gettimeofday(&t_end, NULL);

		printf("Time when run on %d threads: %d us\n", nthreads, ((t_end.tv_sec * 1000000 + t_end.tv_usec) - (t_start.tv_sec * 1000000 + t_start.tv_usec)));
	}
	return 0;
}

void initialize_board (int **curptr){
	int i, j;
		for (i=0; i<width; i++) for (j=0; j<height; j++)
			curptr[i][j] = 0;

}
//Function for reading the starting position from file
void read_file (int **curptr, char *name) {
	FILE	*f;
	int	i, j;
	char	s[100];

	f = fopen (name, "r");
	for (j=0; j<height; j++) {

		fgets (s, 100, f);

		for (i=0; i<width; i++){
			curptr[i][j] = s[i] == 'x'; //The cells are marked as x otherwise empty
		}
	}
	fclose (f);
}

//Function for copying the cells
void copy_region (int **curptr, int **nextptr){
	int i,j;
	for(j=0; j<height; j++)
		for(i=0; i<width; i++)
			if((i==0)|(j==0)|(j==height-1)|(i==width-1))
				nextptr[i][j]=curptr[i][j];
}

//Function to check how many alive neighbors
int alive_neighbors (int **curptr, int i, int j) {
	int row, col, count;

	count = 0;

	//Iterate through all neighbors
	for (row=-1; row<=1; row++)
		for (col=-1; col<=1; col++){
			//The self is not a neighbor
			if (row || col)
				if (curptr[i+row][j+col]) count++; //Add up the alive neighbors

			if(count>3){
				break;
				break;
			}
		}
	return count;
}

//Function for applying the rules of the game
void play (int **curptr, int **nextptr, int start, int finish)
{
	int i, j, alive;

       //Implementation of the rules for the Game of Life
	for (i=1; i<width-1; i++) for (j=start; j<finish; j++){
		alive = alive_neighbors (curptr, i, j);
                //If I am dead and I have three live neighbors, I should come alive.
		if (alive == 2) nextptr[i][j] = curptr[i][j];

	        //If I am alive and I have two or three live neighbors, I should stay alive
		if (alive == 3) nextptr[i][j] = 1;

                // In all other cases, I am dead.
		if (alive < 2) nextptr[i][j] = 0;
		if (alive > 3) nextptr[i][j] = 0;
	}
}
//Function to print the board
void print (int **curptr){
	int i, j;

	for (j=0; j<height; j++){

		for (i=0; i<width; i++){
			printf ("%c", curptr[i][j] ? 'x' : ' ');
		}
		printf ("\n");
	}
}


//Function for assigning the threads their load and running the game itself.
void *do_game_of_life(void *t_id){
	int *thread_id=(int*)t_id;

	//Assign and divide the load among the threads
	int bound = height / nthreads;
	int start = *thread_id * bound;
	int finish = start + bound;

	int i,rc;

	if(*thread_id==0) start++;
	if(*thread_id==nthreads-1) finish=height-1;

	//Run the game for i generations
	for (i=0; i<100; i++){
		play (curptr, nextptr, start, finish);
	        //Synch the threads

    		rc = pthread_barrier_wait(&barr);
    		if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD){
        		pthread_barrier_destroy(&barr);
	        	return NULL;
   		}

	       //Master process keeping track and visualizing the current generation
               //Swap pointers
		if(rc==PTHREAD_BARRIER_SERIAL_THREAD){
			if((i==0) && (!bflag)){
					system("clear");
					print (curptr);

					getchar();
					system("clear");
				}
				else{
					temp=curptr;
					curptr=nextptr;
					nextptr=temp;

					if(!bflag){
						print (curptr);
						printf("----------------------------");
						printf("Current Generation: %d",i);
						printf("---------------------------\n");
						usleep(TIME);
						system("clear");
					}
				}
		}

		//Synch the threads to make sure that the pointers have been swapped
		rc = pthread_barrier_wait(&barr);
    		if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD){
        		pthread_barrier_destroy(&barr);
	        	return NULL;
   		}
	}
	return 0;
}
//Function for checking the arguments
int arg_check(int argc, char *argv[]){
	int opt_char;
	if(argc == 1){
		print_help();
		return -1;
	}
	while ((opt_char = getopt(argc, argv, "n:h:w:f:b")) != -1) {

		 switch(opt_char) {
                        case 'n':
				nthreads=atoi(optarg);
                                break;
                        case 'f':
				filename=optarg;
                                break ;
			case 'b':
				bflag=1;
				break ;
                        default:
                                print_help();
                                break;
                }
        }

		width=WIDTH;
		height=HEIGHT;


	return 0;
}

//Print example on how to use
void print_help(){
        printf("\nExample usage:\t ./gol_threads -f Random -n 4\n");

}
