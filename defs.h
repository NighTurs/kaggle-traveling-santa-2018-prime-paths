/******************************************************************************\
*								 Definitions							 *
\******************************************************************************/
#include <stdio.h>
#include <math.h>

#define NUM_CITIES 197769
#define BUFF_SIZE 255
#define PENALTY 1.1
#define WRITE_BUFF_SIZE 4096

typedef struct {
    double x, y;
} Point;

/* Global variables */
extern int primes[NUM_CITIES];
extern Point cities[NUM_CITIES];
extern double *coord_x, *coord_y;                // coordinates of the cities
extern int **W;                                    // weight matrix for the ATSP
extern int n_cities;                            // number of cities
extern char *prob_name;                            // name of the file for the weight matrix

/* Function declaration */
double dist(double x1, double y1, double x2, double y2);

double distCity(int a, int b);

double gpx(int *solution_blue, int *solution_red, int *offspring);

int *aloc_vectori(int lines);

double *aloc_vectord(int lines);

int **aloc_matrixi(int lines, int collums);

void desaloc_matrixi(int **Matrix, int lines);

void read_problem(const char* filename);

void rand_perm(int *inp, int *out, int size);

int weight(int i, int j);
