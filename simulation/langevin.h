#include <math.h>
#include <stdlib.h>
#include "mtrand.h"

/*
* get the i,j th value from a 2D array stored as a single 1D array with N columns
*/
#ifndef getValue
#define getValue(array, N, i, j) ( array[ (i)*(N) + (j) ] )
#endif

double rnorm();

#ifndef _PropensityFunction
#define _PropensityFunction
typedef void (*PropensityFunction)(double time,double* y,double* rates,void* params);
#endif

#ifndef _EventFunction
#define _EventFunction
typedef int (*EventFunction)(int i, double time, double* y, void* params);
#endif

#ifndef _ResponseFunction
#define _ResponseFunction
typedef void (*ResponseFunction)(int i, double* y, void* params);
#endif

/*! \brief Langevin simulation
 * \ingroup gillespie
*/
double * Langevin(int n, int m, double * N, PropensityFunction propensity, double * inits, double endTime, double dt, void * params, int numEvents, EventFunction eventFunction, ResponseFunction responseFunction);
