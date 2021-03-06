/****************************************************************************
 **
 ** Copyright (C) 2008 Deepak Chandran
 ** Contact: Deepak Chandran (dchandran1@gmail.com)
 **
 ****************************************************************************/

#ifndef DEEPAK_GILLESPIE_IMPLEMENTATION
#define DEEPAK_GILLESPIE_IMPLEMENTATION

#include <math.h>
#include <stdlib.h>
#include "mtrand.h"

#ifndef getValue
/*! \brief
* get the i,j th value from a 2D array stored as a single 1D array with N columns
*/
#define getValue(array, N, i, j) ( array[ (i)*(N) + (j) ] )
#endif

/*! \brief Stochastic simulation using Gillespie algorithm
* \param int number of species (rows of stoichiometry matrix)
* \param int number of reactions (columns of stoichiometry matrix)
* \param double* stoichiometry matrix as one array, arranged as { row1, row2, row3... }, where each row is for one species
* \param (*f)(time, y-values, rates) pointer to propensity function -- f(time, y-values, rates) --- assign values to the rates array
* \param double* initial values for species
* \param double start time
* \param double end time 
* \param int max size of array to allocate (will stop even if end time is not reached)
* \param int returns the size of the final array here
* \param void* any external data
* \return double* one dimentional array -- { row1, row2...}, where each row contains {time1,x1,x2,....}. Use getValue(y,n,i,j) if needed
* \ingroup gillespie
*/
double * SSA(int, int, double *, void (*f)(double,double*,double*,void*), double*, double, double,int, int*,void*);

/*! \brief Get rates from the simulated data
* \param double* simulated data
* \param int number of rows in the simulated data
* \param int number of species
* \param int number of reactions
* \param int the number of columns to skip (IMPORTANT: index starts at 1), i.e. first column is time, use 1. Use 0 for none.
* \param (*f)(time, y-values, rates) pointer to propensity function -- f(time, y-values, rates) --- assign values to the rates array
* \param void* any external data
* \return double* one dimentional array -- { row1, row2...}, where each row contains {time1,x1,x2,....}. Use getValue(y,n,i,j) if needed
* \ingroup gillespie
*/
double * getRatesFromSimulatedData(double* data, int rows, int cols1, int cols2, int skip, void (*f)(double,double*,double*,void*), void* param);


#endif
