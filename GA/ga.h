 /*!
  \file    ga.h
  \author: Deepak Chandran (dchandran1@gmail.com)
  \brief   A genetic genetic algorithm

	This library provides the basic functions for running a Genetic Algorithm (GA), but
	the user is required to setup the fitness function and related functions.
	
	The user MUST define:
		1) a struct that represents an "individual"
		2) a function that returns the fitness of an individual
		3) a mutation function that randomly alters an individual
		4) a function to free an individual
		5) a function to clone an individual
	
	The following function definition is optional but highly recommended:
		1) a crossover function to make a new individual from two individuals
	
	The following functions definitions are entirely optional:
		1) A function that selects individuals using fitness as probabilities (library provides one by default)
		2) A callback function can be used to examine or terminate the GA at any iteration
		
	The main functions are: GAinit and GArun
	
**/
#ifndef GA_MAIN_LOOP
#define GA_MAIN_LOOP

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mtrand.h"

/*! \brief An individual is represented by a user defined struct (void*) */
typedef void* GAindividual;

/*! \brief GApopulation array of individuals, where an individual is represented by a user defined struct (void*) */
typedef GAindividual* GApopulation;

/*!
    \name Required externally defined functions
    The following functions MUST be defined somewhere
	use initGA to initialize the GA with the functions
	\{
*/

/*! \brief
 * Free an individual from memory
 * \param GAindividual a single individual
 * \ingroup ga
*/
typedef void (*GADeleteFunc)(GAindividual);
/*! \brief
 * Make a copy of an individual and return the memory pointer
 * \param GAindividual target individual
 * \ingroup ga
*/
typedef GAindividual (*GACloneFunc)(GAindividual);
/*! \brief
 * Compute fitness of an individual. Fitness must be positive if default selection function is used
 * \param GAindividual target individual
 * \return double fitness of the individual (MUST BE POSITIVE if default selection function is used)
 * \ingroup ga
*/
typedef double(*GAFitnessFunc)(GAindividual);

/*! \}
    \name Mutation and crossover functions
    At least ONE of the following two functions SHOULD
	be defined, otherwise the GA will run, but nothing will evolve.
	\{
*/

/*! \brief
 * combine two individuals to generate a new individual. The function must not delete the parent individuals.
 * \param GAindividual parent individual 1. DO NOT DELETE parents
 * \param GAindividual parent individual 2. DO NOT DELETE parents
 * \return pointer to an individual (can be the same as one of the parents)
 * \ingroup ga
*/
typedef GAindividual (*GACrossoverFunc)(GAindividual, GAindividual);
/*! \brief
 * Change an individual randomly. If a new individual is created, then this function must delete (free) the old one.
 * \param GAindividual parent individual
 * \return GAindividual pointer to the same individual or a new one. If a new individual is created, 
        the original individual must be deleted inside the mutation function
 * \ingroup ga
*/
typedef GAindividual (*GAMutateFunc)(GAindividual);

/*!
  \}
  \name Optional functions
  The following two functions are entirely optional. 
  The GA library provides default Selection Function, 
  which can be overwritten
  They may or may not affect the GA performance 
  \{
*/

/*! \brief
 * Selection function. If null, then a default selection function is provided that linearly converts fitness values to probabilities
 * \param GApopulation individuals
 * \param GApopulation array of fitness values for the individuals
 * \param double total fitness (sum of all fitness values)
 * \param int number of individuals in the population
 * \param int index of the current population where this new individual will be placed
 * \return int index (in population vector) of the individual to select
 * \ingroup ga
*/
typedef int(*GASelectionFunc)(GApopulation , double * , double , int , int);
/*! \brief
 * Callback function. If not null, then this function is called during each iteration of the GA. 
 * This function can be used to terminate the GA at any step
 * \param int current generation (iteration)
 * \param int number of individuals in the population
 * \param GApopulation of individuals
 * \param double* fitness for the individuals
 * \param int*** parents for each generation. use GAgetOriginalParents or GAgetImmediateParents to get values from this array
 * \return int 0 = continue GA, 1 = stop GA. This can be used to stop the GA before it reaches max iterations
 * \ingroup ga
*/
typedef int(*GACallbackFunc)(int iter,int popSz, GApopulation,double* fitnesses,int*** parents);

/*! \}
  \name The main GA functions
  The central functions defined in this genetic algorithm library
  \{
*/

/*! \brief Initialize the GA. This function MUST be called before GArun
 * \param GADeleteFunc deletion function (cannot be 0)
 * \param GACloneFunc cloning function (cannot be 0) 
 * \param GAFitnessFunc fitness function pointer (cannot be 0)
 * \param GACrossoverFunc crossover function pointer (can be 0, but not recommended)
 * \param GAMutateFunc mutation function pointer (can be 0, but not recommended)
 * \param GASelectionFunc selection function pointer (can be 0)
 * \param GACallbackFunc callback function pointer (can be 0).
          This function can be used to monitor the GA progress or stopping the GA before reaching maximum iterations.
 * \ingroup ga
*/
void GAinit(GADeleteFunc, GACloneFunc ,GAFitnessFunc, GACrossoverFunc, GAMutateFunc, GASelectionFunc, GACallbackFunc);

/*! \brief The main GA loop. Must call GAinit before calling GArun. Uses GAnextGen to make new generation of individuals.
 * \param GApopulation initial population (array of individuals)
 * \param int number of individuals in the initial population
 * \param int number of individuals to be kept in the successive populations (will affect speed of GA)
 * \param int maximum number of generations (iterations). Callback can be used to stop the GA at any iteration.
 * \return GApopulation null terminated array of individuals (sorted by increasing fitness)
 * \ingroup ga
*/
GApopulation GArun(GApopulation,int sz0,int sz1,int maxIter);

/*! \}
  \name Available selection functions
  \{
*/

/*! \brief Selects an individual at random, with probability of selection ~ fitness
 * \param GApopulation array of individuals
 * \param GApopulation array of corresponding fitness values
 * \param double sum of all fitness values
 * \param int number of individual
 * \ingroup ga
*/
int GArouletteWheelSelection(GApopulation , double * , double , int , int );
/*! \brief Selects the best of two random individuals
 * \param GApopulation array of individuals
 * \param GApopulation array of corresponding fitness values
 * \param double sum of all fitness values
 * \param int number of individual
 * \ingroup ga
*/
int GAtournamentSelection(GApopulation , double * , double , int , int );
/*! \brief Selects the individual with highest fitness
 * \param GApopulation array of individuals
 * \param GApopulation array of corresponding fitness values
 * \param double sum of all fitness values
 * \param int number of individual
 * \ingroup ga
*/
int GAeliteSelection(GApopulation , double * , double , int , int );
/*! \brief Selects an individual at random, with probability of selection based on a hyperbolic CDF of rank
 * \param GApopulation array of individuals
 * \param GApopulation array of corresponding fitness values
 * \param double sum of all fitness values
 * \param int number of individual
 * \ingroup ga
*/
int GAhyperbolicSelection(GApopulation , double * , double , int , int );
/*! \brief set steepness parameter for the hyperbolic CDF based selection function and the strictness for the elitism selection.
 * \param double must be beteen 0 and 1. lower number = more steep/ more stringent
 * \ingroup ga
*/
void GAsetParameterForSelection( double );

/*! \}
  \name Convenience functions -- substitute for GAinit
  \{
*/
/*! \brief Initialize how to create and remove the struct defining an "individual"
 * \param GADeleteFunc function pointer (cannot be 0)
 * \param GACloneFunc function pointer (cannot be 0)
 * \ingroup ga
*/
void GAsetupNewStruct(GADeleteFunc, GACloneFunc);
/*! \brief set the fitness function for the GA
 * \param GAFitnessFunc function pointer (cannot be 0)
 * \ingroup ga
*/
void GAsetFitnessFunction(GAFitnessFunc);
/*! \brief set the crossover function for the GA
 * \param GACrossoverFunc function pointer (cannot be 0)
 * \ingroup ga
*/
void GAsetCrossoverFunction(GACrossoverFunc);
/*! \brief set the probability of a crossover event happening
 * \param double probability value between 0 and 1
 * \ingroup ga
*/
void GAsetCrossoverProb(double);
/*! \brief set the mutation function for the GA
 * \param GAMutateFunc function pointer (cannot be 0)
 * \ingroup ga
*/
void GAsetMutationFunction(GAMutateFunc);
/*! \brief set the selection function for the GA
 * \param GASelectionFunc function pointer (cannot be 0)
 * \ingroup ga
*/
void GAsetSelectionFunction(GASelectionFunc);
/*! \brief set the callback function for the GA. 
           This function can be used to monitor the GA progress or stopping the GA before reaching maximum iterations.
 * \param GACallbackFunc function pointer (can be 0)
 * \ingroup ga
*/
void GAsetCallbackFunction(GACallbackFunc);

/*! \}
  \name functions that are being used by the GA
  \{
*/

/*! \brief get the fitness function for the GA
 * \return GAFitnessFunc function pointer
 * \ingroup ga
*/
GAFitnessFunc GAgetFitnessFunction();
/*! \brief get the crossover function for the GA
 * \return GACrossoverFunc function pointer
 * \ingroup ga
*/
GACrossoverFunc GAgetCrossoverFunction();
/*! \brief get the probability of crossover happening
 * \return double probability
 * \ingroup ga
*/
double GAgetCrossoverProb();
/*! \brief get the mutation function for the GA
 * \return GAMutateFunc function pointer
 * \ingroup ga
*/
GAMutateFunc GAgetMutationFunction();
/*! \brief get the selection function for the GA
 * \return GASelectionFunc function pointer
 * \ingroup ga
*/
GASelectionFunc GAgetSelectionFunction();
/*! \brief get the callback function for the GA
 * \return GACallbackFunc function pointer (can be 0)
 * \ingroup ga
*/
GACallbackFunc GAgetCallbackFunction();
/*! \}
  \name Helper functions used by GArun.
  \{
*/

/*! \brief Generates the next population from current population
 * \param int generation
 * \param GApopulation array of individuals
 * \param int number of individual in population currently
 * \param int number of individual in the new population (returned array)
 * \param int 0 = delete old population, 1 = keep old population (warning: user must delete it later)
 * \param double* fitness scores to be calculated
 * \param int*** parents to be assigned
 * \return int new array of individual (size = 3rd parameter)
 * \ingroup ga
*/
GApopulation GAnextGen(int,GApopulation,int,int,short,double*,int***);

/*! \brief sort (Quicksort) a population by its fitness (used at the end of GArun)
 * \param GApopulation population to sort
 * \param GAFitnessFunc fitness function
 * \param double* fitness values
 * \param int** parents from a single generation (these need to be swapped too)
 * \param int size of population
 * \return void
 * \ingroup ga
*/
void GAsort(GApopulation, double *, int **, int);

/*! \}
  \name Convenience functions
  \{
*/

/*! \brief deallocate a population of individuals. 
	All populations returned by GArun will be null terminated. 
	The population is ASSUMED to be null terminated.
 * \return void
 * \ingroup ga
*/

void GAfree(GApopulation population);

/*!
  \}
  @name Related to lineage tracking
  \{
*/
/*! \brief turn on lineage tracking. This will track the parents of each individual when crossover occurs
 \ingroup ga
*/
void GAlineageTrackingON();
/*! \brief turn on lineage tracking. 
	This will prevent tracking of the parents of each individual when crossover occurs
	ReactioNetwork's parent field will be 0.
 \ingroup ga
*/
void GAlineageTrackingOFF();
/*! \brief check if lineage tracking is on
 \return int 0 if off, 1 if on
 \ingroup ga
*/
int GAisLineageTrackingOn();
/*! \brief get the IDs for all the original parents of this individual. IDs = (1 + indices of the original population)
 \param int index of an individual
 \param int generation from which this individual is selected
 \param int*** array of parents
 \return int* NULL TERMINATED array with IDs of parents. Free this array after use.
 \ingroup ga
*/
int* GAgetOriginalParents(int, int, int***);
/*! \brief get the ID for all immediate parents of this individual. IDs = (1 + indices of the previous population)
 \param int index of an individual
 \param int generation from which this individual is selected
 \param int*** array of parents
 \return int* NULL TERMINATED array with IDs of parents. Free this array after use.
 \ingroup ga
*/
int* GAgetImmediateParents(int, int, int***);

/*!
  \}
*/

#endif

