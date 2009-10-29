/****************************************************************************
**
** Copyright (C) 2008 Deepak Chandran
** see ga.h
**
****************************************************************************/

#include "ga.h"
#include <stdio.h>

/****************************
* functions required for the GA
****************************/
static GADeleteFnc deleteGAindividual = 0;
static GACloneFnc clone = 0;
static GAFitnessFnc fitness = 0;
static GACrossoverFnc crossover = 0;
static GAMutateFnc mutate = 0;
static GASelectionFnc selection = 0;
static GACallbackFnc callback = 0;

/****************************
* lineage tracking on or  off
*****************************/

static int _TRACK_PARENTS = 1;

/*********************************************
* get and set the above function pointers
**********************************************/

void GAsetupNewStruct(GADeleteFnc f, GACloneFnc g)
{
	deleteGAindividual = f;
	clone = g;
}

void GAsetFitnessFunction(GAFitnessFnc f)
{
	fitness = f;
}

GAFitnessFnc GAgetFitnessFunction()
{
	return fitness;
}

void GAsetCallbackFunction(GACallbackFnc f)
{
	callback = f;
}

GACallbackFnc GAgetCallbackFunction()
{
	return callback;
}

void GAsetCrossoverFunction(GACrossoverFnc f)
{
	crossover = f;
}

GACrossoverFnc GAgetCrossoverFunction()
{
	return crossover;
}

void GAsetMutationFunction(GAMutateFnc f)
{
	mutate = f;
}

GAMutateFnc GAgetMutationFunction()
{
	return mutate;
}

void GAsetSelectionFunction(GASelectionFnc f)
{
	selection = f;
}

GASelectionFnc GAgetSelectionFunction()
{
	return selection;
}

/*************************************************/

/*! \brief Selects an individual at random, with probability of selection ~ fitness
* \param array of individuals
* \param array of corresponding fitness values
* \param sum of all fitness values
* \param number of individual
* \ingroup ga
*/
int GArouletteWheelSelection(GApopulation population, double * fitnessValues, double sumOfFitness, int popSz)
{
	int i;
	double randNum = mtrand() * sumOfFitness, 
		total = 0;
	for (i=0; i < popSz-1; ++i)
		if (total < randNum && randNum < (total+fitnessValues[i]))
			return (i);
		else
			total += fitnessValues[i];
	return (i);
}
/*! \brief Selects the best of two random individuals
* \param array of individuals
* \param array of corresponding fitness values
* \param sum of all fitness values
* \param number of individual
* \ingroup ga
*/
int GAtournamentSelection(GApopulation population, double * fitnessValues, double sumOfFitness, int popSz)
{
	int i = (int)(mtrand() * popSz);
	int j = (int)(mtrand() * popSz);
	int k;

	while (popSz != 1 && i==j)  //prevent tournament with self
		j = (int)(mtrand() * popSz);

	if (fitnessValues[i] > fitnessValues[j])
		k = i;
	else
		k = j;
	//fitnessValues[k] = 0;   //do not pick this individual again?
	return (k);
}
/*! \brief Selects the individual with highest fitness
* \param array of individuals
* \param array of corresponding fitness values
* \param sum of all fitness values
* \param number of individual
* \ingroup ga
*/
int GAeliteSelection(GApopulation population, double * fitnessValues, double sumOfFitness, int popSz)
{
	int i, best = 0;
	for (i=0; i < popSz; ++i)
		if (fitnessValues[i] > fitnessValues[best])
			best = i;
	fitnessValues[best] = 0;
	return (best);
}
/*
* Get next population from current population
* \param: array of individuals
* \param: number of individual in population currently
* \param: number of individual in the new population (returned array)
* \param: 0 = delete old population, 1 = keep old population (warning: user must delete it later)
* @ret: new array of individual (size = 3rd parameter)
*/
GApopulation GAnextGen(int gen, GApopulation currentGApopulation, int oldPopSz, int newPopSz,
					   short keepOldGApopulation, double * fitnessArray, int*** parents)
{
	int i,k,best,k2;
	void * x1 = NULL, * x2 = NULL;
	double totalFitness;
	GApopulation nextGApopulation;

	//allocate memory for next generation

	nextGApopulation = (void*)malloc( (1+newPopSz) * sizeof(void*) );
	nextGApopulation[newPopSz] = 0;

	if (nextGApopulation == NULL) 
	{
		return (0);
	}

	//make array of fitness values
	//fitnessArray = (double*) malloc ( oldPopSz * sizeof(double) );
	totalFitness = 0;
	best = 0;  //save best's index

	for (i = 0; i < oldPopSz; ++i)
	{
		fitnessArray[i] = fitness(currentGApopulation[i]);
		if (fitnessArray[i] < 0) fitnessArray[i] = 0;   //negative fitness not allowed

		totalFitness += fitnessArray[i];
		if (fitnessArray[i] > fitnessArray[best]) 
			best = i;
	}

	//keep the best
	nextGApopulation[0] = clone(currentGApopulation[best]);
	
	if (parents)
	{
		parents[gen][0][0] = best;
		parents[gen][0][1] = 0;
	}
	
	//select the fit individuals
	x1 = NULL;
	x2 = NULL;
	for (i = 1; i < newPopSz; ++i)
	{
		k = selection(currentGApopulation,fitnessArray,totalFitness,oldPopSz);

		x1 = currentGApopulation[k];
		
		if (parents)
		{
			parents[gen][i][0] = k;
			parents[gen][i][1] = 0;
		}
		
		if (crossover != NULL) 
		{
			double temp = fitnessArray[k];
			fitnessArray[k] = 0;   //this is to prevent self-self crossover
			k2 = selection(currentGApopulation,fitnessArray,totalFitness,oldPopSz);
			fitnessArray[k] = temp;
			x2 = currentGApopulation[k2];
			x1 = crossover(x1,x2);

			if (x1 == currentGApopulation[k] || x1 == currentGApopulation[k2])
			{
				x1 = clone(x1); //cannot allow the same x1
			}
			
			if (parents)
			{
				parents[gen][i][1] = k2;
			}
		}
		else
		{
			x1 = clone(x1);
		}

		if (mutate != NULL) 
		{
			x1 = mutate(x1);
		}

		nextGApopulation[i] = x1; //add to the new population
	}
	/*free the memory from the old population*/
	if (keepOldGApopulation == 0)
	{
		for (i = 0; i < oldPopSz; ++i)
		{
			if (currentGApopulation[i] != NULL)
				deleteGAindividual(currentGApopulation[i]);
		}

		free(currentGApopulation);
	}
	return (nextGApopulation);
}

/*
* Initialize the GA. This function MUST be called before GArun
* \param: cloning function (cannot be 0)
* \param: deletion function (cannot be 0)
* \param: fitness function pointer (cannot be 0)
* \param: crossover function pointer (can be 0, but not recommended)
* \param: mutation function pointer (can bt 0, but not recommended)
* \param: selection function pointer (can be 0)
*/
void GAinit(GADeleteFnc deleteGAindividualPtr, 
			GACloneFnc clonePtr,
			GAFitnessFnc fitnessPtr, 
			GACrossoverFnc crossoverPtr, 
			GAMutateFnc mutatePtr, 
			GASelectionFnc selectionPtr,
			GACallbackFnc callbackPtr)
{
	deleteGAindividual = deleteGAindividualPtr;
	clone = clonePtr;
	fitness = fitnessPtr;
	crossover = crossoverPtr;
	mutate = mutatePtr;
	if (selectionPtr == 0)
		selection = &(GArouletteWheelSelection);
	else
		selection = selectionPtr;
	callback = callbackPtr;
}

/*
* The main GA loop
* \param: array of individuals
* \param: number of individual initially
* \param: number of individual in successive populations
* \param: total number of generations
* \param: callback function pointer
* @ret: final array of individuals (sorted)
*/
GApopulation GArun(GApopulation initialGApopulation, int initPopSz, int popSz, int numGenerations)
{
	int i = 0, j = 0, stop = 0;
	GApopulation population = initialGApopulation;
	int*** parents = 0;
	double* fitnessScores = 0;
	
	FILE * errfile = freopen("GArun_errors.log", "w", stderr);
	
	/*function pointers*/
	if (!deleteGAindividual || !clone || !fitness || (!crossover && !mutate) || !selection) return 0;

	if (!MTrandHasBeenInitialized())
		initMTrand(); /*initialize seeds for MT random number generator*/
		
	fitnessScores = (double*)malloc(initPopSz * sizeof(double));
	
	if (_TRACK_PARENTS)
	{
		parents = (int***)malloc((1+numGenerations) * sizeof(int**));
		for (i=0; i <= numGenerations; ++i)
		{
			parents[i] = (int**)malloc((1+popSz) * sizeof(int*));
			parents[i][popSz] = 0;
			
			for (j=0; j < popSz; ++j)
			{
				parents[i][j] = (int*)malloc(3 * sizeof(int));
				parents[i][j][0] = 0;
				parents[i][j][1] = 0;
				parents[i][j][2] = 0;
			}
		}
	}
	
	i = 0;
	while (stop == 0) //keep going until max iterations or until callback function signals a stop
	{
		if (i == 0)  //initial population
			population = GAnextGen(i, population, initPopSz, popSz, 0, fitnessScores, parents);
		else        //successive populations
			population = GAnextGen(i, population, popSz, popSz, 0, fitnessScores, parents);
		
		GAsort(population,fitnessScores,parents[i],popSz);  //sort by fitness (Quicksort)

		if (callback != NULL)
			stop = callback(i,popSz, population,fitnessScores,parents);   //callback function can stop the GA

		++i;

		if (i > numGenerations) 
			stop = 1;  //max number of iterations
	}
	
	if (population[popSz-1])
	{
		deleteGAindividual(population[popSz-1]);
		population[popSz-1] = 0;  //null terminate
	}

	fclose(errfile);
	
	if (parents)
	{
		for (i=0; i <= numGenerations; ++i)
		{
			for (j=0; j < popSz; ++j)
			{
				free(parents[i][j]);
			}
			free(parents[i]);
		}
		free(parents);
	}
	free(fitnessScores);

	return (population);
}

void GAfree(GApopulation pop)
{
	int i=0;
	if (pop)
	{
		while (pop[i]) 
		{
			deleteGAindividual(pop[i]);
			++i;
		}
	}
}


/***********************************************************************
*  Quicksort code from Sedgewick 7.1, 7.2.
***********************************************************************/

// is x < y ?
static int less(double x, double y) 
{
	return (x > y);
}

// exchange a[i] and a[j]
static void exch(double* a, int i, int j, GApopulation population, int ** parents) 
{
	int * p;
	void * temp;
	double swap;
	
	swap = a[i];
	a[i] = a[j];
	a[j] = swap;

	temp = population[i];
	population[i] = population[j];
	population[j] = temp;
	
	p = parents[i];
	parents[i] = parents[j];
	parents[j] = p;
}

// partition a[left] to a[right], assumes left < right
static int partition(double* a, int left, int right, GApopulation population, int ** parents) 
{
	int i = left - 1;
	int j = right;
	while (1) {
		while (less(a[++i], a[right]))      // find item on left to swap
			;                               // a[right] acts as sentinel
		while (less(a[right], a[--j]))      // find item on right to swap
			if (j <= left) break;           // don't go out-of-bounds
		if (i >= j) break;                  // check if pointers cross
		exch(a, i, j, population, parents);         // swap two elements into place
	}
	exch(a, i, right, population, parents);                      // swap with partition element
	return i;
}

// quicksort helper a[left] to a[right]
static void quicksort(double* a, int left, int right, GApopulation population, int ** parents) 
{
	int i = partition(a, left, right, population, parents);

	if (right <= left) return;
	if (left < (i-1))
		quicksort(a, left, i-1, population, parents);
	if ((i+1) < right)
		quicksort(a, i+1, right, population, parents);
}

//quicksort
void GAsort(GApopulation population, double * a, int** parents, int populationSz) 
{
	/*
	int i;
	double biggest = -1.0;
	GAindividual best = 0;
	
	for (i=0; i < populationSz; ++i)
	{
		if (a[i] > biggest || best == 0)
		{
			best = population[i];
			biggest = a[i];
		}
	}*/
	if (a != NULL)
	{
		quicksort(a, 0, populationSz - 1, population, parents);
	}
//	population[0] = best;
}

/*******************************
  Related to lineage tracking
*******************************/

void GAlineageTrackingON()
{
	_TRACK_PARENTS = 1;
}

void GAlineageTrackingOFF()
{
	_TRACK_PARENTS = 0;
}

int GAisLineageTrackingOn()
{
	return _TRACK_PARENTS;
}

int* GAgetOriginalParents(int individual, int generation, int *** _PARENTS)
{
	int maxSz = 1, i=0, j=0;
	int sz = 0;
	int * parents, * clone;
	
	if (_PARENTS == 0 || generation < 0 || individual <= 0)
	{
		parents = (int*)malloc(1 * sizeof(int));
		parents[0] = 0;
		return parents;
	}
	
	maxSz = 100;
	parents = (int*)malloc(maxSz * sizeof(int));
	sz = 0;
	parents[sz] = 0;
	if (_PARENTS[generation][individual][0])
	{
		parents[sz] = _PARENTS[generation][individual][0];
		++sz;
	}
	if (_PARENTS[generation][individual][1])
	{
		parents[sz] = _PARENTS[generation][individual][1];
		++sz;
	}
	
	i = 0;
	while (i < sz && parents[i] && generation > 0)
	{
		--generation;
		
		if (sz >= (maxSz-1))
		{
			clone = parents;
			parents = (int*)malloc((2*sz) * sizeof(int));
			for (j=0; j < sz; ++j)
				parents[j] = clone[j];
			free(clone);
			maxSz = sz*2;
		}
		
		individual = parents[i];
		if (_PARENTS[generation][individual][0])
		{
			parents[sz] = _PARENTS[generation][individual][0];
			++sz;
		}
		
		if (_PARENTS[generation][individual][1])
		{
			parents[sz] = _PARENTS[generation][individual][1];
			++sz;
		}
	}
	
	clone = parents;
	parents = (int*)malloc((1+sz) * sizeof(int));
	parents[sz] = 0;
	for (j=0; j < sz; ++j)
		parents[j] = clone[j];
	
	free(clone);
	return parents;
}

int* GAgetImmediateParents(int individual, int generation, int *** _PARENTS)
{
	int p1=0, p2=0;
	int * parents, * clone;

	if (_PARENTS == 0 || generation < 0 || individual <= 0)
	{
		parents = (int*)malloc(1 * sizeof(int));
		parents[0] = 0;
		return parents;
	}
	
	p1 = _PARENTS[generation][individual][0];
	p2 = _PARENTS[generation][individual][1];
	
	if (p1 > 0 && p2 > 0)
	{
		parents = (int*)malloc(3 * sizeof(int));
		parents[0] = p1;
		parents[1] = p2;
		parents[2] = 0;
		return parents;
	}
	
	if (p1 > 0)
	{
		parents = (int*)malloc(2 * sizeof(int));
		parents[0] = p1;
		parents[1] = 0;
		return parents;
	}
	
	if (p2 > 0)
	{
		parents = (int*)malloc(2 * sizeof(int));
		parents[0] = p2;
		parents[1] = 0;
		return parents;
	}
	
	parents = (int*)malloc(1 * sizeof(int));
	parents[0] = 0;
	
	return parents;
}
