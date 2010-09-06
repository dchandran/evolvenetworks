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
static GADeleteFunc deleteGAindividual = 0;
static GACloneFunc clone = 0;
static GAFitnessFunc fitness = 0;
static GACrossoverFunc crossover = 0;
static GAMutateFunc mutate = 0;
static GASelectionFunc selection = 0;
static GACallbackFunc callback = 0;
static GAPrintSummaryFunc printSummary = 0;
static GAPrintFunc printIndividual = 0;

/********************************************************
		Parameters for keeping log file
*********************************************************/

static FILE * _LOGFILE = 0;

static int _PRINT_SEEDS = 1;

static int _PRINT_EACH_FITNESS = 1;
static int _PRINT_FINAL_FITNESS = 1;

static int _PRINT_EACH_SCRIPT = 0;
static int _PRINT_FINAL_SCRIPT = 1;

static int _PRINT_EACH_SUMMARY = 1;
static int _PRINT_FINAL_SUMMARY = 1;

static int _PRINT_EACH_BEST_LINEAGE = 1;
static int _PRINT_FINAL_BEST_LINEAGE = 1;

static int _PRINT_EACH_ALL_FITNESS = 0;
static int _PRINT_FINAL_ALL_FITNESS = 1;

static int _PRINT_EACH_ALL_LINEAGE = 1;
static int _PRINT_FINAL_ALL_LINEAGE = 1;

/****************************
* options
*****************************/

static int _TRACK_PARENTS = 1;
static double _CROSSOVER_PROB = 0.5;

/*********************************************
* get and set the above function pointers
**********************************************/

void GAsetupNewStruct(GADeleteFunc f, GACloneFunc g)
{
	deleteGAindividual = f;
	clone = g;
}

void GAsetFitnessFunction(GAFitnessFunc f)
{
	fitness = f;
}

GAFitnessFunc GAgetFitnessFunction()
{
	return fitness;
}

void GAsetCallbackFunction(GACallbackFunc f)
{
	callback = f;
}

void GAsetPrintSummaryFunction(GAPrintSummaryFunc p1)
{
    printSummary = p1;
}

void GAsetPrintFunction(GAPrintFunc p2)
{
    printIndividual = p2;
}

GACallbackFunc GAgetCallbackFunction()
{
	return callback;
}

void GAsetCrossoverFunction(GACrossoverFunc f)
{
	crossover = f;
}

GACrossoverFunc GAgetCrossoverFunction()
{
	return crossover;
}

void GAsetMutationFunction(GAMutateFunc f)
{
	mutate = f;
}

GAMutateFunc GAgetMutationFunction()
{
	return mutate;
}

void GAsetSelectionFunction(GASelectionFunc f)
{
	selection = f;
}

GASelectionFunc GAgetSelectionFunction()
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
int GArouletteWheelSelection(GApopulation population, double * fitnessValues, double sumOfFitness, int popSz, int k)
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
int GAtournamentSelection(GApopulation population, double * fitnessValues, double sumOfFitness, int popSz, int n)
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
double _SELECTION_PARAMETER = 0.1;
/*! \brief Selects the individual with highest fitness
* \param array of individuals
* \param array of corresponding fitness values
* \param sum of all fitness values
* \param number of individual
* \ingroup ga
*/
int GAeliteSelection(GApopulation population, double * fitnessValues, double sumOfFitness, int popSz, int k)
{
	double p = (double)k/(double)popSz;
	return (int)(p * (double)popSz * _SELECTION_PARAMETER);
}
/*! \brief Selects an individual at random, with probability of selection based on a hyperbolic CDF of rank
 * \param GApopulation array of individuals
 * \param GApopulation array of corresponding fitness values
 * \param double sum of all fitness values
 * \param int number of individual
 * \ingroup ga
*/
int GAhyperbolicSelection(GApopulation population, double * fitnessValues, double sumOfFitness, int popSz, int k)
{
	double u = mtrand(),
			b = _SELECTION_PARAMETER * (double)popSz, a;

	a = (b+popSz)/popSz;
	return (int)( b*u/(a-u) );
}

/*! \brief set steepness parameter for the hyperbolic CDF based selection function
 * \param double steepness of the hyperbola. must be [0,1] lower number = more steep
 * \ingroup ga
*/
void GAsetParameterForSelection( double p )
{
	_SELECTION_PARAMETER = p;
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
	double temp;
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
		k = selection(currentGApopulation,fitnessArray,totalFitness,oldPopSz,i);

		x1 = currentGApopulation[k];

		if (parents)
		{
			parents[gen][i][0] = k;
			parents[gen][i][1] = 0;
		}

		if (crossover != NULL && mtrand() < _CROSSOVER_PROB)
		{
			temp = fitnessArray[k];
			fitnessArray[k] = 0;   //this is to prevent self-self crossover
			k2 = selection(currentGApopulation,fitnessArray,totalFitness,oldPopSz,i);
			fitnessArray[k] = temp;
			x2 = currentGApopulation[k2];
			x1 = crossover(x1,x2);

			if (!x1)
				x1 = currentGApopulation[k];

			if (x1 == currentGApopulation[k] || x1 == currentGApopulation[k2])
				x1 = clone(x1); //cannot allow the same x1

			if (parents)
				parents[gen][i][1] = k2;
		}
		else
		{
			x1 = clone(x1);
		}

		if (mutate != NULL)
			x1 = mutate(x1);

		if (!x1)
			x1 = clone( currentGApopulation[k] );

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

/*********************
* Log file
**********************/

void GAenableLog(FILE * file)
{
	_LOGFILE = file;
}

void GAdisableLog()
{
	GAenableLog(0);
}

void GAconfigureContinuousLog(int bestNetworkFitness,
							int bestNetworkScript,
							int bestNetworkInfo,
							int bestNetworkLineage,
							int allFitness,
							int allNetworkLineage )
{

	_PRINT_EACH_FITNESS = bestNetworkFitness > 0;
	_PRINT_EACH_SCRIPT = bestNetworkScript > 0;
	_PRINT_EACH_SUMMARY = bestNetworkInfo > 0;
	_PRINT_EACH_BEST_LINEAGE = bestNetworkLineage > 0;
	_PRINT_EACH_ALL_FITNESS = allFitness > 0;
	_PRINT_EACH_ALL_LINEAGE = allNetworkLineage > 0;
	if (_PRINT_EACH_ALL_LINEAGE || _PRINT_EACH_BEST_LINEAGE)
		_TRACK_PARENTS = 1;
}

void GAconfigureFinalLog(int bestNetworkFitness,
							int bestNetworkScript,
							int bestNetworkInfo,
							int bestNetworkLineage,
							int allFitness,
							int allNetworkLineage,
							int seeds)
{
	_PRINT_FINAL_FITNESS = bestNetworkFitness > 0;
	_PRINT_FINAL_SCRIPT = bestNetworkScript > 0;
	_PRINT_FINAL_SUMMARY = bestNetworkInfo > 0;
	_PRINT_FINAL_BEST_LINEAGE = bestNetworkLineage > 0;
	_PRINT_FINAL_ALL_FITNESS = allFitness > 0;
	_PRINT_FINAL_ALL_LINEAGE = allNetworkLineage > 0;
	_PRINT_SEEDS = seeds > 0;
	if (_PRINT_FINAL_BEST_LINEAGE || _PRINT_FINAL_ALL_LINEAGE)
		_TRACK_PARENTS = 1;
}

static void GAkeepLog(int iter, int stop, int popSz, GApopulation pop, double * fitnessArray, int *** parentsArray)
{
	unsigned long long * seeds;
	double f;
	int i,j,k,*parents, num = 10*popSz, max = 0;
	int * temp = 0, * ids = 0;
	GAindividual * p;
	//save
	int each_fitness = _PRINT_EACH_FITNESS,
		each_script = _PRINT_EACH_SCRIPT,
		each_size = _PRINT_EACH_SUMMARY,
		each_best_lineage = _PRINT_EACH_BEST_LINEAGE,
		each_all_fitness = _PRINT_EACH_ALL_FITNESS,
		each_all_lineage = _PRINT_EACH_ALL_LINEAGE;

	if (stop)
	{
		//cheat
		_PRINT_EACH_FITNESS = _PRINT_FINAL_FITNESS;
		_PRINT_EACH_SCRIPT = _PRINT_FINAL_SCRIPT;
		_PRINT_EACH_SUMMARY = _PRINT_FINAL_SUMMARY;
		_PRINT_EACH_BEST_LINEAGE = _PRINT_FINAL_BEST_LINEAGE;
		_PRINT_EACH_ALL_FITNESS = _PRINT_FINAL_ALL_FITNESS;
		_PRINT_EACH_ALL_LINEAGE = _PRINT_FINAL_ALL_LINEAGE;

		//printf("\n========final results=======\n");
		fprintf(_LOGFILE,"\n========final results=======\n");

		if (_LOGFILE && _PRINT_SEEDS)
		{
			seeds = getMTseeds();
			fprintf(_LOGFILE,"random number generator seeds: %llf,%llf,%llf,%llf\n",seeds[0],seeds[1],seeds[2],seeds[3]);
		}
	}

	if (iter == 0) //header
	{
		fprintf(_LOGFILE,"gen");
		if (_PRINT_EACH_FITNESS && !_PRINT_EACH_ALL_FITNESS)
		{
			fprintf(_LOGFILE,"\tfitness ");
		}

		if (_PRINT_EACH_ALL_FITNESS)
		{
			for (i=0; i < popSz; ++i)
			{
				fprintf(_LOGFILE,"\tfitness_%i",i);
			}
		}

		if (_PRINT_EACH_SUMMARY)
		{
			fprintf(_LOGFILE,"\tinfo");
		}

		if (_TRACK_PARENTS && (_PRINT_EACH_BEST_LINEAGE || _PRINT_EACH_ALL_LINEAGE))
		{
			fprintf(_LOGFILE,"\tparents");
		}

		fprintf(_LOGFILE,"\n");

		fprintf(_LOGFILE,"---");

		if (_PRINT_EACH_FITNESS && !_PRINT_EACH_ALL_FITNESS)
		{
			fprintf(_LOGFILE,"\t------- ");
		}

		if (_PRINT_EACH_ALL_FITNESS)
		{
			for (i=0; i < popSz; ++i)
			{
				fprintf(_LOGFILE,"\t -------- ",i);
			}
		}

		if (_PRINT_EACH_SUMMARY)
		{
			fprintf(_LOGFILE,"\t-------\t---------");
		}

		if (_TRACK_PARENTS && (_PRINT_EACH_BEST_LINEAGE || _PRINT_EACH_ALL_LINEAGE))
		{
			fprintf(_LOGFILE,"\t-------");
		}

		fprintf(_LOGFILE,"\n");
	}

	if (_TRACK_PARENTS && (_PRINT_EACH_BEST_LINEAGE || _PRINT_EACH_ALL_LINEAGE))
	{
		ids = (int*)malloc(num * sizeof(int));

		for (i=0; i < num; ++i)
			ids[i] = 0;

		for (i=0; i < popSz; ++i)
		{
			p = pop[i];

			if (!p) continue;

			parents = GAgetOriginalParents(i,iter,parentsArray);
			if (parents)
			{
				for (j=0; parents[j] > 0; ++j)
				{
					if (parents[j] >= max)
						max = parents[j];

					if (parents[j] >= num)
					{
						temp = ids;
						ids = (int*)malloc( num * 2 * sizeof(int) );

						for (k=0; k < num; ++k)
							ids[k] = temp[k];

						num *= 2;

						for (; k < num; ++k)
							ids[k] = 0;

						free(temp);
					}

					ids[ parents[j] ] += 1;
				}
			}
			else
			{
			}

			if (_PRINT_EACH_BEST_LINEAGE && !_PRINT_EACH_ALL_LINEAGE)
				break;
		}
	}

	fprintf(_LOGFILE,"%i",iter);
	if (_PRINT_EACH_FITNESS && !_PRINT_EACH_ALL_FITNESS)
	{
		f = fitnessArray[0];
		fprintf(_LOGFILE,"\t%lf",f);
	}

	if (_PRINT_EACH_ALL_FITNESS)
	{
		for (i=0; i < popSz; ++i)
		{
			f = fitnessArray[i];
			fprintf(_LOGFILE,"\t%lf",f);
		}
	}

	if (_PRINT_EACH_SUMMARY && printSummary)
	{
		fprintf(_LOGFILE,"\t");
		printSummary(_LOGFILE,pop[0]);
	}

	if (_TRACK_PARENTS && (_PRINT_EACH_BEST_LINEAGE || _PRINT_EACH_ALL_LINEAGE))
	{
		for (i=0; i < max; ++i)
		{
			fprintf(_LOGFILE,"\t%i",ids[i]);
		}
	}

	if (_PRINT_EACH_SCRIPT && printIndividual)
	{
		fprintf(_LOGFILE,"\n======== best individual =======\n");

		printIndividual(_LOGFILE,pop[0]);

		fprintf(_LOGFILE,"\n=====================\n");
	}

	fprintf(_LOGFILE,"\n");

	if (ids && (num > 0))
		free(ids);

	if (stop)
	{
		//restore
		_PRINT_EACH_FITNESS = each_fitness;
		_PRINT_EACH_SCRIPT = each_script;
		_PRINT_EACH_SUMMARY = each_size;
		_PRINT_EACH_BEST_LINEAGE = each_best_lineage;
		_PRINT_EACH_ALL_FITNESS = each_all_fitness;
		_PRINT_EACH_ALL_LINEAGE = each_all_lineage;
	}
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
void GAinit(GADeleteFunc deleteGAindividualPtr,
			GACloneFunc clonePtr,
			GAFitnessFunc fitnessPtr,
			GACrossoverFunc crossoverPtr,
			GAMutateFunc mutatePtr,
			GASelectionFunc selectionPtr,
			GACallbackFunc callbackPtr)
{
	deleteGAindividual = deleteGAindividualPtr;
	clone = clonePtr;
	fitness = fitnessPtr;
	crossover = crossoverPtr;
	mutate = mutatePtr;
	if (selectionPtr == 0)
		selection = &GAhyperbolicSelection;//&GArouletteWheelSelection;
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

	for (j=0; j < initPopSz; ++j)
		if (population[j])
			fitnessScores[j] = fitness(population[j]);
	GAsort(population,fitnessScores,0,initPopSz);
	i = 0;

	while (stop == 0) //keep going until max iterations or until callback function signals a stop
	{
		if (i == 0)  //initial population
		{
			population = GAnextGen(i, population, initPopSz, popSz, 0, fitnessScores, parents);
			free(fitnessScores);
			fitnessScores = (double*)malloc(popSz*sizeof(double));
		}
		else        //successive populations
		{
			population = GAnextGen(i, population, popSz, popSz, 0, fitnessScores, parents);
		}
		for (j=0; j < popSz; ++j)
			if (population[j])
				fitnessScores[j] = fitness(population[j]);

		GAsort(population,fitnessScores,parents[i],popSz);  //sort by fitness (Quicksort)

		if (callback != NULL)
			stop = callback(i,popSz, population,fitnessScores,parents);   //callback function can stop the GA

        if (_LOGFILE != NULL)
            GAkeepLog(i, 0, popSz, population, fitnessScores, parents);

		++i;

		if (i > numGenerations)
			stop = 1;  //max number of iterations
	}

    if (_LOGFILE != NULL)
        GAkeepLog(numGenerations, 1, popSz, population, fitnessScores, parents);

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

	if (population)
	{
		temp = population[i];
		population[i] = population[j];
		population[j] = temp;
	}
	if (parents)
	{
		p = parents[i];
		parents[i] = parents[j];
		parents[j] = p;
	}
}

// partition a[left] to a[right], assumes left < right
static int partition(double* a, int left, int right, GApopulation population, int ** parents)
{
	int i = left - 1;
	int j = right;
	while (1)
	{
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
	int * parents;

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

double GAgetCrossoverProb()
{
	return _CROSSOVER_PROB;
}

void GAsetCrossoverProb(double p)
{
	_CROSSOVER_PROB = p;

	if (_CROSSOVER_PROB < 0.0) _CROSSOVER_PROB = 0.0;
	if (_CROSSOVER_PROB > 1.0) _CROSSOVER_PROB = 1.0;
}

