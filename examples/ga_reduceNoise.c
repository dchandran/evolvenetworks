/****************************************************
	This file uses either one of the three network types
	included in the network evolution library to
	evolve a network that minimizes the coefficient
	of variation at steady state

	If you already made the libode library using CMake, then use the following
	command to compile this file:

	gcc -I../simulation -I../GA -L../lib ga_reduceNoise.c -o evolveNoiseDamper -lode

	If you do not wish to use CMake, then you need to make the cvode library by compiling
	everything in the cvode_src directory and making the library using ar *.o -o libcvode.a

	After building cvode use the following to compile this file:

	gcc -I../simulation -I../GA mtrand.c ga.c cvodesim.c ssa.c reactionNetwork.c ga_reduceNoise.c -lm -lcvode

****************************************************/

#include "blocks.h"

/****************************************************/

#define INITIAL_POPULATION_SIZE 1000
#define SUCCESSIVE_POPULATION_SIZE 100
#define NUM_GENERATIONS 50

/* fitness that calculates the coefficient of variation (CV) */
double fitness(System * p);

/* print the generation number and fitness of best network during each iteration */
int callback(int iter,int popSz, GApopulation pop, double * fitnessArray, int *** );

/*Main*/
int main()
{
	int n, sz;
	double * y;
	GApopulation pop;
	System * best;

	setSizeRange(2,10);
	setMutationRate(5);
	GAsetCrossoverProb(0.5);
	GAconfigureContinuousLog(1,0,1,0,0,0);
	GAconfigureFinalLog(1,1,1,0,1,1,1);
	GAlineageTrackingON();
	GAenableLog(stdout);

	printf ("Evolving Noise-Reducing Network\n\n");

	pop = evolveNetworks(&fitness, INITIAL_POPULATION_SIZE, SUCCESSIVE_POPULATION_SIZE, NUM_GENERATIONS, &callback);

	best = (System*)pop[0];  //get the best network

	/******simulate the best network and write the result to a file************/

	n = numSpeciesTotal(best);
	y = simulateStochastic(best,500.0,&sz);  //stochastic simulation

	writeToFile("dat.txt",y,sz,n+1);  //write table to file

	free(y);

	/****** free all the networks returned by the genetic algorithm ************/

	GAfree(pop);

	return 0;
}

/* fitness that calculates the coefficient of variation (CV) */
double fitness(System * p)
{
	int i,r,n,sz;
	double f, sd, dt, time, * y, mXY = 0,mX = 0, mY = 0, mX2 = 0, mY2 = 0;

	n = numSpeciesTotal(p);
	r = numReactionsTotal(p);

	time = 500.0;

	y = simulateStochastic(p,time,&sz);  //stochastic simulation

	f = 0;
	if (y != 0)         //compute the variance
	{
		mXY = mX = mY = mX2 = mY2 = 0;
		for (i = 0; i < (sz-1); ++i)
		{
			dt = getValue(y,n+1,i+1,0) - getValue(y,n+1,i,0);
			mX += getValue(y,n+1,i,1) * dt;
			mX2 += getValue(y,n+1,i,1)*getValue(y,n+1,i,1)*dt;
		}

		mX /= time;
		mX2 /= time;

		sd = sqrt(mX2 - mX*mX);  //standard deviation

		if (sd <= 0 || mX <= 0 || mX > 5.0)
			f = 0.0;
		else
			f = mX / sd;   // CV = sdev/mean, but the fitness = 1/CV = mean/sdev

		free(y);
	}

	return (f);
}

/* print the generation number and fitness of best network during each iteration */
int callback(int iter,int popSz, GApopulation pop, double * fitnessArray, int *** parents)
{
	int i,j;
	double f = fitnessArray[0];

	if (iter > 50 && f < 0.5)
	{
		for (i=1; i < popSz; ++i)
		{
			for (j=0; j < 10; ++j)
				pop[i] = mutateSystem(pop[i]);
		}
	}
	return 0;
}


