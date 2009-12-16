/****************************************************
    This file uses either one of the three network types
	included in the network evolution library to
	evolve an oscillator.

	If you already made the libode library using CMake, then use the following
	command to compile this file:

	gcc -I../simulation -I../GA -L../lib ga_oscillator.c -o evolveOscillator -lode

	If you do not wish to use CMake, then you need to make the cvode library by compiling
	everything in the cvode_src directory and making the library using ar *.o -o libcvode.a

	After building cvode use the following to compile this file:

	gcc -I../simulation -I../GA mtrand.c ga.c cvodesim.c ssa.c reactionNetwork.c ga_oscillator.c -lm -lcvode

****************************************************/

#include "blocks.h"

/****************************************************/

void init()
{
    setSizeRange(2,10);
    GAsetCrossoverRate(0.5);
    GAconfigureContinuousLog(1,0,1,0,0,1);
    GAconfigureFinalLog(1,1,1,0,1,1,1);
    GAenableLogFile("evolution.log");
    GAlineageTrackingON();

    GAenableLog(stdout);
}

/* Fitness function that tests for oscillations by using correlation to a sine wave */
double fitness(GAindividual p);

#define INITIAL_POPULATION_SIZE 1000
#define SUCCESSIVE_POPULATION_SIZE 200
#define NUM_GENERATIONS 80

int callback(int iter,int popSz, GApopulation P, double * fitnessArray, int *** parents)
{
	return (fitnessArray[0] > 0.45);
}

/* main */
int main()
{
	int N;
	double * y;
	GApopulation pop;
	GAindividual * best;

	init();

	printf ("Oscillator Evolution\n\n");

	pop = evolveNetworks(&fitness, int INITIAL_POPULATION_SIZE, SUCCESSIVE_POPULATION_SIZE, NUM_GENERATIONS, &callback);

	best = pop[0]; // Get the best network

	/******simulate the best network and write the result to a file************/

	N = getNumSpecies(best);                // Number of variables in the network

	y = simulateODE(best, 500, 1);   // Simulate

	writeToFile("dat.txt",y,500,N+1);       // Save to file

	free(y);

	/****** free all the networks returned by the genetic algorithm ************/
	GAfree(pop);

	printf("Press a key to exit\n");
	getchar();

	return 0; // Done
}

/* Fitness function that tests for oscillations by counting the number of peaks*/
double fitness(GAindividual net)
{
	int i, N, n;
	double x, * y, time, f, mXY = 0,mX = 0, mY = 0, mX2 = 0, mY2 = 0, dx = 0.01;

	N = getNumSpecies(net);

	time = 500.0;

	y = simulateNetworkODE(net,time,1);  //simulate

	f = 0;   // Calculate correlation to sine wave
	if (y != 0)
	{
		n = 0;
		mXY = mX = mY = mX2 = mY2 = 0;
		for (i = 0; i < (time); ++i)
		{
			x = sin((double)i/4.0) > 0;
			mX += getValue(y,N+1,i,1);
			mY += x;
			mXY += x * getValue(y,N+1,i,1);
			mX2 += getValue(y,N+1,i,1)*getValue(y,N+1,i,1);
			mY2 += x*x;
			++n;
		}

		mX /= (double)(n);
		mY /= (double)(n);
		mXY /= (double)(n);
		mX2 /= (double)(n);
		mY2 /= (double)(n);

		if (((mX2 - mX*mX)) < 0.0001)
		{
			f = 0.0;
		}
		else
		{
			f = ( (mXY - mX*mY)/(sqrt(mX2 - mX*mX)*sqrt(mY2 - mY*mY)) );   // Correlation formula
			if (f < 0) f = -f; // Negative correlation is just as good as positive (for oscillations)
		}
		free(y);
	}

	return (f);
}
