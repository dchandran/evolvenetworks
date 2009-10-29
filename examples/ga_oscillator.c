/****************************************************
    This file uses either one of the three network types
	included in the network evolution library to 
	evolve an oscillator.
	
	Build the cvode library first by compiling all the source files in cvode_src and 
	using ar *.o -o libcvode.a
	
	use the following to compile:
	
	gcc mtrand.c ga.c cvodesim.c ssa.c reactionNetwork.c ga_oscillator.c -lm -lcvode
	
	Uncomment one of the following pairs (which pairs?):
****************************************************/

#include "reactionNetwork.h"

/****************************************************/

void init()
{
    setNetworkTypeProbability(0,1);
    setNetworkTypeProbability(1,0);
    setNetworkTypeProbability(2,0);
    setNetworkTypeProbability(3,0);
	setDistributionOfMassActionNetwork(0.2,0.2,0.2,0.2,0.2,0.2);
    setRateConstantForMassActionNetwork(0.01,100);
    setMutationRatesForMassActionNetwork(0.5,0.25,0.25);
    setRateConstantsForProteinInteractionNetwork(0.001,100,0.1,100,0.1,100);
    setMutationRatesForProteinInteractionNetwork(0.2,0.4,0.2,0.2);
    setRateConstantsForEnzymeNetwork(0.001,100,-4,4,-4,4,1,4,0.001,20,0.001,20);
    setMutationRatesForEnzymeNetwork(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1);
    setResourceRestriction(0,0);
    setRateConstantsForGeneRegulationNetwork(1,4,0.001,100,0.1,20,0.1,10);
    setMutationRatesForGeneRegulationNetwork(0.2,0.2,0.2,0.2,0.2);
    setNetworkSize(4,12,5,24);
    setCrossoverRate(1.0);
    setAverageInitialValue(1.0);
    setMutationRateOfInitialValues(0.05);
    configureContinuousLog(1,0,1,0,0,1);
    configureFinalLog(1,1,1,0,1,1,1);
    enableLogFile("evolution.log");
    GAlineageTrackingON();
}

/* Fitness function that tests for oscillations by using correlation to a sine wave */
double fitness(GAindividual p);

#define INITIAL_POPULATION_SIZE 100
#define SUCCESSIVE_POPULATION_SIZE 50
#define NUM_GENERATIONS 20

int callback(int iter, int popSz, GApopulation P, double * fitnessArray, int *** parents)
{
	return 0;
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
	pop = evolveNetworks(INITIAL_POPULATION_SIZE, SUCCESSIVE_POPULATION_SIZE, NUM_GENERATIONS,  &fitness, &callback);  
	
	best = pop[0]; // Get the best network
	
	/******simulate the best network and write the result to a file************/
	
	N = getNumSpecies(best);                // Number of variables in the network
	
	y = simulateNetworkODE(best, 500, 1);   // Simulate

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
	int i, N;
	double peaks, troughs, * y, time, f, mXY = 0,mX = 0, mY = 0, mX2 = 0, mY2 = 0, dx = 0.01;

	N = getNumSpecies(net);
	
	time = 500.0;

	y = simulateNetworkODE(net,time,1);  //simulate

	N = getNumSpecies(net);
	
	time = 500.0;

	y = simulateNetworkODE(net,time,1);  //simulate

	f = 0;   // Calculate correlation to sine wave
	if (y != 0)
	{
		if (getValue(y,(1+N),20,N) < getValue(y,(1+N),50,N) &&
			getValue(y,(1+N),70,N) < getValue(y,(1+N),50,N))
		{
			for (i=10; (i+50) < 500; i+=50)
				if ( getValue(y,(1+N),i-10,N) < ( getValue(y,(1+N),i,N) - 1.0 ) &&
					 getValue(y,(1+N),i+10,N) < ( getValue(y,(1+N),i,N) - 1.0 ))
					f += (20 - getValue(y,(1+N),i,N))*(20 - getValue(y,(1+N),i,N));
				else
					f += 100.0;

			f = 100.0/f;
		}
		else
		{
			f = 0.0;
		}
		free(y);
	}

	return (f);
}
