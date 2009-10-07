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

/* Fitness function that tests for oscillations by using correlation to a sine wave */
double fitness(GAindividual p);

#define INITIAL_POPULATION_SIZE 800
#define SUCCESSIVE_POPULATION_SIZE 100
#define NUM_GENERATIONS 50

/* main */
int main()
{	
	int N;
	double * y;
	GApopulation pop;
	GAindividual * best;
	
	lineageTrackingON();
	setFitnessFunction( &fitness );    // Set the fitness function	
	
	setNetworkType( MASS_ACTION_NETWORK );  // Use this network type
	setMutationRatesForMassActionNetwork(0.5,0.2,0.2);
	setCrossoverRate(1.0);
	setDistributionOfMassActionNetwork(0.5,0.25,0.25,0.0,0.1,0.1);
	setRateConstantForMassActionNetwork(0.001,100.0);
	
	setNetworkSize(3,16,3,30);  //network size
	
	printf ("Oscillator Evolution\n\n");

	enableLogFile("log.txt");
	pop = evolveNetworks(INITIAL_POPULATION_SIZE, SUCCESSIVE_POPULATION_SIZE, NUM_GENERATIONS, 0);  
	
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
	int i, N, peaks, troughs;
	double * y, time, f, mX = 0, mX2 = 0, dx = 0.01;
	
	N = getNumSpecies(net);
	
	time = 500.0;

	y = simulateNetworkODE(net,time,1);  //simulate

	f = 0;   // Calculate correlation to sine wave
	if (y != 0)
	{
		peaks = 0;
		troughs = 0;
		for (i = 5; i < (time-3); ++i)
		{
			if ( (getValue(y,N+1,i,1) > 0.1) &&
				 (getValue(y,N+1,i,1) < 1000.0) &&
				 (getValue(y,N+1,i,1) > (dx + getValue(y,N+1,i-3,1))) &&
				 (getValue(y,N+1,i,1) > (dx + getValue(y,N+1,i+3,1))) &&
				 (getValue(y,N+1,i-3,1) < getValue(y,N+1,i-2,1)) && 
				 (getValue(y,N+1,i-2,1) < getValue(y,N+1,i-1,1)) && 
				 (getValue(y,N+1,i-1,1) < getValue(y,N+1,i,1)) && 
				 (getValue(y,N+1,i+1,1) < getValue(y,N+1,i,1)) && 
				 (getValue(y,N+1,i+2,1) < getValue(y,N+1,i+1,1)) && 
				 (getValue(y,N+1,i+3,1) < getValue(y,N+1,i+2,1))
				)
			{
				 ++peaks;
				 mX += y[i];
				 mX2 += y[i]*y[i];
			}

			if ( (getValue(y,N+1,i,1) > 0.1) &&
				 (getValue(y,N+1,i,1) < 1000.0) &&
				 (getValue(y,N+1,i,1) < (getValue(y,N+1,i-3,1) - dx)) &&
				 (getValue(y,N+1,i,1) < (getValue(y,N+1,i+3,1) - dx)) &&
				 (getValue(y,N+1,i-3,1) > getValue(y,N+1,i-2,1)) && 
				 (getValue(y,N+1,i-2,1) > getValue(y,N+1,i-1,1)) && 
				 (getValue(y,N+1,i-1,1) > getValue(y,N+1,i,1)) && 
				 (getValue(y,N+1,i+1,1) > getValue(y,N+1,i,1)) && 
				 (getValue(y,N+1,i+2,1) > getValue(y,N+1,i+1,1)) && 
				 (getValue(y,N+1,i+3,1) > getValue(y,N+1,i+2,1))
				)
			{
				 ++troughs;
			}
		}
		
		if (peaks < 1 || ((mX2 - mX*mX) < 0.0))
			f = 0.0;
		else
		{
			if ((troughs+peaks) > 10)
			{
				f = (double)10.0 + 1.0/(1.0 + mX2 - mX*mX);
			}
			else
			{
				f = (double)troughs + (double)peaks;
			}
		}
		free(y);
	}

	return (f);
}
