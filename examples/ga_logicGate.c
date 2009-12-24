/****************************************************
    This file uses either one of the three network types
	included in the network evolution library to
	evolve a logic gate. Demonstrates the use of a predefined fitness function called
	compareSteadyStates

	If you already made the libode library using CMake, then use the following
	command to compile this file:

	gcc -I../simulation -I../GA -L../lib ga_logicGate.c -o evolveXOR -lode

	If you do not wish to use CMake, then you need to make the cvode library by compiling
	everything in the cvode_src directory and making the library using ar *.o -o libcvode.a

	After building cvode use the following to compile this file:

	gcc -I../simulation -I../GA mtrand.c ga.c cvodesim.c ssa.c reactionNetwork.c ga_logicGate.c -lm -lcvode

****************************************************/

#include "blocks.h"

/****************************************************/

#define INITIAL_POPULATION_SIZE 1000
#define SUCCESSIVE_POPULATION_SIZE 200
#define NUM_GENERATIONS 50

double ** XORtable(); //make XOR table

/* fitness function that tests for oscillations by using correlation to a sine wave */
double fitness(GAindividual p);

/* print the number of each generation and the fitness of the best network */
int callback(int iter,int popSz, GApopulation pop,double *, int ***);

/*main*/
int main()
{
	int i, N, num_rows, num_inputs, num_outputs;
	double *iv;
	GApopulation pop;
	GAindividual * best;
	double ** table, ** table2, f;

	setSizeRange(2,10);
	setMutationRate(5);
	GAsetCrossoverProb(0.5);
	GAconfigureContinuousLog(1,0,1,0,0,0);
	GAconfigureFinalLog(1,1,1,0,1,1,1);
	GAlineageTrackingON();
	GAenableLog(stdout);

	printf ("XOR Gate Evolution\n\n");

	pop = evolveNetworks(&fitness, INITIAL_POPULATION_SIZE, SUCCESSIVE_POPULATION_SIZE, NUM_GENERATIONS, &callback);

	best = pop[0]; //get the best network

	N = numSpeciesTotal(best);
	iv = (double*)malloc( N * sizeof(double) );
	for (i=0; i < N; ++i)
		iv[i] = 0.0;
	iv[0] = 0.1;
	iv[1] = 10.0;

	table = XORtable();
	table2 = XORtable();

	num_rows = 4;
	num_inputs = 2;
	num_outputs = 1;

	f = compareSteadyStates(best, table, num_rows, num_inputs, num_outputs, 1, table2);

	printf("\n%lf\n",f);

	//free table
	for (i=0; i < 4; ++i)
	{
		printf("%lf\t%lf\t%lf\n",table2[i][0],table2[i][1],table2[i][2]);
	}

	for (i=0; i < 4; ++i)
	{
		free(table[i]);
		free(table2[i]);
	}

	free(table);
	free(table2);

	/****** free all the networks returned by the genetic algorithm ************/
	GAfree(pop);

	printf("Press a key to exit\n");
	getchar();

	return 0; //done
}

/*
fitness function that tests how well the steady state outputs match the
given input/output table
*/
double fitness(GAindividual p)
{
	int i, num_rows, num_inputs, num_outputs;
	double score;

	double ** table = XORtable();

	num_rows = 4;
	num_inputs = 2;
	num_outputs = 1;

	score = compareSteadyStates(p, table, num_rows, num_inputs, num_outputs, 1, 0);

	//free table
	for (i=0; i < 4; ++i)
	{
		free(table[i]);
	}

	free(table);

	return score;
}

/* print the number of each generation and the fitness of the best network */
int callback(int iter,int popSz,GApopulation pop,double * fitnessArray, int *** parents)
{
	double f = fitnessArray[0];

	if (f > 0.9)
	{
		fitness(pop[0]);
		return 1;
	}
	return (int)(f > 0.9);
}

double** XORtable() //make XOR table
{
	int i,j;
	//XOR logic
	double XOR[][3] =
	{
	   {  0.1,  0.1,  0.0 },  //input = low low,   output = low
	   {  0.1, 10.0, 10.0 },  //input = low high,  output = high
	   { 10.0,  0.1, 10.0 },  //input = high low,  output = high
	   { 10.0, 10.0,  0.0 },  //input = high high, output = low
	};

	double ** table = (double**)malloc( 4 * sizeof(double*) );

	for (i=0; i < 4; ++i)
	{
		table[i] = (double*)malloc( 3 * sizeof(double) );
		for (j=0; j < 3; ++j)
			table[i][j] = XOR[i][j];
	}

	return table;
}
