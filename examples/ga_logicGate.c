/****************************************************
    This file uses either one of the three network types
	included in the network evolution library to 
	evolve a logic gate. Demonstrates the use of a predefined fitness function called
	compareSteadyStates
	
	Build the cvode library first by compiling all the source files in cvode_src and 
	using ar *.o -o libcvode.a
	
	use the following to compile:
	
	gcc mtrand.c ga.c cvodesim.c ssa.c reactionNetwork.c ga_logicGate.c -lm -lcvode
	
	Uncomment one of the following pairs:
****************************************************/

#include "reactionNetwork.h"

/****************************************************/

double ** XORtable(); //make XOR table

/* fitness function that tests for oscillations by using correlation to a sine wave */
double fitness(GAindividual p);

/* print the number of each generation and the fitness of the best network */
int callback(int iter,GApopulation pop,int popSz);

/*main*/
int main()
{	
	int i, N;
	double *iv, * y;
	GApopulation pop;
	GAindividual * best;
	
	setFitnessFunction( &fitness );  //set the fitness function	
	
	setNetworkType( MASS_ACTION_NETWORK );  //use this network type
	
	setInitialNetworkSize(5,8);  //network size

	printf("generation\tbest fitness\tnetwork size\n");
	
	//evolve using 1000 initial networks, 200 neworks during each successive generation, for 20 generations
	pop = evolveNetworks(1000,200,10,&callback);  
	
	best = pop[0]; //get the best network
	
	printNetwork(best); //print the best network
	
	printNetworkToFile(best,"network.txt"); //print the best network
	
	N = getNumSpecies(best);
	iv = (double*)malloc( N * sizeof(double) );
	for (i=0; i < N; ++i)
		iv[i] = 0.0;
	iv[0] = 0.1;
	iv[1] = 10.0;
	
	y = networkSteadyState(best,iv);	
	
	if (y)
	{
		for (i=0; i < N; ++i)
			printf("%lf\t",y[i]);
		free(y);
	}
	else
		printf("no ss\n");
	
	/****** free all the networks returned by the genetic algorithm ************/
	for (i=0; i < 50; ++i)
		deleteNetwork(pop[i]);
	
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
	
	score = compareSteadyStates(p, table, num_rows, num_inputs, num_outputs);
	
	//free table
	for (i=0; i < 4; ++i)
	{
		free(table[i]);
	}
	
	free(table);
	
	return score;
}

/* print the number of each generation and the fitness of the best network */
int callback(int iter,GApopulation pop,int popSz)
{
	double f = fitness(pop[0]);
	
	printf("%i\t%lf\n",iter,f);
	
	if (f >= 0.9) return 1;  //stop if good enough
	
	return 0;
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
