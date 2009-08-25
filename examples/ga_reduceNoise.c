#include "ga.h"

/****************************************************
This file uses either one of the three network types
included in the network evolution library to 
evolve a network that minimizes the coefficient
of variation at steady state

Build the cvode library first by compiling all the source files in cvode_src and 
using ar *.o -o libcvode.a

use the following to compile:

gcc mtrand.c ga.c cvodesim.c ssa.c reactionNetwork.c ga_reduceNoise.c -lm -lcvode

Uncomment one of the following pairs:
****************************************************/

#include "reactionNetwork.h"

/****************************************************/

/* fitness that calculates the coefficient of variation (CV) */
double fitness(void * p);

/* print the generation number and fitness of best network during each iteration */
int callback(int iter,GApopulation pop,int popSz);

/*Main*/
int main()
{	
	int i, r, n, sz;
	double * iv, * y;
	GApopulation pop;
	ReactionNetwork * net;

	setFitnessFunction( &fitness );  //set the fitness function	
	setNetworkType( GENE_REGULATION_NETWORK );  //use this network type 
	setInitialNetworkSize(4,3);  //network size

	//evolve using 1000 initial networks, 500 neworks during each successive generation, for 20 generations
	pop = evolveNetworks(200,100,20,&callback);  

	net = (ReactionNetwork*)pop[0];  //get the best network

	printNetwork(net); //print the best network

	/******simulate the best network and write the result to a file************/

	r = getNumReactions(net);
	n = getNumSpecies(net);
	iv = malloc( n * sizeof(double));
	for (i = 0; i < n; ++i) iv[i] = 0.0;

	y = simulateNetworkStochastically(net,iv,500,&sz);  //stochastic simulation
	free(iv);

	writeToFile("dat.txt",y,sz,n+1);  //write table to file

	free(y);

	/****** free all the networks returned by the genetic algorithm ************/

	for (i=0; i < 100; ++i)
		deleteNetwork(pop[i]);

	return 0;
}

/* fitness that calculates the coefficient of variation (CV) */
double fitness(void * p)
{
	int i,r,n,sz;
	double f, sd, dt, time, * iv, * y, mXY = 0,mX = 0, mY = 0, mX2 = 0, mY2 = 0;
	ReactionNetwork * net = (ReactionNetwork*)p;  //the network to test

	n = getNumSpecies(net);
	r = getNumReactions(net);

	iv = malloc( n * sizeof(double));  //initial concentrations
	for (i = 0; i < n; ++i) iv[i] = 0.0;

	time = 500.0;

	y = simulateNetworkStochastically(net,iv,time,&sz);  //stochastic simulation
	free(iv);

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
	if(getNumSpecies(net) > 5)        //disallow large networks
		f = 0.0;

	return (f);
}

/* print the generation number and fitness of best network during each iteration */
int callback(int iter,GApopulation pop,int popSz)
{
	int i,j;
	double f = fitness(pop[0]);
	ReactionNetwork * net = (ReactionNetwork*)(pop[0]);

	printf("%i\t%i\t%lf\n",iter,getNumSpecies(net),f);
	if (iter > 50 && f < 0.5)
	{
		for (i=1; i < popSz; ++i)
		{
			for (j=0; j < 10; ++j)
				pop[i] = mutateNetwork(pop[i]);
		}
	}
	return 0;
}


