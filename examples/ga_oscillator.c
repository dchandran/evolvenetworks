#include "ga.h"

/****************************************************
    This file uses either one of the three network types
	included in the network evolution library to 
	evolve an oscillator.
	
	Build the cvode library first by compiling all the source files in cvode_src and 
	using ar *.o -o libcvode.a
	
	use the following to compile:
	
	gcc mtrand.c ga.c cvodesim.c ssa.c reactionNetwork.c ga_oscillator.c -lm -lcvode
	
	Uncomment one of the following pairs:
****************************************************/

#include "reactionNetwork.h"

/****************************************************/

/* fitness function that tests for oscillations by using correlation to a sine wave */
double fitness(void * p);

/* print the number of each generation and the fitness of the best network */
int callback(int iter,GApopulation pop,int popSz);

/*main*/
int main()
{	
	setFitnessFunction( &fitness );  //set the fitness function	
	setNetworkType( GENE_REGULATION_NETWORK );  //use this network type
	setInitialNetworkSize(4,3);  //network size
	
	//evolve using 1000 initial networks, 500 neworks during each successive generation, for 20 generations
	GApopulation pop = evolveNetworks(1000,300,30,&callback);  
	
	ReactionNetwork * net = pop[0]; //get the best network
	
	printNetwork(net); //print the best network
	
	/******simulate the best network and write the result to a file************/
	
	int i;
	int N = getNumSpecies(net);    //number of variables in the network
	double * iv = malloc( N * sizeof(double));  
	for (i = 0; i < N; ++i) iv[i] = 0.0; //initial values
	
	double * y = simulateNetworkODE(net, iv, 500, 1); //simulate
	free(iv);
	
	writeToFile("dat.txt",y,500,N+1);  //print to file
	
	free(y);
	
	/****** free all the networks returned by the genetic algorithm ************/
	for (i=0; i < 50; ++i)
		deleteNetwork(pop[i]);
		
	
	return 0; //done
}


/* fitness function that tests for oscillations by using correlation to a sine wave */
double fitness(void * p)
{
	int i;
	ReactionNetwork * net = (ReactionNetwork*)p;  //get the network
	
	int N = getNumSpecies(net);
	
	double * iv = malloc( N * sizeof(double));  //initial concentrations
	for (i = 0; i < N; ++i) iv[i] = 0.0;

	double time = 500.0;

	double * y = simulateNetworkODE(net,iv,time,1);  //simulate
	free(iv);

	double f = 0;   //calculate correlation to sine wave
	if (y != 0)
	{
		double mXY = 0,mX = 0, mY = 0, mX2 = 0, mY2 = 0;
		for (i = 0; i < time; ++i)
		{
			mX += getValue(y,N+1,i,1);
			mY += sin(i/4.0);
			mXY += sin(i/4.0) * getValue(y,N+1,i,1);
			mX2 += getValue(y,N+1,i,1)*getValue(y,N+1,i,1);
			mY2 += sin(i/4.0)*sin(i/4.0);
		}
		mX /= time;
		mY /= time;
		mXY /= time;
		mX2 /= time;
		mY2 /= time;

		f = ( (mXY - mX*mY)/( 0.01 + sqrt(mX2 - mX*mX)*sqrt(mY2 - mY*mY)) );   //correlation formula
		if (f < 0) f = -f; //negative correlation is just as good as positive (for oscillations)
		free(y);
	}
	if(getNumSpecies(net) > 30)        //disallow large networks
	  return (0);
	return (f);
}

/* print the number of each generation and the fitness of the best network */
int callback(int iter,GApopulation pop,int popSz)
{
	int i,j;
	double f = fitness(pop[0]);
	
	printf("%i\t%lf\n",iter,f);
	if (iter > 50 && f < 0.5)
	{
		for (i=1; i < popSz; ++i)
		{
			for (j=0; j < 10; ++j)
				pop[i] = mutateNetwork(pop[i]);
		}
	}
	
	if (f >= 0.5) return 1;  //stop
	
	return 0;
}
