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
double fitness(GAindividual p);

/* print the number of each generation and the fitness of the best network */
int callback(int iter,GApopulation pop,int popSz);

FILE * lineageFile;

/*main*/
int main()
{	
	int i, N;
	double * iv, * y;
	GApopulation pop;
	ReactionNetwork * net;
	
	lineageFile = fopen("lineage.txt","w");
	setFitnessFunction( &fitness );  //set the fitness function	
	
	//setNetworkType( MASS_ACTION_NETWORK );  //use this network type
	setNetworkType( PROTEIN_INTERACTION_NETWORK );  //use this network type
	//setNetworkType( GENE_REGULATION_NETWORK );  //use this network type
	
	setInitialNetworkSize(6,2);  //network size
	
	//evolve using 1000 initial networks, 200 neworks during each successive generation, for 20 generations
	pop = evolveNetworks(1000,300,20,&callback);  
	
	net = pop[0]; //get the best network
	
	printNetworkToFile(net,"network.txt"); //print the best network
	
	/******simulate the best network and write the result to a file************/
	
	N = getNumSpecies(net);    //number of variables in the network
	iv = malloc( N * sizeof(double));  
	for (i = 0; i < N; ++i) iv[i] = 0.0; //initial values
	
	y = simulateNetworkODE(net, iv, 500, 1); //simulate
	free(iv);
	
	writeToFile("dat.txt",y,500,N+1);  //print to file
	
	free(y);
	
	/****** free all the networks returned by the genetic algorithm ************/
	for (i=0; i < 50; ++i)
		deleteNetwork(pop[i]);
	
	fclose(lineageFile);

	return 0; //done
}

void printLineage(GApopulation pop, int popSz, int num)
{
	int i,j;
	int * ids = malloc(num * sizeof(int)), * ids2 = malloc(num * sizeof(int));
	ReactionNetwork * r;

	for (i=0; i < num; ++i)
		ids[i] = 0;

	for (i=0; i < popSz; ++i)
	{
		r = (ReactionNetwork*)pop[i];
		
		if (r->parents)
		{
			for (j=0; j < num; ++j)
				ids2[j] = 0;

			for (j=0; r->parents[j] != 0; ++j)
				if (r->parents[j] < num)
					ids2[ r->parents[j] ] = 1;
			
			for (j=0; j < num; ++j)
				ids[ j ] += ids2[ j ];
		}
		else
		{
			if (r->id < num)
				ids[r->id]++;
		}
	}
	for (i=0; i < num; ++i)
		if (i==0)
		{
			fprintf(lineageFile,"%i",ids[i]);
		}
		else
		{
			fprintf(lineageFile,",%i",ids[i]);
		}
	fprintf(lineageFile,"\n");
	free(ids);
	free(ids2);
}

/* fitness function that tests for oscillations by using correlation to a sine wave */
double fitness(GAindividual p)
{
	int i, N;
	double * y, * iv, time, f, mXY = 0,mX = 0, mY = 0, mX2 = 0, mY2 = 0;
	ReactionNetwork * net = (ReactionNetwork*)p;  //get the network
	
	N = getNumSpecies(net);
	
	iv = malloc( N * sizeof(double));  //initial concentrations
	for (i = 0; i < N; ++i) iv[i] = 0.0;

	time = 500.0;

	y = simulateNetworkODE(net,iv,time,1);  //simulate
	free(iv);

	f = 0;   //calculate correlation to sine wave
	if (y != 0)
	{
		mXY = mX = mY = mX2 = mY2 = 0;
		
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

	printf("%i\t%lf\t%i\n",iter,f,getNumSpecies((ReactionNetwork*)(pop[0])));
	
	if (iter > 1 && (iter % 10 == 0) && f < 0.5)
	{
		for (i=1; i < popSz; ++i)
		{
			for (j=0; j < 10; ++j)
				pop[i] = mutateNetwork(pop[i]);
		}
	}
	
	printLineage(pop,popSz,1000);

	if (f >= 0.5) return 1;  //stop
	
	return 0;
}
