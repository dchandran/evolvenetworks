/********************************************************

	Copyright (C) 2009 Deepak Chandran
	see header file
	
	Build the cvode library first by compiling all the source files in cvode_src and 
	using ar *.o -o libcvode.a
	
	Use the following command to compile:
	
	gcc mtrand.c ssa.c cvodeim.c ga.c loops.c mass_action_network.c -lm -lcvode

********************************************************/
#include "loops.h"
#include "evolve_chemotaxis.h"

/**********************************************************************
  The following global parameters are used to simulate chemotaxis
**********************************************************************/

#define PI 3.14159265
#define DOUBLEPI 6.283185
static int STOCHASTIC;
static double MAX, X0, Y0, SPEED, GRIDSZ, SPREAD, TIME;
static FILE * FOUT;
static int FOOD, ANGLE, THRUST;
static double MEASURING_TIME;

/**********************************************************************/

//Gaussian function defining attractor
double gaussian(double x0, double y0, double x, double y)
{
	double distsq = (x0-x)*(x0-x) + (y0-y)*(y0-y);
	return exp(- distsq/SPREAD );
}

//ODE function for simulating chemotaxis
void chemotaxisODE(double time,double* u,double* du,void * p)
{
	int i,j,h;
	double angle, * rates, * N;
	int numReactions, numVars;
	chemotaxis_network * net;
	
	net = (chemotaxis_network*)p;
	u[FOOD] = MAX*gaussian(X0,Y0,(*net).x,(*net).y);  //update food
	
	(*net).angle += /*(mtrand()-0.5) */ SPEED * u[ANGLE]; //update angle
	
	//make sure angle is between 0 and 2*PI
	if ((*net).angle > DOUBLEPI)
		(*net).angle -= DOUBLEPI;
		
	if ((*net).angle < DOUBLEPI)
		(*net).angle += DOUBLEPI;
	
	angle = (*net).angle;
	
	//update coordinates
	if (THRUST > 0)
	{
		(*net).x += SPEED * (u[THRUST]) * cos(angle);
		(*net).y += SPEED * (u[THRUST]) * sin(angle);
	}
	else
	{
		(*net).x += SPEED * cos(angle);
		(*net).y += SPEED * sin(angle);
	}
	
	//make sure location is within the grid
	if ((*net).x > GRIDSZ) (*net).x = -GRIDSZ;
	if ((*net).y > GRIDSZ) (*net).y = -GRIDSZ;
	if ((*net).x < -GRIDSZ) (*net).x = GRIDSZ;
	if ((*net).y < -GRIDSZ) (*net).y = GRIDSZ;
	
	//This is for recording the variability in the path, which is used to compute fitness
	if (time > MEASURING_TIME)
	{
		if ((*net).x > (*net).xmax || (*net).xmax == 0.0) (*net).xmax =(*net).x;
		if ((*net).x < (*net).xmin || (*net).xmin == 0.0) (*net).xmin =(*net).x;
		if ((*net).y > (*net).ymax || (*net).ymax == 0.0) (*net).ymax =(*net).y;
		if ((*net).y < (*net).ymin || (*net).ymin == 0.0) (*net).ymin =(*net).y;
	}
	
	//This is for printing the path
	if (FOUT)
		fprintf(FOUT, "%lf\t%lf\t%lf\t%lf\n",u[0],(*net).x,(*net).y,(*net).angle);
	
	//ODEfunction(time,u,du,(*net).network);  //update species
	
	//use rates and stoichiometry functions from the network evolution  library
	
	ReactionNetwork * rnet = (*net).network;
	
	numReactions = getNumReactions(rnet);
	numVars = getNumSpecies(rnet);
	
	rates = getReactionRates(rnet, u);
	N = getStoichiometryMatrix(rnet);
	
	for (i=0; i < numVars; ++i)
	{
		du[i] = 0;
		for (j=0; j < numReactions; ++j)
		{
			if (getValue(N,numReactions,i,j) != 0)
				du[i] += rates[j]*getValue(N,numReactions,i,j);
		}
	}
	
	free(rates);
	free(N);
   
    du[0] = 0.0;
}

//Propensity function for simulating chemotaxis stochastically
void chemotaxisPropensity(double time,double* u,double* v,void * p)
{
	int i;
	double angle;
	int numReactions;
	double * rates;
	chemotaxis_network * net;
	ReactionNetwork * rnet;

	net = (chemotaxis_network*)p;
	
	u[FOOD] = MAX*gaussian(X0,Y0,(*net).x,(*net).y);  //update food
	
	(*net).angle += /*(mtrand()-0.5) */ SPEED * u[ANGLE]; 
	if ((*net).angle > DOUBLEPI)
		(*net).angle -= DOUBLEPI;
		
	if ((*net).angle < DOUBLEPI)
		(*net).angle += DOUBLEPI;
	
	angle = (*net).angle;
	
	//update coordinates
	if (THRUST > 0)
	{
		(*net).x += SPEED * (u[THRUST]) * cos(angle);
		(*net).y += SPEED * (u[THRUST]) * sin(angle);
	}
	else
	{
		(*net).x += SPEED * cos(angle);
		(*net).y += SPEED * sin(angle);
	}
	
	if ((*net).x > GRIDSZ) (*net).x = -GRIDSZ;
	if ((*net).y > GRIDSZ) (*net).y = -GRIDSZ;
	if ((*net).x < -GRIDSZ) (*net).x = GRIDSZ;
	if ((*net).y < -GRIDSZ) (*net).y = GRIDSZ;
	
	//used by fitness function
	if (time > MEASURING_TIME)
	{
		if ((*net).x > (*net).xmax || (*net).xmax == 0.0) (*net).xmax =(*net).x;
		if ((*net).x < (*net).xmin || (*net).xmin == 0.0) (*net).xmin =(*net).x;
		if ((*net).y > (*net).ymax || (*net).ymax == 0.0) (*net).ymax =(*net).y;
		if ((*net).y < (*net).ymin || (*net).ymin == 0.0) (*net).ymin =(*net).y;
	}
	
	//print path to file
	if (FOUT)
		fprintf(FOUT, "%lf\t%lf\t%lf\t%lf\n",u[0],(*net).x,(*net).y,(*net).angle);
	
	//call propensity function from the network evolution library
	rnet = (*net).network;
	
	numReactions = getNumReactions(rnet);
	
	rates = getReactionRates(rnet, u);
	
	for (i=0; i < numReactions; ++i) v[i] = rates[i];
}

//compute fitness of a network
double fitness(void * p)
{
	int i,j,n,r,invalid, sz;
	double * N, * init, total, score, xinit, yinit, * y;
	chemotaxis_network cnet;
	ReactionNetwork * rnet = (ReactionNetwork*)p;  //make chemotaxis network
	
	n = getNumSpecies(rnet);
	r = getNumReactions(rnet);
	cnet.network = rnet;
	
	if (n < 4 || n > 10 || r > 12)
	{
		return 0.0; //minimum size: 1 species for input, 1 for turning
	}
	
	N = getStoichiometryMatrix(rnet);	
	init = malloc(n * sizeof(double));
	score = total = 0.0;
	
	for (j=0; j < 10; ++j)  //take average fitness from 10 different initial positions
	{
		cnet.angle = 0.0;
		//cnet.x = xinit = (2.0*mtrand()-1.0) * GRIDSZ;
		//cnet.y = yinit = (2.0*mtrand()-1.0) * GRIDSZ;
		
		cnet.x = xinit = cos(j*2*PI/10.0) * GRIDSZ;  //initial position is around the origin
		cnet.y = yinit = sin(j*2*PI/10.0) * GRIDSZ;
		MEASURING_TIME = 4.0 * TIME/5.0;
		cnet.xmax = cnet.xmin = cnet.ymax = cnet.ymin = 0.0;
		
		for (i=0; i < n; ++i)
			init[i] = 0.0;
		
		sz = TIME * 10;
		
		y = 0;
		//if (STOCHASTIC)
			//y = SSA( n, r, N, &(chemotaxisPropensity), init, 0, TIME, 100000, &sz, &cnet); //simulate stochastically
		//else
			y = ODEsim( n , init, &(chemotaxisODE), 0, TIME, 0.1, &cnet); //simulate deterministically
		score = 0.0;
		
		if (y && sz > TIME)
		{
			invalid = 1;
			
			if (getValue(y,1+n,sz-1,2) != 0.0)
			{
				invalid = 0;
			}
			
			if (invalid)
				score = 0.0;
			else
			{
				for (i=(int)(sz/2); i < sz; ++i)
				{
					score += getValue(y,1+n,i,1+FOOD);
				}
				score /= (int)(sz/2);
			}
			
			free(y);
		}
		
		total += score/(1.0 + sqrt((cnet.xmax - cnet.xmin)*(cnet.xmax - cnet.xmin) + (cnet.ymax - cnet.ymin)*(cnet.ymax - cnet.ymin)));
		//total += score;
	}
	free(N);
	free(init);
		
	return total/10.0;  //fitness
}

/* prints the fitness of the best network during each iteration of the genetic algorithm */
int callback(int iter,GApopulation pop,int popSz)
{
	ReactionNetwork * net = (ReactionNetwork*)pop[0];
	
	printf("%i\t%i\t%lf\n",iter,getNumSpecies(net),fitness(pop[0]));
	return 0;
}

int main(int args, char ** argv)
{
	char * outfile = "out.tab";
	double * N = 0, * J = 0, * init, *y, * p;
	int sz, initsz = 1000, iter = 30, i,j,n,r;
	ReactionNetwork * rnet;
	LoopsInformation info;
	
	//parameters used in simulating chemotaxis
	TIME = 500;
	FOUT = 0;
	MAX = 100.0;  //global settings
	X0 = 0.0;
	Y0 = 0.0;
	SPEED = 0.1;
	GRIDSZ = 20.0;
	SPREAD = 200.0;
	STOCHASTIC = 0;
	
	FOOD = 0;
	ANGLE = 1;
	THRUST = 2;
	
	/**** get parameters from command line arguments, if they exist ****/
	if (args > 1)
		STOCHASTIC = (int)(argv[1][0] == 's' || argv[1][0] == 'S');
	
	if (args > 2)
		iter = atoi(argv[2]);
	
	if (args > 3)
		initsz = atoi(argv[3]);
	
	int numspecs = 5, numreacs = 2;
	
	if (args > 5)
	{
		numspecs = atoi(argv[4]);
		numreacs = atoi(argv[5]);
	}
	
	if (args > 6)
	{
		MAX = atof(argv[6]);
		SPREAD = MAX*2.0;
	}

	if (args > 7)
		TIME = atof(argv[7]);
	
	if (args > 8)
		outfile = argv[8];
	
	//GApopulation pop = randomNetworks(initsz,numspecs,numreacs);
	
	/**** print the parameters being used ****/
	
	printf("Evolve chemotaxis using:\n\t");
	if (STOCHASTIC)
		printf("stochastic");
	else
		printf("deterministic");

	printf(" simulation\n\tgenerations = %i\n\tinitial population size = %i\n\tinitial num. species = %i\n\tinitial num. reactions = %i\n\tmaximum fitness = %lf\n\tsimulation time = %lf\n\toutput file = %s\n",iter,initsz,numspecs,numreacs,MAX,TIME,outfile);
	
	
	
	/*****************************************************************************************
	initialize and run the genetic algorithm with functions from the network evolution library 
	******************************************************************************************/
	
	//GAinit(&deleteNetwork, &cloneNetwork ,&fitness, &crossover, &mutate, 0);
	
	setNetworkType(MASS_ACTION_NETWORK);
	
	setFitnessFunction(&fitness);
	
	setInitialNetworkSize(numspecs, numreacs);
	
	GApopulation pop = evolveNetworks(initsz,initsz/2,iter, &callback);
	
	//pop = GArun(pop,initsz,initsz/2,iter, &callback);
	
	
	/*****************************************************************************************/
	
	//print the chemotaxis path
	
	FOUT = fopen(outfile,"w");
	
	chemotaxis_network cnet;
	rnet = (ReactionNetwork*)(pop[0]);
	
	n = getNumSpecies(rnet);
	r = getNumReactions(rnet);
	
	printNetwork(rnet);
	
	cnet.network = rnet;	
	
	init = malloc(n * sizeof(double));
	
	cnet.angle = 0.0;
	cnet.x = GRIDSZ;
	cnet.y = GRIDSZ;
	for (i=0; i < n; ++i)
		init[i] = 0.0;
	
	sz = TIME * 10;
	y = 0;
	
	if (STOCHASTIC)	
	{
		N = getStoichiometryMatrix(rnet);
		y = SSA( n, r, N, &(chemotaxisPropensity), init, 0, TIME, 100000, &sz, &cnet);	
		free(N);
	}
	else
	{
		y = ODEsim( n , init, &(chemotaxisODE), 0, TIME, 0.1, &cnet);
	}
	
	if (!y)
	{
		for (i=0; i < (initsz/2); ++i)
			deleteNetwork(pop[i]);
		
		return 0;
	}
	
	/*********   find feedback loops   ***********/
	
	p = malloc( n * sizeof(double) );
	
	for (i=0; i < n; ++i)
		p[i] = getValue(y,(1+n),sz-1,i+1);
		
	J = jacobian(n, p, &(chemotaxisODE), &cnet);
	
	info = getLoops(J,n);
	
	for (i=0; i < info.numLoops; ++i)
	{
		printf("length = %i\ttype = %i\t",info.loopLengths[i],info.loopTypes[i]);
		for (j=0; j < info.loopLengths[i]; ++j)
		{
			printf("%i\t",info.nodes[i][j]);
		}
		printf("\n");
	}

	free(p);
	freeLoopsInfo(info);
	
	/********* free memory  *********/	
	
	free(init);
	
	free(y);
	
	fclose(FOUT);
	FOUT = 0;
	
	for (i=0; i < (initsz/2); ++i)
		deleteNetwork(pop[i]);
	
	return 0;
}
