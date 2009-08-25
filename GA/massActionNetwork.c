/*******************************************************

	Copyright (C) 2009 Deepak Chandran
	see header file

********************************************************/

#include "massActionNetwork.h"

/************************
    global variables
*************************/

static double CROSSOVER_PROB = 1.0;
static double MUTATE_COEFF_PROB = 0.5;
static double PROB_SECOND_REACTANT = 0.3;
static double PROB_SECOND_PRODUCT = 0.3;
static double AVG_RATE_CONSTANT = 1.0;
static int AVG_NUM_SPECIES = 1.0;
static int AVG_NUM_REACTIONS = 1.0;

void setParametersForMassActionNetwork(double prob1, double prob2, double k)
{
	AVG_RATE_CONSTANT = k;
	PROB_SECOND_REACTANT = prob1;
	PROB_SECOND_PRODUCT = prob2;
}

void setSizeForMassActionNetwork(int n, int s)
{
	AVG_NUM_SPECIES = n;
	AVG_NUM_REACTIONS = s;
}

void setMutationAndCrossoverRatesForMassActionNetwork(double a, double b)
{
	if (a > 1.0) a /= 100.0;
	if (b > 1.0) b /= 100.0;
	MUTATE_COEFF_PROB = a;
	CROSSOVER_PROB = b;
}

/*******************************************************
    Clone, delete, mutate, crossover  (required by GA)
*********************************************************/

void deleteMassActionNetwork(void * individual)
{
	MassActionNetwork * net;
	
	if (!individual) return;
	net = (MassActionNetwork*)(individual);
	
	if ((*net).r1) free((*net).r1);  //free each array
	if ((*net).r2) free((*net).r2);
	if ((*net).p1) free((*net).p1);
	if ((*net).p2) free((*net).p2);
	if ((*net).k) free((*net).k);
}

void* cloneMassActionNetwork(void * individual)
{
	int i,m,n;
	MassActionNetwork * net, * net2;
	
	if (!individual) return 0;
	
	net = (MassActionNetwork*)(individual);   //original
	net2 = malloc(sizeof(MassActionNetwork)); //soon to be clone
	
	m = (*net).reactions;    //number of reactions
	n = (*net).species;    //number of species
	(*net2).reactions = m;
	(*net2).species = n;
	
	(*net2).k = malloc(m * sizeof(double));   //allocate space
	(*net2).r1 = malloc(m * sizeof(int));
	(*net2).r2 = malloc(m * sizeof(int));
	(*net2).p1 = malloc(m * sizeof(int));
	(*net2).p2 = malloc(m * sizeof(int));
	
	for (i=0; i < m; ++i)   //copy values
	{
		(*net2).k[i] = (*net).k[i];
		(*net2).r1[i] = (*net).r1[i];
		(*net2).r2[i] = (*net).r2[i];
		(*net2).p1[i] = (*net).p1[i]; 
		(*net2).p2[i] = (*net).p2[i];
	}
	
	return (void*)(net2);  //done
}

void* crossoverMassActionNetwork(void * individualA, void * individualB)
{
	int i, j, i1, i2, n = 0, m = 0; 
	MassActionNetwork * net1, * net2, * net3;
	
	if (mtrand() > CROSSOVER_PROB) return mutateMassActionNetwork(cloneMassActionNetwork(individualA));  //do crossover?
	
	if (!individualA) return mutateMassActionNetwork(cloneMassActionNetwork(individualB));
	if (!individualB) return mutateMassActionNetwork(cloneMassActionNetwork(individualA));
	
	net1 = (MassActionNetwork*)(individualA);  //parents
	net2 = (MassActionNetwork*)(individualB);
	
	if ((*net1).reactions < 3) return mutateMassActionNetwork(cloneMassActionNetwork(net2));  //if parents are too small
	if ((*net2).reactions < 3) return mutateMassActionNetwork(cloneMassActionNetwork(net1));
	
	i1 = (int)(mtrand() * ((*net1).reactions - 1) + 1.0);	//crossover point in net1
	i2 = (int)(mtrand() * ((*net2).reactions - 2) + 1.0);	//crossover point in net2
	
	n = i1 + (*net2).reactions - i2;
	
	net3 = newMassActionNetwork(m,n);  //child network
	
	for (i=0; i < i1; ++i)
	{
		(*net3).k[i] = (*net1).k[i];
		(*net3).r1[i] = (*net1).r1[i];
		(*net3).r2[i] = (*net1).r2[i];
		(*net3).p1[i] = (*net1).p1[i]; 
		(*net3).p2[i] = (*net1).p2[i];
		
		if ((*net1).r1[i] > m) m = (*net1).r1[i];
		if ((*net1).r2[i] > m) m = (*net1).r2[i];
		if ((*net1).p1[i] > m) m = (*net1).p1[i];
		if ((*net1).p2[i] > m) m = (*net1).p2[i];
	}
	
	for (i=i2; i < (*net2).reactions; ++i)
	{
		(*net3).k[i+i1-i2] = (*net2).k[i];
		(*net3).r1[i+i1-i2] = (*net2).r1[i];
		(*net3).r2[i+i1-i2] = (*net2).r2[i];
		(*net3).p1[i+i1-i2] = (*net2).p1[i]; 
		(*net3).p2[i+i1-i2] = (*net2).p2[i];
		
		if ((*net2).r1[i] > m) m = (*net2).r1[i];
		if ((*net2).r2[i] > m) m = (*net2).r2[i];
		if ((*net2).p1[i] > m) m = (*net2).p1[i];
		if ((*net2).p2[i] > m) m = (*net2).p2[i];
	}
	
	(*net3).species = m + 1;
	return (void*)(net3);
}

void* mutateMassActionNetwork(void * individual)
{
	int i,j,j2,k,m,n;
	MassActionNetwork * net, * net2;
	
	net = (MassActionNetwork*)individual;

	m = (*net).reactions;
	n = (*net).species;

	i = (int)(mtrand() * m);  //pick random reaction

	if (mtrand() < MUTATE_COEFF_PROB)   //mutate coefficient
	{
		(*net).k[i] *= (mtrand() * 2.0);
		return (void*)(net);
	}
	else              //add or remove a new reaction and/or species to the network
	{
		if (mtrand() < 0.5 && m > 2)     //remove a reaction
		{
			net2 = newMassActionNetwork( n, (m - 1) );
			for (j=0,j2=0; j < m; ++j)
			{
				if (j != i)
				{
					(*net2).k[j2] = (*net).k[j];
					(*net2).r1[j2] = (*net).r1[j];
					(*net2).r2[j2] = (*net).r2[j];
					(*net2).p1[j2] = (*net).p1[j];
					(*net2).p2[j2] = (*net).p2[j];
					++j2;
				}
			}
			deleteMassActionNetwork(net);
			return (void*)(net2);
		}
		else
		{
			net2 = newMassActionNetwork( n, (m + 1) );
			for (j=0; j < m; ++j)
			{
				(*net2).k[j] = (*net).k[j];
				(*net2).r1[j] = (*net).r1[j];
				(*net2).r2[j] = (*net).r2[j];
				(*net2).p1[j] = (*net).p1[j];
				(*net2).p2[j] = (*net).p2[j];
			}
			(*net2).k[j] = 2.0 * mtrand();   //reaction rate constant
			
			if (mtrand() < 0.3) 
			{
				(*net2).r1[j] = -1;
				(*net2).r2[j] = (int)( (1 + (*net2).species) * mtrand());  //if r1 is empty, r2 must not be empty
				if ((*net2).r2[j] >= (*net2).species)
					(*net2).species = (*net2).r2[j] + 1;
			}
			else
			{
				(*net2).r1[j] = (int)( (1 + (*net2).species) * mtrand());  //if r1 is not empty, r2 can be empty
				if (mtrand() < 0.3) 
					(*net2).r2[j] = -1;
				else
					(*net2).r2[j] = (int)( (1 + (*net2).species) * mtrand());
				
				if ((*net2).r1[j] >= (*net2).species)
					(*net2).species = (*net2).r1[j] + 1;
					
				if ((*net2).r2[j] >= (*net2).species)
					(*net2).species = (*net2).r2[j] + 1;
			}
			
			if (mtrand() < 0.3) 
			{
				(*net2).p1[j] = -1;
				(*net2).p2[j] = (int)( (1 + (*net2).species) * mtrand());
				if ((*net2).p2[j] >= (*net2).species)
					(*net2).species = (*net2).p2[j] + 1;
			}
			else
			{
				(*net2).p1[j] = (int)( (1 + (*net2).species) * mtrand());
				if (mtrand() < 0.3) 
					(*net2).p2[j] = -1;
				else
					(*net2).p2[j] = (int)( (1 + (*net2).species) * mtrand());
				
				if ((*net2).p1[j] >= (*net2).species)
					(*net2).species = (*net2).p1[j] + 1;
				
				if ((*net2).p2[j] >= (*net2).species)
					(*net2).species = (*net2).p2[j] + 1;
			}
			
			deleteMassActionNetwork(net);
			return (void*)(net2);
		}
	}
    return (net);
}

/*****************************************************
   Functions for simulating and printing
******************************************************/


void ratesForMassActionNetwork(double time,double* u,double* rate,void * p)
{
	int i;
	MassActionNetwork * net;
	
	net = (MassActionNetwork*)(p);
	
	for (i=0; i < (*net).reactions; ++i)
	{
		rate[i] = (*net).k[i];
		if ((*net).r1[i] > -1) rate[i] *= u[ (*net).r1[i] ];
		if ((*net).r2[i] > -1) rate[i] *= u[ (*net).r2[i] ];
	}
}


double * stoichiometryForMassActionNetwork(void * p)
{
	int i,j,n;
	double * N;
	MassActionNetwork * net;
	
	net = (MassActionNetwork*)(p);
	n = (*net).reactions;
	N = malloc((*net).species * n * sizeof(double));
	for (i=0; i < n; ++i)
	{
		for (j=0; j < (*net).species; ++j)
			getValue(N,n,j,i) = 0;
		
		if ((*net).r1[i] >= 0) getValue(N,n,(*net).r1[i],i) = -1;
		if ((*net).r2[i] >= 0) getValue(N,n,(*net).r2[i],i) = -1;
		if ((*net).p1[i] >= 0) getValue(N,n,(*net).p1[i],i) = 1;
		if ((*net).p2[i] >= 0) getValue(N,n,(*net).p2[i],i) = 1;
	}
	return N;
}

void printMassActionNetwork(void * individual)
{
	int i,j,k;
	MassActionNetwork * net;
	
	if (!individual) return;
	net = (MassActionNetwork*)individual;
	
	for (i=0; i < (*net).reactions; ++i)
	{
		printf("%i + %i -> %i + %i;\t%lf\n",(*net).r1[i]+1,(*net).r2[i]+1,(*net).p1[i]+1,(*net).p2[i]+1,(*net).k[i]);
	}
}

/***********************
  GA related functions
***********************/

GApopulation randomMassActionNetworks(int num)
{
	int s = AVG_NUM_SPECIES;
	int r = AVG_NUM_REACTIONS;
	int i,j,n;
	MassActionNetwork * net;
	MassActionNetwork ** array;
	
	initMTrand(); /*initialize seeds for MT random number generator*/
	array = malloc(num * sizeof(MassActionNetwork*));
	for (i=0; i < num; ++i)
	{
		n = (int)(1 + s * 2.0 * mtrand());
		net = newMassActionNetwork(n,(int)(2 + r * 2.0 * mtrand()));
		
		for (j=0; j < (*net).reactions; ++j)
		{
			(*net).k[j] = AVG_RATE_CONSTANT * 2.0 * mtrand();   //reaction rate constant
			
			(*net).r1[j] = (int)((*net).species * mtrand());  //first reactant
			
			(*net).p1[j] = (int)((*net).species * mtrand());  //first product
			
			if (mtrand() < PROB_SECOND_REACTANT) 
				(*net).r2[j] = (int)((*net).species * mtrand()); //second reactant
				
			if (mtrand() < PROB_SECOND_PRODUCT) 
				(*net).p2[j] = (int)((*net).species * mtrand()); //second product
		}
		
		array[i] = net;
	}
	
	return (GApopulation)(array);
}

MassActionNetwork * newMassActionNetwork(int s,int r)
{
	MassActionNetwork * net;

	net = malloc(sizeof(MassActionNetwork));
	(*net).species = s;
	(*net).reactions = r;
		
	(*net).k = malloc((*net).reactions * sizeof(double));
	(*net).r1 = malloc((*net).reactions * sizeof(int));
	(*net).r2 = malloc((*net).reactions * sizeof(int));
	(*net).p1 = malloc((*net).reactions * sizeof(int));
	(*net).p2 = malloc((*net).reactions * sizeof(int));
	
	return net;
}

/*
int main()  //just for testing
{
	int i;
	GApopulation pop = randomMassActionNetworks(5,2);
	
	for (i=0; i < 2; ++i)
	{
		printMassActionNetwork(pop[i]);
		printf("\n");
	}
	
	for (i=0; i < 100; ++i)
	{
		void * net = crossoverMassActionNetwork(pop[0],pop[1]);
		printMassActionNetwork(net);
		net = mutateMassActionNetwork(net);
		printf("\n");
		deleteMassActionNetwork(net);
	}
	
	for (i=0; i < 2; ++i)
	{
		deleteMassActionNetwork(pop[i]);
	}
}
*/
