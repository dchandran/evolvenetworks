/*******************************************************

	Copyright (C) 2009 Deepak Chandran
	see header file

********************************************************/

#include "enzymeNetwork.h"

/************************
    global variables
*************************/

static double MUTATE_KM_PROB = 0.2;
static double KM_RANGE = 20.0;
static double MUTATE_ENZYME_PROB = 0.2;
static double CROSSOVER_PROB = 0.2;

void setParametersForEnzymeNetwork(double uni_uni, double uni_bi, double bi_uni, double bi_bi, double no_reactant, double no_product, double avg_rate_constant, double avg_km)
{
	if (avg_km > 0)
		KM_RANGE = avg_km;
	setParametersForMassActionNetwork(uni_uni, uni_bi, bi_uni, bi_bi, no_reactant, no_product, avg_rate_constant);
}

void setSizeForEnzymeNetwork(int n, int s)
{
	setSizeForMassActionNetwork(n,s);
}

void setMutationAndCrossoverRatesForEnzymeNetwork(double mutateE, double mutatek, double remove, double add, double crossover)
{
	double total;

	if (mutatek < 0) mutatek = 0;
	if (remove < 0) remove = 0;
	if (add < 0) add = 0;
	if (mutateE < 0) mutateE = 0;
		
	total = mutatek + remove + add + mutateE;
	if (total == 0) 
		mutatek = remove = mutateE = 0.33333;
	else
	{
		mutateE /= total;
		mutatek /= total;
		remove /= total;
	}
	setMutationAndCrossoverRatesForMassActionNetwork(mutatek/2.0,remove,add,1.0);

	MUTATE_KM_PROB = mutatek/2.0;
	MUTATE_ENZYME_PROB = mutateE;
	CROSSOVER_PROB = crossover;
}

/********************************************************
    Clone, delete, mutate, crossover  (required by GA)
*********************************************************/

void deleteEnzymeCatalyzedNetwork(GAindividual individual)
{
	EnzymeCatalyzedNetwork * net;
	
	if (!individual) return;
	
	net = (EnzymeCatalyzedNetwork*)(individual);
	
	if (net->massActionNetwork)
	{
		deleteMassActionNetwork(net->massActionNetwork);
	}

	if (net->enzymes)
		free (net->enzymes);
	

	if (net->Km)
		free(net->Km);
}

GAindividual cloneEnzymeCatalyzedNetwork(GAindividual individual)
{
	int i,n;
	EnzymeCatalyzedNetwork * net, * net2;
	
	if (!individual) return 0;
	
	net = (EnzymeCatalyzedNetwork*)(individual);   //original
	net2 = (EnzymeCatalyzedNetwork*) malloc(sizeof(EnzymeCatalyzedNetwork)); //soon to be clone
	
	net2->massActionNetwork = (MassActionNetwork*) cloneMassActionNetwork(net->massActionNetwork);

	n = net->massActionNetwork->reactions;
	net2->enzymes = (int*) malloc (n * sizeof(int));
	net2->Km = (double*) malloc (n * sizeof(double));

	for (i=0; i < n; ++i)
	{
		net2->enzymes[i] = net->enzymes[i];
		net2->Km[i] = net->Km[i];
	}

	return (GAindividual)(net2);  //done
}

GAindividual crossoverEnzymeCatalyzedNetwork(GAindividual individualA, GAindividual individualB)  //crossover between complexes in two networks
{
	int i,j;
	EnzymeCatalyzedNetwork * net1, * net2, * net3;
	
	if (mtrand() > CROSSOVER_PROB) return mutateEnzymeCatalyzedNetwork(cloneEnzymeCatalyzedNetwork(individualA));
	
	if (!individualA) return mutateEnzymeCatalyzedNetwork(cloneEnzymeCatalyzedNetwork(individualB));
	if (!individualB) return mutateEnzymeCatalyzedNetwork(cloneEnzymeCatalyzedNetwork(individualA));
	
	net1 = (EnzymeCatalyzedNetwork*)(individualA);  //parents
	net2 = (EnzymeCatalyzedNetwork*)(individualB);
	
	if (net1->massActionNetwork->reactions < 3) return mutateEnzymeCatalyzedNetwork(cloneEnzymeCatalyzedNetwork(net2));  //if parents are too small
	if (net2->massActionNetwork->reactions < 3) return mutateEnzymeCatalyzedNetwork(cloneEnzymeCatalyzedNetwork(net1));
	
	net3 =  (EnzymeCatalyzedNetwork*) malloc (sizeof(EnzymeCatalyzedNetwork));

	net3->massActionNetwork = crossoverMassActionNetwork(net1->massActionNetwork,net2->massActionNetwork);

	net3->enzymes = (int*)malloc(net3->massActionNetwork->reactions * sizeof(int));

	net3->Km = (double*)malloc(net3->massActionNetwork->reactions * sizeof(double));

	//get some set of enzymes from one parent and some from the other

	j = (int)(mtrand() * net3->massActionNetwork->reactions);
	for (i=0; i < j && i < net1->massActionNetwork->reactions && i < net3->massActionNetwork->reactions; ++i)
	{
		net3->enzymes[i] = net1->enzymes[i];
		net3->Km[i] = net1->Km[i];
		if (net3->enzymes[i] >= net3->massActionNetwork->species)
			net3->enzymes[i] = (int)(mtrand() * net3->massActionNetwork->species);
	}

	for (i=0; (i+j) < net3->massActionNetwork->reactions && i < net2->massActionNetwork->reactions; ++i)
	{
		net3->enzymes[i+j] = net2->enzymes[i];
		net3->Km[i+j] = net2->Km[i];
		if (net3->enzymes[i+j] >= net3->massActionNetwork->species)
			net3->enzymes[i+j] = (int)(mtrand() * net3->massActionNetwork->species);
	}

	//still leftover?
	j=i;
	for (i=j; i < net3->massActionNetwork->reactions; ++i)
	{
		net3->enzymes[i] = (int)(mtrand() * net3->massActionNetwork->species);
		net3->Km[i] = KM_RANGE * mtrand();
	}
	
	return (GAindividual)(net3);
}

GAindividual mutateEnzymeCatalyzedNetwork(GAindividual individual)
{
	int i;
	double r;
	EnzymeCatalyzedNetwork * net;

	if (!individual) return individual;
	
	net = (EnzymeCatalyzedNetwork*)individual;
	
	r = mtrand();

	if (r < MUTATE_KM_PROB) //mutate km
	{
		i = (int)(mtrand() * net->massActionNetwork->reactions);
		net->Km[i] *= 2.0 * mtrand();
	}
	else
	if (r < (MUTATE_KM_PROB+MUTATE_ENZYME_PROB))   //mutate enzyme
	{
		i = (int)(mtrand() * net->massActionNetwork->reactions);
		net->enzymes[i] = (int)(mtrand() * net->massActionNetwork->species);
		return (GAindividual)(net);
	}
	
	net->massActionNetwork = mutateMassActionNetwork(net->massActionNetwork);
	
    return (GAindividual)(net);
}

/*****************************************************
   Functions for simulating and printing
******************************************************/

int getNumSpeciesForEnzymeCatalyzedNetwork(GAindividual individual)
{
	EnzymeCatalyzedNetwork * net = (EnzymeCatalyzedNetwork*)(individual);
	if (!net) return 0;
	return (getNumSpeciesForMassActionNetwork(net->massActionNetwork));
}

int getNumReactionsForEnzymeCatalyzedNetwork(GAindividual individual)
{
	EnzymeCatalyzedNetwork * net = (EnzymeCatalyzedNetwork*)(individual);
	if (!net) return 0;
	return (getNumReactionsForMassActionNetwork(net->massActionNetwork));
}

void setFixedSpeciesForEnzymeCatalyzedNetwork(GAindividual individual, int i, int value)
{
	EnzymeCatalyzedNetwork * net = (EnzymeCatalyzedNetwork*)(individual);
	if (!net) return;
	setFixedSpeciesForMassActionNetwork(net->massActionNetwork,i,value);
}

void ratesForEnzymeCatalyzedNetwork(double time,double* u,double* rate,GAindividual individual)
{
	int i,k;
	EnzymeCatalyzedNetwork * enet;
	MassActionNetwork * net;
	
	enet = (EnzymeCatalyzedNetwork*)(individual);
	
	if (!enet) return;

	net = enet->massActionNetwork;

	if (!net) return;
	
	for (i=0; i < net->reactions; ++i)
	{
		k =(((net->r1[i] == -1 && net->r2[i] > -1) ||  //uni-uni
			(net->r1[i] > -1 && net->r2[i] == -1))
			&&
			((net->p1[i] == -1 && net->p2[i] > -1) ||
			(net->p1[i] > -1 && net->p2[i] == -1)));
		
		rate[i] = net->k[i];
		if (k)
			rate[i] *= u[ enet->enzymes[i] ];

		if (net->r1[i] > -1)
		{
			rate[i] *= u[ net->r1[i] ];
			if (k)
				rate[i] /= enet->Km[i] + u[ net->r1[i] ];
		}
		if (net->r2[i] > -1) 
		{
			rate[i] *= u[ net->r2[i] ];
			if (k)
				rate[i] /= enet->Km[i] + u[net->r2[i] ];
		}

		
	}
}

double * stoichiometryForEnzymeCatalyzedNetwork(GAindividual individual)
{
	EnzymeCatalyzedNetwork * net = (EnzymeCatalyzedNetwork*)(individual);
	if (!net) return 0;
	return (stoichiometryForMassActionNetwork(net->massActionNetwork));
}

void printEnzymeCatalyzedNetwork(FILE * stream,GAindividual individual)
{
	int i,fix;
	MassActionNetwork * net;
	EnzymeCatalyzedNetwork * enet;
	
	if (!individual) return;
	enet = (EnzymeCatalyzedNetwork*)individual;
	net = enet->massActionNetwork;

	if (!net) return;
	
	fix = 0;
	for (i=0; i < net->species; ++i)
	{
		if (net->fixed[i])
		{
			fix = i+1;
			break;
		}
	}

	if (fix)
	{
		fprintf(stream, "const s%i",fix);
		for (i=0; i < net->species; ++i)
		{
			if (net->fixed[i])			
				fprintf(stream, ", s%i",i+1);			
		}
		fprintf(stream, "\n");
	}
	
	for (i=0; i < net->reactions; ++i)
	{
		if (net->r1[i] > -1 && net->r2[i] > -1)
			fprintf(stream, "s%i + s%i -> ",net->r1[i]+1,net->r2[i]+1);
		else
		{
			if (net->r1[i] > -1 || net->r2[i] > -1)
			{
				if (net->r1[i] > -1)
					fprintf(stream, "s%i ->",net->r1[i]+1);
				else
					fprintf(stream, "s%i ->",net->r2[i]+1);
			}
			else
			{
				fprintf(stream, "$x -> ");
			}
		}

		if (net->p1[i] > -1 && net->p2[i] > -1)
			fprintf(stream, "s%i + s%i; ",net->p1[i]+1,net->p2[i]+1);
		else
		{
			if (net->p1[i] > -1 || net->p2[i] > -1)
			{
				if (net->p1[i] > -1)
					fprintf(stream, "s%i; ",net->p1[i]+1);
				else
					fprintf(stream, "s%i; ",net->p2[i]+1);
			}
			else
			{
				fprintf(stream, "$x; ");
			}
		}
		
		//rate
		if (net->r1[i] > -1 && net->r2[i] > -1)
			fprintf(stream, "k%i * s%i * s%i ",i+1,net->r1[i]+1,net->r2[i]+1);
		else
		{
			if (net->r1[i] > -1 || net->r2[i] > -1)
			{
				if (net->r1[i] > -1)
				{
					if ((net->p1[i] == -1 && net->p2[i] > -1) ||  //uni-uni
						(net->p1[i] > -1 && net->p2[i] == -1))
					{
						fprintf(stream, "k%i * s%i * s%i/ (km%i + s%i)",i+1,enet->enzymes[i]+1,net->r1[i]+1,i+1,net->r1[i]);
					}
					else
						fprintf(stream, "k%i * s%i",i+1,net->r1[i]+1);
				}
				else
				{
					if ((net->p1[i] == -1 && net->p2[i] > -1) ||  //uni-uni
						(net->p1[i] > -1 && net->p2[i] == -1))
					{
						fprintf(stream, "k%i * s%i * s%i/ (km%i + s%i)",i+1,enet->enzymes[i]+1,net->r2[i]+1,i+1,net->r2[i]);
					}
					else
						fprintf(stream, "k%i * s%i",i+1,net->r2[i]+1);
				}
			}
			else
			{
				fprintf(stream, "k%i",i+1);
			}
		}
		fprintf(stream, ";\n");
	}

	fprintf(stream, "\n");
	
	for (i=0; i < net->reactions; ++i)
	{
		fprintf(stream, "k%i = %lf;\n",i+1,net->k[i]);
		if (((net->r1[i] == -1 && net->r2[i] > -1) ||  //uni-uni
			(net->r1[i] > -1 && net->r2[i] == -1))
			&&
			((net->p1[i] == -1 && net->p2[i] > -1) ||
			(net->p1[i] > -1 && net->p2[i] == -1)))
		{
			fprintf(stream, "km%i = %lf;\n",i+1,enet->Km[i]);
		}
	}
}

/***********************
  GA related functions
***********************/

GApopulation randomEnzymeCatalyzedNetworks(int num)
{
	int i,j;
	MassActionNetwork * mnet;
	EnzymeCatalyzedNetwork * enet;
	GApopulation pop = randomMassActionNetworks(num);

	for (i=0; i < num; ++i)
	{
		mnet = (MassActionNetwork*)pop[i];
		enet = (EnzymeCatalyzedNetwork*) malloc(sizeof(EnzymeCatalyzedNetwork));
		enet->enzymes = (int*) malloc( mnet->reactions * sizeof(int) );
		enet->Km = (double*) malloc( mnet->reactions * sizeof(double) );
		for (j=0; j < mnet->reactions; ++j)
		{
			enet->enzymes[j] = (int)(mtrand() * mnet->species);
			enet->Km[j] = KM_RANGE * mtrand();
		}
		enet->massActionNetwork = mnet;
		pop[i] = enet;
	}
	
	return pop;
}

EnzymeCatalyzedNetwork * newEnzymeCatalyzedNetwork(int m,int n)
{
	int i;
	EnzymeCatalyzedNetwork * net;

	net = (EnzymeCatalyzedNetwork*) malloc(sizeof(EnzymeCatalyzedNetwork));

	net->massActionNetwork = newMassActionNetwork(m,n);

	net->enzymes = (int*) malloc( n * sizeof(int) );

	net->Km = (double*) malloc( n * sizeof(double) );

	for (i=0; i < n; ++i)
	{
		net->enzymes[i] = 0;
		net->Km[i] = 0.0;
	}
	
	return net;
}
