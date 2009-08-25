#include "reactionNetwork.h"

/********************************************************

		Pointers to the existing functions

*********************************************************/

static int NUMBER_OF_NETWORK_TYPES = 3;

static double networkProbs[] =
{
	0.33333333333,
	0.33333333333,
	0.33333333333
};

static PropensityFunction rateFunctions[] =
{
	&ratesForMassActionNetwork,
	&ratesForProteinInteractionNetwork, 
	&ratesForGeneRegulationNetwork
};

static GACrossoverFnc crossoverFunctions[] =
{
	&crossoverMassActionNetwork, 
	&crossoverProteinInteractionNetwork, 
	&crossoverGeneRegulationNetwork	
};

static GAMutateFnc mutateFunctions[] =
{
	&mutateMassActionNetwork, 
	&mutateProteinInteractionNetwork,
	&mutateGeneRegulationNetwork	
};

static GADeleteFnc deleteFunctions[] =
{
	&deleteMassActionNetwork,
	&deleteProteinInteractionNetwork,
	&deleteGeneRegulationNetwork
};

static GACloneFnc cloneFunctions[] =
{
	&cloneMassActionNetwork,
	&cloneProteinInteractionNetwork,
	&cloneGeneRegulationNetwork
};

typedef double* (*StoichiometryFunction)(GAindividual);

static StoichiometryFunction stoicFunctions[] =
{
	&stoichiometryForMassActionNetwork,
	&stoichiometryForProteinInteractionNetwork, 
	&stoichiometryForGeneRegulationNetwork
};

typedef int (*GetNumSpeciesFunction)(GAindividual);

static GetNumSpeciesFunction getNumSpeciesFunctions[] =
{
	&getNumSpeciesForMassActionNetwork,
	&getNumSpeciesForProteinInteractionNetwork, 
	&getNumSpeciesForGeneRegulationNetwork
};

typedef int (*GetNumReactionsFunction)(GAindividual);

static GetNumReactionsFunction getNumReactionsFunctions[] =
{
	&getNumReactionsForMassActionNetwork,
	&getNumReactionsForProteinInteractionNetwork, 
	&getNumReactionsForGeneRegulationNetwork
};

typedef void (*PrintNetworkFunction)(GAindividual);

static PrintNetworkFunction printNetworkFunctions[] =
{
	&printMassActionNetwork,
	&printProteinInteractionNetwork, 
	&printGeneRegulationNetwork
};

typedef void (*SetFixedSpeciesFunction)(GAindividual,int,int);

static SetFixedSpeciesFunction setFixedSpeciesFunctions[] =
{
	&setFixedSpeciesForMassActionNetwork,
	&setFixedSpeciesForProteinInteractionNetwork, 
	&setFixedSpeciesForGeneRegulationNetwork
};

/***********************************************************/

GApopulation randomNetworks(int sz0)
{
	int i, total = 0, k = 0;
	ReactionNetwork * rnet;
	double r;
	GApopulation P1 = 0, P2 = 0, P3 = 0, P;

	P = malloc( sz0 * sizeof (GAindividual) );

	r = networkProbs[MASS_ACTION_NETWORK] * sz0;

	if (r > 0 && r < sz0)
	{
		P1 = randomMassActionNetworks((int)r);
		for (i=0; i < r; ++i)
		{
			rnet = malloc(sizeof(ReactionNetwork));
			rnet->type = MASS_ACTION_NETWORK;
			rnet->network = P1[i];

			P[i] = rnet;
		}
		total += (int)r;
	}

	r = networkProbs[PROTEIN_INTERACTION_NETWORK] * sz0;
	if (r > 0 && r < sz0)
	{
		P2 = randomProteinInteractionNetworks((int)r);
		for (i=0; i < r; ++i)
		{
			rnet = malloc(sizeof(ReactionNetwork));
			rnet->type = PROTEIN_INTERACTION_NETWORK;
			rnet->network = P2[i];

			P[i] = rnet;
		}
		total += (int)r;
	}

	if (total < sz0)
	{
		k = sz0 - total;
		P3 = randomGeneRegulationNetworks(k);
		for (i=0; i < k; ++i)
		{
			rnet = malloc(sizeof(ReactionNetwork));
			rnet->type = GENE_REGULATION_NETWORK;
			rnet->network = P3[i];

			P[i] = rnet;
		}
	}

	return P;
}

void setCrossoverFunction( int i, GACrossoverFnc f)
{
	if (i < NUMBER_OF_NETWORK_TYPES)
	crossoverFunctions[i] = f;
}

void setMutationFunction( int i, GAMutateFnc f)
{
	if (i < NUMBER_OF_NETWORK_TYPES)
	mutateFunctions[i] = f;
}

void setRatesFunction( int i, PropensityFunction f)
{
	if (i < NUMBER_OF_NETWORK_TYPES)
	rateFunctions[i] = f;
}

void setStoichiometryFunction( int i, double* (*f)(GAindividual) )
{
	if (i < NUMBER_OF_NETWORK_TYPES)
	stoicFunctions[i] = f;
}

void setNetworkType(int p)
{
	int i;
	for (i=0; i < NUMBER_OF_NETWORK_TYPES; ++i) networkProbs[i] = 0.0;
	setNetworkTypeProbability(p, 1.0 );
}

void setNetworkTypeProbability(int i, double p)
{
	double total;
	
	if (i < NUMBER_OF_NETWORK_TYPES)
		networkProbs[ i ] = p;

	total = 0.0;

	for (i=0; i < NUMBER_OF_NETWORK_TYPES; ++i) total += networkProbs[i];

	if (total <= 0.0)
		for (i=0; i < NUMBER_OF_NETWORK_TYPES; ++i) networkProbs[i] = 0.33333;
	else
		for (i=0; i < NUMBER_OF_NETWORK_TYPES; ++i) networkProbs[i] /= total;
}

double * simulateNetworkODE( ReactionNetwork * r, double* iv, double time, double dt)
{	
	StoichiometryFunction stoic;
	PropensityFunction rate;
	double * N, * y;
	int species = 0, reactions = 0;
	void * p = 0;
	MassActionNetwork * net1;
	ProteinInteractionNetwork * net2;
	GeneRegulationNetwork * net3;

	if (!r || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;

	stoic = stoicFunctions[r->type];
	N = stoic(r->network);

	rate = rateFunctions[r->type];
	
	p = r->network;
	species = getNumSpecies(r);
	reactions = getNumReactions(r);

	y = ODEsim2(species, reactions,	N, rate, iv, 0, time, dt, p);

	free(N);

	return y;
}

double * networkSteadyState( ReactionNetwork * r, double* iv)
{
	StoichiometryFunction stoic;
	double * N;
	PropensityFunction rate;
	double * y = 0;
	MassActionNetwork * net1;
	ProteinInteractionNetwork * net2;
	GeneRegulationNetwork * net3;
	void * p = 0;
	int species = 0, reactions = 0;

	if (!r || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	stoic = stoicFunctions[r->type];
	N = stoic(r->network);

	rate = rateFunctions[r->type];
	
	p = r->network;
	species = getNumSpecies(r);
	reactions = getNumReactions(r);

	y = steadyState2(species, reactions, N, rate, iv, p, 1.0E-3,10000.0,0.1);
	free(N);
	
	return y;
}

double * simulateNetworkStochastically( ReactionNetwork * r, double* iv, double time, int* sz)
{
	StoichiometryFunction stoic;
	double * N;
	PropensityFunction rate;
	double * y = 0;
	int species = 0, reactions = 0;
	void * p = 0;
	MassActionNetwork * net1;
	ProteinInteractionNetwork * net2;
	GeneRegulationNetwork * net3;

	if (!r || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	stoic = stoicFunctions[r->type];
	N = stoic(r->network);

	rate = rateFunctions[r->type];
	
	p = r->network;
	species = getNumSpecies(r);
	reactions = getNumReactions(r);

	y = SSA(species, reactions,	N, rate, iv, 0, time, 1.0E5, sz, p);

	free(N);
	return y;
}

void printNetwork(ReactionNetwork * r)
{
	PrintNetworkFunction f;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return;
	
	f = printNetworkFunctions[r->type];
	
	f(r->network);
}

void printNetworkToFile(ReactionNetwork * r, char * filename)
{
	FILE *stream ;
	if ((stream = freopen(filename, "w", stdout)) == 0)
		return;

	printNetwork(r);
	
	stream = freopen("CON", "w", stdout);
}

int getNumSpecies(ReactionNetwork * r)
{
	GetNumSpeciesFunction f;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	f = getNumSpeciesFunctions[r->type];
	
	return (f(r->network));
}

int getNumReactions(ReactionNetwork * r)
{
	GetNumReactionsFunction f;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	f = getNumReactionsFunctions[r->type];
	
	return (f(r->network));
}

double* getReactionRates(ReactionNetwork * r, double* u)
{
	PropensityFunction f;
	double * rates;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	f = rateFunctions[r->type];
	rates = malloc( getNumReactions(r) );
	f(0,u,rates,r->network);

	return rates;
}

double* getStoichiometryMatrix(ReactionNetwork * r)
{
	StoichiometryFunction stoic;
	double * N;

	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	stoic = stoicFunctions[r->type];

	N = stoic(r->network);
	return N;
}

GAindividual mutateNetwork(GAindividual p)
{
	ReactionNetwork * r = (ReactionNetwork*)(p);
	GAMutateFnc f;

	GAindividual net;
	ReactionNetwork * r2;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return p;
	
	f = mutateFunctions[r->type];

	if (f)
	{
		net = f(r->network);
		r2 = malloc(sizeof(ReactionNetwork));
		r2->type = r->type;
		r2->network = net;
		return r;
	}
	return 0;
}

GAindividual crossoverNetwork(GAindividual p1, GAindividual p2)
{
	ReactionNetwork * r1 = (ReactionNetwork*)(p1);
	ReactionNetwork * r2 = (ReactionNetwork*)(p2);
	GACrossoverFnc f;
	GAindividual net;
	ReactionNetwork * r;

	if (r1->type != r2->type || r2->type < 0 || r2->type >= NUMBER_OF_NETWORK_TYPES)
		return mutateNetwork(p1);

	f = crossoverFunctions[r2->type];

	if (f)
	{
		net = f(r1->network,r2->network);
		r = malloc(sizeof(ReactionNetwork));
		r->type = r1->type;
		r->network = net;
		return r;
	}

	return mutateNetwork(p1);
}

void deleteNetwork(GAindividual p)
{
	ReactionNetwork * r = (ReactionNetwork*)(p);
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return;
	
	deleteFunctions[r->type](r->network);
	free(r);
}

GAindividual cloneNetwork(GAindividual p)
{
	ReactionNetwork * r = (ReactionNetwork*)(p);
	ReactionNetwork * r2;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return p;
	
	r2 = malloc(sizeof(ReactionNetwork));
	
	r2->type = r->type;
	r2->network = cloneFunctions[r->type](r->network);
	return r2;
}

void setFitnessFunction(GAFitnessFnc f)
{
	GAsetFitnessFunction(f);
}

void setInitialNetworkSize(int a,int b)
{
	setSizeForMassActionNetwork(a,b);
	setSizeForProteinInteractionNetwork(a,b);
	setSizeForGeneRegulationNetwork(a,b);
}

GApopulation evolveNetworks(int sz0,int sz1,int maxIter, GACallbackFnc callbackFunc)
{
	GApopulation P;
	if (!GAgetFitnessFunction()) return 0;

	P = randomNetworks(sz0);

	GAinit(&deleteNetwork, &cloneNetwork ,GAgetFitnessFunction(), &crossoverNetwork, &mutateNetwork, GAgetSelectionFunction());

	P = GArun(P,sz0,sz1,maxIter,callbackFunc);
	return P;
}

double compareSteadyStates(GAindividual p, double ** table, int rows, int inputs, int outputs)
{
	int i, j, k, l, m, cols, n;
	double * ss, * iv, closest, temp, sumOfSq;
	ReactionNetwork * r = (ReactionNetwork*)(p);
	MassActionNetwork * net1;
	ProteinInteractionNetwork * net2;
	GeneRegulationNetwork * net3;
	
	cols = inputs + outputs;
	
	n = getNumSpecies(r);
	
	if (n < cols) return 0.0; // not enough species
	
	SetFixedSpeciesFunction setFixed = setFixedSpeciesFunctions[r->type];
	
	for (i=0; i < inputs; ++i)
	{
		setFixed(r->network,i,1);
	}

	sumOfSq = 0.0;
	iv = malloc( n * sizeof(double) );
	
	for (m=0; m < rows; ++m)
	{	
		for (i=0; i < inputs; ++i)
			iv[i] = table[m][i];
		
		for (i=inputs; i < n; ++i)
			iv[i] = 0.0;
		
		ss = networkSteadyState(r,iv);
		
		if (ss) //error in simulation?
		{			
			for (i=0; i < outputs; ++i) //for each target output
			{
				closest = -1.0;
				
				for (j=inputs; j < n; ++j) //find best match
				{
					temp = (ss[j] - table[m][inputs+i]);
					if ((closest < 0.0) || ((temp*temp) < closest))
						closest = temp*temp;
				}
				
				sumOfSq += closest;
				
				if (closest < 0.0)
				{
					sumOfSq = -1.0;
					break;
				}
			}
			
			free(ss);		
		}
		else
		{
			sumOfSq = -1.0;
			break;
		}
	}
	
	free(iv);
	
	//restore the fixed
	for (i=0; i < inputs; ++i)
	{
		setFixed(r->network,i,0);
	}
	
	if (sumOfSq < 0.0) return 0.0;
	
	return (1.0 / (1.0 + sumOfSq));
}