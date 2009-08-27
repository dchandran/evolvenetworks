#include "reactionNetwork.h"

/********************************************************

		Pointers to the existing functions

*********************************************************/

static int TRACK_NETWORK_PARENTS = 1;
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
typedef int (*GetNumSpeciesFunction)(GAindividual);
typedef int (*GetNumReactionsFunction)(GAindividual);

static StoichiometryFunction stoicFunctions[] =
{
	&stoichiometryForMassActionNetwork,
	&stoichiometryForProteinInteractionNetwork, 
	&stoichiometryForGeneRegulationNetwork
};

static GetNumSpeciesFunction getNumSpeciesFunctions[] =
{
	&getNumSpeciesForMassActionNetwork,
	&getNumSpeciesForProteinInteractionNetwork, 
	&getNumSpeciesForGeneRegulationNetwork
};

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
	int i, total = 0, r, k;
	ReactionNetwork * rnet;
	GApopulation P1 = 0, P2 = 0, P3 = 0, P;

	P = (GAindividual*) malloc( sz0 * sizeof (GAindividual) );

	r = (int)(networkProbs[MASS_ACTION_NETWORK] * sz0);
	
	if (r > 0 && r <= sz0)
	{
		P1 = randomMassActionNetworks(r);
		for (i=0; i < r; ++i)
		{
			rnet = (ReactionNetwork*) malloc(sizeof(ReactionNetwork));
			rnet->type = MASS_ACTION_NETWORK;
			rnet->network = P1[i];
			rnet->id = i;
			rnet->parents = 0;

			P[i] = rnet;
		}
		total += r;
	}

	r = (int)(networkProbs[PROTEIN_INTERACTION_NETWORK] * sz0);
	if (r > 0 && r <= sz0)
	{
		P2 = randomProteinInteractionNetworks(r);
		for (i=0; i < r; ++i)
		{
			rnet = (ReactionNetwork*) malloc(sizeof(ReactionNetwork));
			rnet->type = PROTEIN_INTERACTION_NETWORK;
			rnet->network = P2[i];
			rnet->parents = 0;

			rnet->id = i+total;
			P[i+total] = rnet;
		}
		total += r;
	}

	if (total < sz0)
	{
		k = sz0 - total;
		P3 = randomGeneRegulationNetworks(k);
		for (i=0; i < k; ++i)
		{
			rnet = (ReactionNetwork*) malloc(sizeof(ReactionNetwork));
			rnet->type = GENE_REGULATION_NETWORK;
			rnet->network = P3[i];
			rnet->parents = 0;

			rnet->id = i+total;
			P[i+total] = rnet;
		}
	}

	return P;
}

/*******************************************************************

The following functions rely entirely on 
the array of function pointers defined at 
the beginning of this file

********************************************************************/

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

double * simulateNetworkODE( GAindividual individual, double* iv, double time, double dt)
{	
	ReactionNetwork * r = (ReactionNetwork*)individual;
	StoichiometryFunction stoic;
	PropensityFunction rate;
	double * N, * y;
	int species = 0, reactions = 0;
	void * p = 0;

	if (!r || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	stoic = stoicFunctions[r->type];
	rate = rateFunctions[r->type];
	
	p = r->network;
	species = getNumSpecies(r);
	reactions = getNumReactions(r);

	if (species == 0 || reactions == 0) return 0;
	N = stoic(r->network);
	y = ODEsim2(species, reactions,	N, rate, iv, 0, time, dt, p);

	free(N);

	return y;
}

double * networkSteadyState( GAindividual individual, double* iv)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	StoichiometryFunction stoic;
	double * N;
	PropensityFunction rate;
	double * y = 0;
	void * p = 0;
	int species = 0, reactions = 0;

	if (!r || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	stoic = stoicFunctions[r->type];
	
	rate = rateFunctions[r->type];
	
	p = r->network;
	species = getNumSpecies(r);
	reactions = getNumReactions(r);

	if (species == 0 || reactions == 0) return 0;
	N = stoic(r->network);
	y = steadyState2(species, reactions, N, rate, iv, p, 1.0E-3,10000.0,0.1);
	free(N);
	
	return y;
}

double * simulateNetworkStochastically( GAindividual individual, double* iv, double time, int* sz)
{	
	ReactionNetwork * r = (ReactionNetwork*)individual;
	StoichiometryFunction stoic;
	double * N;
	PropensityFunction rate;
	double * y = 0;
	int species = 0, reactions = 0;
	void * p = 0;

	if (!r || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	stoic = stoicFunctions[r->type];
	rate = rateFunctions[r->type];
	
	p = r->network;
	species = getNumSpecies(r);
	reactions = getNumReactions(r);

	if (species == 0 || reactions == 0) return 0;
	N = stoic(r->network);
	y = SSA(species, reactions,	N, rate, iv, 0.0, time, 1000000, sz, p);
	free(N);
	return y;
}

void printNetwork(GAindividual individual)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	PrintNetworkFunction f;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return;
	
	f = printNetworkFunctions[r->type];
	
	f(r->network);
}

void printNetworkToFile(GAindividual individual, char * filename)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	FILE *stream ;
	if ((stream = freopen(filename, "w", stdout)) == 0)
		return;

	printNetwork(r);
	
	stream = freopen("CON", "w", stdout);
}

int getNumSpecies(GAindividual individual)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	GetNumSpeciesFunction f;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	f = getNumSpeciesFunctions[r->type];
	
	return (f(r->network));
}

int getNumReactions(GAindividual individual)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	GetNumReactionsFunction f;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	f = getNumReactionsFunctions[r->type];
	
	return (f(r->network));
}

double* getReactionRates(GAindividual individual, double* u)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	PropensityFunction f;
	double * rates;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	f = rateFunctions[r->type];
	rates = (double*) malloc( getNumReactions(r) );
	f(0,u,rates,r->network);

	return rates;
}

double* getStoichiometryMatrix(GAindividual individual)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
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
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return p;
	
	f = mutateFunctions[r->type];

	if (f)
	{
		net = f(r->network);
		r->network = net;
		
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
	int i,j,k,sz1,sz2;

	if (r1->type != r2->type || r2->type < 0 || r2->type >= NUMBER_OF_NETWORK_TYPES)
		return mutateNetwork(p1);

	f = crossoverFunctions[r2->type];

	if (f)
	{
		net = f(r1->network,r2->network);
		r = (ReactionNetwork*)malloc(sizeof(ReactionNetwork));
		r->type = r1->type;
		r->network = net;
		r->id = r1->id;
		i = j = sz1 = sz2 = 0;
		if (r1->parents)
		{
			while (r1->parents[i]) 
			{
				j = 0;
				while (j < i && r1->parents[i] != r1->parents[j])
					++j;
				if (j == i)
					++sz1;
				++i;
			}
		}
		i = 0;
		if (r2->parents)
		{
			while (r2->parents[i]) 
			{
				j = 0;
				while (j < i && r2->parents[i] != r2->parents[j])
					++j;
				if (j == i)
					++sz2;
				++i;
			}
		}

		r->parents = 0;
		if (TRACK_NETWORK_PARENTS)
		{
			if ((sz1+sz2) > 0)
			{
				r->parents = (int*) malloc((sz1+sz2+1)*sizeof(int));
				r->parents[sz1+sz2] = 0;

				if (sz1 > 0)
					for (k=0; k < sz1; ++k)
						r->parents[k] = r1->parents[k];
				if (sz2 > 0)
					for (k=0; k < sz2; ++k)
						r->parents[k+sz1] = r2->parents[k];
			}
			else
			{
				r->parents = (int*) malloc((3)*sizeof(int));
				r->parents[2] = 0;

				r->parents[0] = r1->id;
				r->parents[1] = r2->id;
			}
		}
		return r;
	}

	return mutateNetwork(p1);
}

void deleteNetwork(GAindividual p)
{
	ReactionNetwork * r = (ReactionNetwork*)(p);
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return;
	
	deleteFunctions[r->type](r->network);
	
	if (r->parents)
		free(r->parents);

	free(r);
}

GAindividual cloneNetwork(GAindividual p)
{
	ReactionNetwork * r = (ReactionNetwork*)(p);
	ReactionNetwork * r2;
	int i,j;

	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return p;
	
	r2 = (ReactionNetwork*) malloc(sizeof(ReactionNetwork));
	
	i = 0;
	if (r->parents)
	{
		while (r->parents[i]) ++i;
	}

	r2->parents = 0;
	if (TRACK_NETWORK_PARENTS && (i > 0))
	{
		r2->parents = (int*) malloc((i+1)*sizeof(int));
		for (j=0; j < i; ++j)
			r2->parents[j] = r->parents[j];
		r2->parents[i] = 0;
	}
	r2->id = r->id;
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

/*******************************
  Related to lineage tracking
*******************************/

void lineageTrackingON()
{
	TRACK_NETWORK_PARENTS = 1;
}

void lineageTrackingOFF()
{
	TRACK_NETWORK_PARENTS = 0;
}

void setID(GAindividual individual,int i)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	if (r)
		r->id = i;
}

int getID(GAindividual individual)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	if (!r) return -1;

	return r->id;
}

int* getParentIDs(GAindividual individual)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;

	if (!r) return 0;
	return r->parents;
}

/*******************************
   Special fitness function
*******************************/

double compareSteadyStates(GAindividual p, double ** table, int rows, int inputs, int outputs)
{
	int i, j, m, cols, n, best;
	double * ss, * iv, closest, temp, sumOfSq;
	ReactionNetwork * r = (ReactionNetwork*)(p);
	SetFixedSpeciesFunction setFixed;
	
	cols = inputs + outputs;
	
	n = getNumSpecies(r);
	
	if (n < cols) return 0.0; // not enough species
	
	setFixed = setFixedSpeciesFunctions[r->type];
	
	for (i=0; i < inputs; ++i)
	{
		setFixed(r->network,i,1);
	}

	sumOfSq = 0.0;
	iv = (double*) malloc( n * sizeof(double) );

	best = -1;
	
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
				if (best < 0)
				{
					closest = -1.0;
					
					for (j=inputs; j < n; ++j) //find best match
					{
						temp = (ss[j] - table[m][inputs+i]);
						if ((closest < 0.0) || ((temp*temp) < closest))
						{
							closest = temp*temp;
							best = j;
						}
					}
				}
				else
				{
					j = best;
					temp = (ss[j] - table[m][inputs+i]);
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
