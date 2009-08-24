#include "reactionNetwork.h"

/********************************************************

     Pointers to the existing function pointers
  
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

typedef double* (*StoichiometryFunction)(GAindividual);

static StoichiometryFunction stoicFunctions[] =
{
	&stoichiometryForMassActionNetwork,
	&stoichiometryForProteinInteractionNetwork, 
	&stoichiometryForGeneRegulationNetwork
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

/***********************************************************/

GApopulation randomNetworks(int sz0)
{
	int i, total = 0, k = 0;
	ReactionNetwork * rnet;
	double r;
	
	GApopulation P = malloc( sz0 * sizeof (GAindividual) );
	
	r = networkProbs[MASS_ACTION_NETWORK_INDEX] * sz0;
	
	GApopulation P1 = 0, P2 = 0, P3 = 0;
	
	if (r > 0 && r < 1.0)
	{
		P1 = randomMassActionNetworks((int)r);
		for (i=0; i < r; ++i)
		{
			rnet = malloc(sizeof(ReactionNetwork));
			(*rnet).type = MASS_ACTION_NETWORK_INDEX;
			(*rnet).network = P1[i];
			
			P[i] = rnet;
		}
		total += (int)r;
	}
	
	r = networkProbs[PROTEIN_INTERACTION_NETWORK_INDEX] * sz0;
	if (r > 0 && r < 1.0)
	{
		P2 = randomProteinInteractionNetworks((int)r);
		for (i=0; i < r; ++i)
		{
			rnet = malloc(sizeof(ReactionNetwork));
			(*rnet).type = PROTEIN_INTERACTION_NETWORK_INDEX;
			(*rnet).network = P2[i];
			
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
			(*rnet).type = GENE_REGULATION_NETWORK_INDEX;
			(*rnet).network = P3[i];
			
			P[i] = rnet;
		}
	}
	
	return P;
}

void setCrossoverFunction( int i, GACrossoverFnc f)
{
	if (i == MASS_ACTION_NETWORK)
	{
		crossoverFunctions[MASS_ACTION_NETWORK_INDEX] = f;
	}
	else
	if (i == PROTEIN_INTERACTION_NETWORK)
	{
		crossoverFunctions[PROTEIN_INTERACTION_NETWORK_INDEX] = f;
	}
	else
	if (i == GENE_REGULATION_NETWORK)
	{
		crossoverFunctions[GENE_REGULATION_NETWORK_INDEX] = f;
	}
}

void setMutationFunction( int i, GAMutateFnc f)
{
	if (i == MASS_ACTION_NETWORK)
	{
		mutateFunctions[MASS_ACTION_NETWORK_INDEX] = f;
	}
	else
	if (i == PROTEIN_INTERACTION_NETWORK)
	{
		mutateFunctions[PROTEIN_INTERACTION_NETWORK_INDEX] = f;
	}
	else
	if (i == GENE_REGULATION_NETWORK)
	{
		mutateFunctions[GENE_REGULATION_NETWORK_INDEX] = f;
	}
}

void setRatesFunction( int i, PropensityFunction f)
{
	if (i == MASS_ACTION_NETWORK)
	{
		rateFunctions[MASS_ACTION_NETWORK_INDEX] = f;
	}
	else
	if (i == PROTEIN_INTERACTION_NETWORK)
	{
		rateFunctions[PROTEIN_INTERACTION_NETWORK_INDEX] = f;
	}
	else
	if (i == GENE_REGULATION_NETWORK)
	{
		rateFunctions[GENE_REGULATION_NETWORK_INDEX] = f;
	}
}

void setStoichiometryFunction( int i, double* (*f)(GAindividual) )
{
	if (i == MASS_ACTION_NETWORK)
	{
		stoicFunctions[MASS_ACTION_NETWORK_INDEX] = f;
	}
	else
	if (i == PROTEIN_INTERACTION_NETWORK)
	{
		stoicFunctions[PROTEIN_INTERACTION_NETWORK_INDEX] = f;
	}
	else
	if (i == GENE_REGULATION_NETWORK)
	{
		stoicFunctions[GENE_REGULATION_NETWORK_INDEX] = f;
	}
}

void setNetworkType(int i)
{
	for (i=0; i < NUMBER_OF_NETWORK_TYPES; ++i) networkProbs[i] = 0.0;
	setNetworkTypeProbability(i, 1.0 );
}

void setNetworkTypeProbability(int i, double p)
{
	if (i | MASS_ACTION_NETWORK)
		networkProbs[ MASS_ACTION_NETWORK_INDEX ] = p;
		
	if (i | PROTEIN_INTERACTION_NETWORK)
		networkProbs[ PROTEIN_INTERACTION_NETWORK_INDEX ] = p;
		
	if (i | GENE_REGULATION_NETWORK)
		networkProbs[ GENE_REGULATION_NETWORK_INDEX ] = p;
		
	double total = 0.0;
	
	for (i=0; i < NUMBER_OF_NETWORK_TYPES; ++i) total += networkProbs[i];
	
	if (total <= 0.0)
		for (i=0; i < NUMBER_OF_NETWORK_TYPES; ++i) networkProbs[i] = 0.33333;
	else
		for (i=0; i < NUMBER_OF_NETWORK_TYPES; ++i) networkProbs[i] /= total;
}

double * simulateNetworkODE( ReactionNetwork * r, double* iv, double time, double dt)
{	
	if (!r) return 0;
	
	StoichiometryFunction stoic = stoicFunctions[(*r).type];
	double * N = stoic((*r).network);
	
	PropensityFunction rate = rateFunctions[(*r).type];
	double * y = 0;
	
	int species = 0, 
		reactions = 0;
	
	void * p = 0;
	
	if ((*r).type == MASS_ACTION_NETWORK_INDEX)
	{	
		MassActionNetwork * net = (MassActionNetwork*)((*r).network);
		species = (*net).species;
		reactions = (*net).reactions;
		p = net;
	}
	
	if ((*r).type == PROTEIN_INTERACTION_NETWORK_INDEX)
	{
		ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)((*r).network);
		species = (*net).species;
		reactions = 2 * (*net).species;
		p = net;
	}
	
	if ((*r).type == GENE_REGULATION_NETWORK_INDEX)
	{
		GeneRegulationNetwork * net = (GeneRegulationNetwork*)((*r).network);
		species = (*net).species;
		reactions = 2 * (*net).species;
		p = net;
	}
	
	y = ODEsim2(species, reactions,	N, rate, iv, 0, time, dt, p);
	
	free(N);
	
	return y;
}

double * networkSteadyState( ReactionNetwork * r, double* iv)
{
	if (!r) return 0;
	
	StoichiometryFunction stoic = stoicFunctions[(*r).type];
	double * N = stoic((*r).network);
	
	PropensityFunction rate = rateFunctions[(*r).type];
	double * y = 0;
	
	void * p = 0;
	
	int species = 0, 
		reactions = 0;
	
	if ((*r).type == MASS_ACTION_NETWORK_INDEX)
	{
		MassActionNetwork * net = (MassActionNetwork*)((*r).network);
		species = (*net).species;
		reactions = (*net).reactions;
		p = net;
	}
	
	if ((*r).type == PROTEIN_INTERACTION_NETWORK_INDEX)
	{
		ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)((*r).network);
		species = (*net).species;
		reactions = 2 * (*net).species;
		p = net;
	}
	
	if ((*r).type == GENE_REGULATION_NETWORK_INDEX)
	{
		GeneRegulationNetwork * net = (GeneRegulationNetwork*)((*r).network);
		species = (*net).species;
		reactions = 2 * (*net).species;
		p = net;
	}
	
	y = steadyState2(species, reactions, N, rate, iv, p, 1.0E-3,10000.0,0.1);
	free(N);
	return y;
}

double * simulateNetworkStochastically( ReactionNetwork * r, double* iv, double time, int* sz)
{
	if (!r) return 0;
	
	StoichiometryFunction stoic = stoicFunctions[(*r).type];
	double * N = stoic((*r).network);
	
	PropensityFunction rate = rateFunctions[(*r).type];
	double * y = 0;

	int species = 0, reactions = 0;
	
	void * p = 0;
	
	if ((*r).type == MASS_ACTION_NETWORK_INDEX)
	{	
		MassActionNetwork * net = (MassActionNetwork*)((*r).network);
		species = (*net).species;
		reactions = (*net).reactions;
		p = net;
	}
	
	if ((*r).type == PROTEIN_INTERACTION_NETWORK_INDEX)
	{
		ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)((*r).network);
		species = (*net).species;
		reactions = 2 * (*net).species;
		p = net;
	}
	
	if ((*r).type == GENE_REGULATION_NETWORK_INDEX)
	{
		GeneRegulationNetwork * net = (GeneRegulationNetwork*)((*r).network);
		species = (*net).species;
		reactions = 2 * (*net).species;
		p = net;
	}
	
	y = SSA(species, reactions,	N, rate, iv, 0, time, 1.0E5, sz, p);
	
	free(N);
	return y;
}

void printNetwork(ReactionNetwork * r)
{
	if ((*r).type == MASS_ACTION_NETWORK_INDEX)
	{	
		printMassActionNetwork((MassActionNetwork*)(*r).network);		
	}
	
	if ((*r).type == PROTEIN_INTERACTION_NETWORK_INDEX)
	{
		printProteinInteractionNetwork((ProteinInteractionNetwork*)(*r).network);
	}
	
	if ((*r).type == GENE_REGULATION_NETWORK_INDEX)
	{
		printGeneRegulationNetwork((GeneRegulationNetwork*)(*r).network);
	}
}

int getNumSpecies(ReactionNetwork * r)
{
	if (!r) return 0;
	if ((*r).type == MASS_ACTION_NETWORK_INDEX)
	{	
		MassActionNetwork * net = (MassActionNetwork*)((*r).network);
		return (*net).species;
	}
	
	if ((*r).type == PROTEIN_INTERACTION_NETWORK_INDEX)
	{
		ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)((*r).network);
		return (*net).species;
	}
	
	if ((*r).type == GENE_REGULATION_NETWORK_INDEX)
	{
		GeneRegulationNetwork * net = (GeneRegulationNetwork*)((*r).network);
		return  (*net).species;
	}
	return 0;
}

int getNumReactions(ReactionNetwork * r)
{
	if (!r) return 0;
	if ((*r).type == MASS_ACTION_NETWORK_INDEX)
	{	
		MassActionNetwork * net = (MassActionNetwork*)((*r).network);
		return (*net).reactions;
	}
	
	if ((*r).type == PROTEIN_INTERACTION_NETWORK_INDEX)
	{
		ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)((*r).network);
		return (2 * (*net).species);
	}
	
	if ((*r).type == GENE_REGULATION_NETWORK_INDEX)
	{
		GeneRegulationNetwork * net = (GeneRegulationNetwork*)((*r).network);
		return (2 * (*net).species);
	}
	return 0;
}

double* getReactionRates(ReactionNetwork * r, double* u)
{
	if (!r) return 0;
	
	PropensityFunction f = rateFunctions[(*r).type];
	double * rates = malloc( getNumReactions(r) );
	
	f(0,u,rates,(*r).network);
	
	return rates;
}

double* getStoichiometryMatrix(ReactionNetwork * r)
{
	if (!r) return 0;
	StoichiometryFunction stoic = stoicFunctions[(*r).type];
	double * N = stoic((*r).network);
	return N;
}

GAindividual mutateNetwork(GAindividual p)
{
	ReactionNetwork * r = (ReactionNetwork*)(p);
	
	GAMutateFnc f = mutateFunctions[(*r).type];
	
	if (f)
	{
		GAindividual net = f((*r).network);
		ReactionNetwork * r2 = malloc(sizeof(ReactionNetwork));
		(*r2).type = (*r).type;
		(*r2).network = net;
		return r;
	}
	return 0;
}

GAindividual crossoverNetwork(GAindividual p1, GAindividual p2)
{
	ReactionNetwork * r1 = (ReactionNetwork*)(p1);
	ReactionNetwork * r2 = (ReactionNetwork*)(p2);
	
	if ((*r1).type != (*r2).type || (*r2).type < 0 || (*r2).type >= NUMBER_OF_NETWORK_TYPES)
		return mutateNetwork(p1);
	
	GACrossoverFnc f = crossoverFunctions[(*r2).type];
	
	if (f)
	{
		GAindividual net = f((*r1).network,(*r2).network);
		ReactionNetwork * r = malloc(sizeof(ReactionNetwork));
		(*r).type = (*r1).type;
		(*r).network = net;
		return r;
	}
	
	return mutateNetwork(p1);
}

void deleteNetwork(GAindividual p)
{
	ReactionNetwork * r = (ReactionNetwork*)(p);
	deleteFunctions[(*r).type]((*r).network);
	free(r);
}

GAindividual cloneNetwork(GAindividual p)
{
	ReactionNetwork * r = (ReactionNetwork*)(p);
	ReactionNetwork * r2 = malloc(sizeof(ReactionNetwork));
	(*r2).type = (*r).type;
	(*r2).network = cloneFunctions[(*r).type]((*r).network);
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
	if (!GAgetFitnessFunction()) return 0;
	
	GApopulation P = randomNetworks(sz0);
	
	GAinit(&deleteNetwork, &cloneNetwork ,GAgetFitnessFunction(), &crossoverNetwork, &mutateNetwork, GAgetSelectionFunction());
	
	P = GArun(P,sz0,sz1,maxIter,callbackFunc);
	return P;
}
