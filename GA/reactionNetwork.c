#include "reactionNetwork.h"

/********************************************************

		For keeping log file

*********************************************************/

static FILE * LOGFILE = 0;

static int PRINT_SEEDS = 1; 

static int PRINT_EACH_FITNESS = 1; 
static int PRINT_FINAL_FITNESS = 1;

static int PRINT_EACH_SCRIPT = 0;
static int PRINT_FINAL_SCRIPT = 1;

static int PRINT_EACH_SIZE = 1;
static int PRINT_FINAL_SIZE = 1;

static int PRINT_EACH_BEST_LINEAGE = 1;
static int PRINT_FINAL_BEST_LINEAGE = 1;

static int PRINT_EACH_ALL_FITNESS = 0;
static int PRINT_FINAL_ALL_FITNESS = 1;

static int PRINT_EACH_ALL_LINEAGE = 1;
static int PRINT_FINAL_ALL_LINEAGE = 1;

static GACallbackFnc USER_CALLBACK_FNC = 0;

/********************************************************

		Global parameters

*********************************************************/

static double MUTATE_INIT_VALUE_PROB = 1.0;
static double AVG_INIT_VALUES = 2.0;
static int NUMBER_OF_NETWORK_TYPES = 4;
static double CROSSOVER_PROB = 1.0;
static int TRACK_NETWORK_PARENTS = 0;

void setCrossoverRate(double d)
{
	CROSSOVER_PROB = d;
}

void setAverageInitialValue(double d)
{
	if (d > 0.0)
		AVG_INIT_VALUES = d;
}

void setMutationRateOfInitialValues(double d)
{
	MUTATE_INIT_VALUE_PROB = d;
}

static double SS_FUNC_ERROR_TOLERANCE = 1.0E-3;
static double SS_FUNC_DELTA_TIME = 0.1;
static double SS_FUNC_MAX_TIME = 10000.0;

void configureSteadyStateFunction(double tolerance, 
									double delta,
									double maxTime)
{
	SS_FUNC_ERROR_TOLERANCE = tolerance;
	SS_FUNC_DELTA_TIME = delta;
	SS_FUNC_MAX_TIME = maxTime;
}

/********************************************************

		Pointers to the existing functions

*********************************************************/

static double networkProbs[] =
{
	0.33333333333,
	0.33333333333,
	0.33333333333,
	0.33333333333
};

static PropensityFunction rateFunctions[] =
{
	&ratesForMassActionNetwork,
	&ratesForEnzymeNetwork,
	&ratesForGeneRegulationNetwork,
	&ratesForProteinInteractionNetwork
};

static GACrossoverFnc crossoverFunctions[] =
{
	&crossoverMassActionNetwork, 
	&crossoverEnzymeNetwork,
	&crossoverGeneRegulationNetwork,
	&crossoverProteinInteractionNetwork
};

static GAMutateFnc mutateFunctions[] =
{
	&mutateMassActionNetwork, 
	&mutateEnzymeNetwork,
	&mutateGeneRegulationNetwork,
	&mutateProteinInteractionNetwork
};

static GADeleteFnc deleteFunctions[] =
{
	&deleteMassActionNetwork,
	&deleteEnzymeNetwork,
	&deleteGeneRegulationNetwork,
	&deleteProteinInteractionNetwork
};

static GACloneFnc cloneFunctions[] =
{
	&cloneMassActionNetwork,
	&cloneEnzymeNetwork,
	&cloneGeneRegulationNetwork,
	&cloneProteinInteractionNetwork
};

typedef double* (*StoichiometryFunction)(GAindividual);
typedef int (*GetNumSpeciesFunction)(GAindividual);
typedef int (*GetNumReactionsFunction)(GAindividual);

static StoichiometryFunction stoicFunctions[] =
{
	&stoichiometryForMassActionNetwork,
	&stoichiometryForEnzymeNetwork,
	&stoichiometryForGeneRegulationNetwork,
	&stoichiometryForProteinInteractionNetwork
};

static GetNumSpeciesFunction getNumSpeciesFunctions[] =
{
	&getNumSpeciesForMassActionNetwork,
	&getNumSpeciesForEnzymeNetwork, 
	&getNumSpeciesForGeneRegulationNetwork,
	&getNumSpeciesForProteinInteractionNetwork
};

static GetNumReactionsFunction getNumReactionsFunctions[] =
{
	&getNumReactionsForMassActionNetwork,
	&getNumReactionsForEnzymeNetwork, 
	&getNumReactionsForGeneRegulationNetwork,
	&getNumReactionsForProteinInteractionNetwork
};

typedef void (*PrintNetworkFunction)(FILE *, GAindividual);

static PrintNetworkFunction printNetworkFunctions[] =
{
	&printMassActionNetwork,
	&printEnzymeNetwork,
	&printGeneRegulationNetwork,
	&printProteinInteractionNetwork
};

typedef void (*SetFixedSpeciesFunction)(GAindividual,int,int);

static SetFixedSpeciesFunction setFixedSpeciesFunctions[] =
{
	&setFixedSpeciesForMassActionNetwork,
	&setFixedSpeciesForEnzymeNetwork,
	&setFixedSpeciesForGeneRegulationNetwork,
	&setFixedSpeciesForProteinInteractionNetwork
};

/***********************************************************

    The following functions need to be modified
	when adding a new network type

************************************************************/

void setNetworkSize(int n0, int n1, int r0, int r1)
{
	setSizeForMassActionNetwork(n0,n1,r0,r1);
	setSizeForEnzymeNetwork(n0,n1,r0,r1);
	setSizeForProteinInteractionNetwork(n0,n1,r0,r1);
	setSizeForGeneRegulationNetwork(n0,n1,r0,r1);
}

GApopulation randomNetworks(int sz0)
{
	int i, j, r, k, n, total = 0;
	ReactionNetwork * rnet;
	GApopulation P0 = 0, P;

	P = (GAindividual*) malloc( sz0 * sizeof (GAindividual) );

	r = (int)(networkProbs[MASS_ACTION_NETWORK] * sz0);
	
	if (r > 0 && r <= sz0)
	{
		P0 = randomMassActionNetworks(r);
		for (i=0; i < r; ++i)
		{
			rnet = (ReactionNetwork*) malloc(sizeof(ReactionNetwork));
			rnet->type = MASS_ACTION_NETWORK;
			rnet->network = P0[i];
			//rnet->id = i+total;
			//rnet->parents = 0;
			n = getNumSpecies(rnet);
			rnet->initialValues = (double*)malloc(n*sizeof(double));
			for (j=0; j < n; ++j)
				rnet->initialValues[j] = AVG_INIT_VALUES * mtrand();

			P[i] = rnet;
		}
		total += r;
	}

	r = (int)(networkProbs[ENZYME_NETWORK] * sz0);
	
	if (r > 0 && r <= sz0)
	{
		P0 = randomEnzymeNetworks(r);
		for (i=0; i < r; ++i)
		{
			rnet = (ReactionNetwork*) malloc(sizeof(ReactionNetwork));
			rnet->type = ENZYME_NETWORK;
			rnet->network = P0[i];
			//rnet->id = i+total;
			//rnet->parents = 0;
			n = getNumSpecies(rnet);
			rnet->initialValues = (double*)malloc(n*sizeof(double));
			for (j=0; j < n; ++j)
				rnet->initialValues[j] = AVG_INIT_VALUES * mtrand();

			P[i] = rnet;
		}
		total += r;
	}

	r = (int)(networkProbs[PROTEIN_INTERACTION_NETWORK] * sz0);
	if (r > 0 && r <= sz0)
	{
		P0 = randomProteinInteractionNetworks(r);
		for (i=0; i < r; ++i)
		{
			rnet = (ReactionNetwork*) malloc(sizeof(ReactionNetwork));
			rnet->type = PROTEIN_INTERACTION_NETWORK;
			//rnet->network = P0[i];
			//rnet->parents = 0;
			n = getNumSpecies(rnet);
			rnet->initialValues = (double*)malloc(n*sizeof(double));
			for (j=0; j < n; ++j)
				rnet->initialValues[j] = AVG_INIT_VALUES * mtrand();

			//rnet->id = i+total;
			P[i+total] = rnet;
		}
		total += r;
	}

	if (total < sz0)
	{
		k = sz0 - total;
		P0 = randomGeneRegulationNetworks(k);
		for (i=0; i < k; ++i)
		{
			rnet = (ReactionNetwork*) malloc(sizeof(ReactionNetwork));
			rnet->type = GENE_REGULATION_NETWORK;
			//rnet->network = P0[i];
			//rnet->parents = 0;
			n = getNumSpecies(rnet);
			rnet->initialValues = (double*)malloc(n*sizeof(double));
			for (j=0; j < n; ++j)
				rnet->initialValues[j] = AVG_INIT_VALUES * mtrand();

			//rnet->id = i+total;
			P[i+total] = rnet;
		}
	}

	return P;
}

/*******************************************************************

	The following functions rely entirely on 
	the array of function pointers defined at 
	the beginning of this file. They do NOT need
	to be modified when adding a new network type

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

double * simulateNetworkODE( GAindividual individual, double time, double dt)
{	
	ReactionNetwork * r = (ReactionNetwork*)individual;
	StoichiometryFunction stoic;
	PropensityFunction rate;
	double * N, * y;
	int species = 0, reactions = 0;
	void * p = 0;
	double* iv;

	if (!r || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	iv = r->initialValues;
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

double * networkSteadyState( GAindividual individual )
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	StoichiometryFunction stoic;
	double * N;
	PropensityFunction rate;
	double * y = 0;
	void * p = 0;
	int species = 0, reactions = 0;
	double* iv;

	if (!r || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	iv = r->initialValues;
	stoic = stoicFunctions[r->type];	
	rate = rateFunctions[r->type];
	
	p = r->network;
	species = getNumSpecies(r);
	reactions = getNumReactions(r);

	if (species == 0 || reactions == 0) return 0;
	N = stoic(r->network);
	y = steadyState2(species, reactions, N, rate, iv, p, SS_FUNC_ERROR_TOLERANCE,SS_FUNC_MAX_TIME,SS_FUNC_DELTA_TIME);
	free(N);
	
	return y;
}

double * simulateNetworkStochastically( GAindividual individual, double time, int* sz)
{	
	ReactionNetwork * r = (ReactionNetwork*)individual;
	StoichiometryFunction stoic;
	double * N;
	PropensityFunction rate;
	double * y = 0, * iv;
	int species = 0, reactions = 0;
	void * p = 0;

	if (!r || (r->type > NUMBER_OF_NETWORK_TYPES)) return 0;
	
	iv = r->initialValues;
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

void printNetwork(FILE * stream, GAindividual individual)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	PrintNetworkFunction f;
	int n,i;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return;
	
	f = printNetworkFunctions[r->type];
	
	f(stream,r->network);

	if (r->initialValues)
	{
		fprintf(stream,"\n");
		n = getNumSpecies(r);
		for (i=0; i < n; ++i)
		{
			fprintf(stream,"s%i = %lf;\n",i+1,r->initialValues[i]);
		}
	}
}

void printNetworkToFile(char * filename, GAindividual individual)
{
	ReactionNetwork * r = (ReactionNetwork*)individual;
	PrintNetworkFunction f;
	FILE *stream;
	int i, n;

	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return;

	stream = fopen(filename, "w");
	if (!stream) return;
	f = printNetworkFunctions[r->type];
	f(stream,r->network);

	if (r->initialValues)
	{
		fprintf(stream,"\n");
		n = getNumSpecies(r);
		for (i=0; i < n; ++i)
		{
			fprintf(stream, "s%i = %lf;\n",i+1,r->initialValues[i]);
		}
	}
	fclose(stream);
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
	int n0, n1, i;
	double * iv;
	GAMutateFnc f;

	GAindividual net;
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return p;
	
	f = mutateFunctions[r->type];
	n0 = getNumSpecies(r);

	if (f)
	{
		net = f(r->network);
		r->network = net;
		n1 = getNumSpecies(r);
		if (n0 != n1)
		{
			iv = r->initialValues;
			r->initialValues = (double*)malloc(n1 * sizeof(double));
			for (i=0; i < n0 && i < n1; ++i)
			{
				r->initialValues[i] = iv[i];
			}
			for (i=n0; i < n1; ++i)
			{
				r->initialValues[i] = AVG_INIT_VALUES * mtrand();
			}
			free(iv);
		}

		if (mtrand() < MUTATE_INIT_VALUE_PROB) 
		{
			r->initialValues[ (int)(mtrand() * getNumSpecies(r)) ] *= 2.0 * mtrand();
		}	
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
	int i,j,k,sz1,sz2,*p;

	if (mtrand() > CROSSOVER_PROB || r1->type != r2->type || r2->type < 0 || r2->type >= NUMBER_OF_NETWORK_TYPES)
		return mutateNetwork(p1);

	f = crossoverFunctions[r2->type];

	if (f)
	{
		net = f(r1->network,r2->network);
		r = (ReactionNetwork*)malloc(sizeof(ReactionNetwork));
		r->type = r1->type;
		r->network = net;
		//r->id = r1->id;

		sz1 = getNumSpecies(r1);
		sz2 = getNumSpecies(r2);
		j = getNumSpecies(r);

		r->initialValues = (double*)malloc(j * sizeof(double));

		for (i=0; i < sz1 && i < j; ++i)
			r->initialValues[i] = r1->initialValues[i];
		for (i=0; i < sz2 && (sz1+i) < j; ++i)
			r->initialValues[sz1+i] = r2->initialValues[i];
		for (i=sz1+sz2; i < j; ++i)
			r->initialValues[i] = AVG_INIT_VALUES * mtrand();

		/*
		i = j = sz1 = sz2 = 0;
		if (r1->parents)
		{
			while (r1->parents[i]) 
			{
				++i;
			}
		}
		
		j = 0;
		if (r2->parents)
		{
			while (r2->parents[j]) 
			{
				++j;
			}
		}
		
		sz1 = i + j;		

		r->parents = 0;
		if (TRACK_NETWORK_PARENTS)
		{
			if (sz1 > 0)
			{
				p = (int*)malloc(sz1 * sizeof(int));
				
				for (k=0; k < i; ++k)
					p[k] = r1->parents[k];
					
				for (k=0; k < j; ++k)
					p[k+i] = r2->parents[k];
				
				sz2 = 0;
				
				for (i=0; i < sz1; ++i)
				{
					for (j=0; j < i; ++j)
						if (p[j] == p[i]) break;
					if (i == j)
						++sz2;
				}
				
				r->parents = (int*) malloc((sz2+1)*sizeof(int));
				r->parents[sz2] = 0;

				k = 0;
				for (i=0; i < sz1; ++i)
				{
					for (j=0; j < i; ++j)
						if (p[j] == p[i]) break;
					if (i == j)
					{
						r->parents[k] = p[i];
						++k;
					}
				}
				free(p);
			}
			else
			{
				r->parents = (int*) malloc((3)*sizeof(int));
				r->parents[2] = 0;

				r->parents[0] = r1->id;
				r->parents[1] = r2->id;
			}
		}
		*/
		return r;
	}

	return mutateNetwork(p1);
}

void deleteNetwork(GAindividual p)
{
	ReactionNetwork * r = (ReactionNetwork*)(p);
	
	if (!r || (r->type < 0) || (r->type > NUMBER_OF_NETWORK_TYPES)) return;
	
	deleteFunctions[r->type](r->network);
	
	//if (r->parents)
		//free(r->parents);

	if (r->initialValues)
		free(r->initialValues);

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
	/*if (r->parents)
		while (r->parents[i]) 
			++i;
	r2->parents = 0;
	if (TRACK_NETWORK_PARENTS && (i > 0))
	{
		r2->parents = (int*) malloc((i+1)*sizeof(int));
		for (j=0; j < i; ++j)
			r2->parents[j] = r->parents[j];
		r2->parents[i] = 0;
	}
	r2->id = r->id;
	*/
	r2->type = r->type;
	
	r2->network = cloneFunctions[r->type](r->network);

	j = getNumSpecies(r2);
	r2->initialValues = (double*)malloc(j * sizeof(double));
	for (i=0; i < j; ++i)
		r2->initialValues[i] = r->initialValues[i];
	

	return r2;
}

void setFitnessFunction(GAFitnessFnc f)
{
	GAsetFitnessFunction(f);
}

void setCallbackFunction(GACallbackFnc f)
{
	GAsetCallbackFunction(f);
}

static int callBackWithLogKeeping(int iter,GApopulation pop,int popSz)
{
	GAFitnessFnc fitness = GAgetFitnessFunction();
	double f;
	int i,j,k,*parents, num = 10*popSz, max = 0;
	int * temp = 0, * ids = 0;
	GAindividual * p;
	
	if (iter == 0) //header
	{
		printf("gen");
		fprintf(LOGFILE,"gen");
		if (PRINT_EACH_FITNESS && !PRINT_EACH_ALL_FITNESS)
		{
			printf("\tfitness ");
			fprintf(LOGFILE,"\tfitness ");
		}
		
		if (PRINT_EACH_ALL_FITNESS)
		{
			for (i=0; i < popSz; ++i)
			{
				//printf("\tfitness_%i",i);
				fprintf(LOGFILE,"\tfitness_%i",i);
			}
		}
	
		if (PRINT_EACH_SIZE)
		{
			printf("\tspecies \treactions");
			fprintf(LOGFILE,"\tspecies \treactions");
		}
	
		if (TRACK_NETWORK_PARENTS && (PRINT_EACH_BEST_LINEAGE || PRINT_EACH_ALL_LINEAGE))
		{
			//printf("\tparents");
			fprintf(LOGFILE,"\tparents");
		}
		
		printf("\n");
		fprintf(LOGFILE,"\n");
		
		printf("---");
		fprintf(LOGFILE,"---");

		if (PRINT_EACH_FITNESS && !PRINT_EACH_ALL_FITNESS)
		{
			printf("\t ------ ");
			fprintf(LOGFILE,"\t ------ ");
		}
		
		if (PRINT_EACH_ALL_FITNESS)
		{
			for (i=0; i < popSz; ++i)
			{
				//printf("\t----------",i);
				fprintf(LOGFILE,"\t -------- ",i);
			}
		}
	
		if (PRINT_EACH_SIZE)
		{
			printf("\t-------\t---------");
			fprintf(LOGFILE,"\t-------\t---------");
		}
	
		if (TRACK_NETWORK_PARENTS && (PRINT_EACH_BEST_LINEAGE || PRINT_EACH_ALL_LINEAGE))
		{
			//printf("\t-------");
			fprintf(LOGFILE,"\t-------");
		}
		
		printf("\n");
		fprintf(LOGFILE,"\n");
	}

	if (TRACK_NETWORK_PARENTS && (PRINT_EACH_BEST_LINEAGE || PRINT_EACH_ALL_LINEAGE))
	{
		ids = (int*)malloc(num * sizeof(int));

		for (i=0; i < num; ++i)
			ids[i] = 0;

		for (i=0; i < popSz; ++i)
		{
			p = pop[i];

			if (!p) continue;

			parents = GAgetOriginalParents(i,iter);
			if (parents)
			{
				for (j=0; parents[j] > 0; ++j)
				{
					if (parents[j] >= max)
						max = parents[j];
					if (parents[j] >= num)
					{
						temp = ids;
						ids = (int*)malloc( num * 2 * sizeof(int) );
						
						for (k=0; k < num; ++k)
							ids[k] = temp[k];
						
						num *= 2;
						
						for (; k < num; ++k)
							ids[k] = 0;
							
						free(temp);
					}
					
					ids[ parents[j] ] += 1;
				}
			}
			else
			{
				/*j = getID(p);
				
				if (j < 0) continue;

				if (j >= max)
					max = j;
				if (j >= num)
				{
					temp = ids;
					ids = (int*)malloc( num * 2 * sizeof(int) );
					for (k=0; k < num; ++k)
						ids[k] = temp[k];
					
					num *= 2;
					
					for (; k < num; ++k)
						ids[k] = 0;
						
					free(temp);
				}
				
				ids[j]++;*/
			}
			
			if (PRINT_EACH_BEST_LINEAGE && !PRINT_EACH_ALL_LINEAGE)
				break;
		}
	}
	
	printf("%i",iter);
	fprintf(LOGFILE,"%i",iter);
	if (PRINT_EACH_FITNESS && !PRINT_EACH_ALL_FITNESS)
	{
		f = fitness(pop[0]);
		printf("\t%lf",f);
		fprintf(LOGFILE,"\t%lf",f);
	}
	
	if (PRINT_EACH_ALL_FITNESS)
	{
		for (i=0; i < popSz; ++i)
		{
			f = fitness(pop[i]);
			//printf("\t%lf",f);
			fprintf(LOGFILE,"\t%lf",f);
		}
	}
	
	if (PRINT_EACH_SIZE)
	{
		printf("\t%i\t%i",getNumSpecies(pop[0]),getNumReactions(pop[0]));
		fprintf(LOGFILE,"\t%i\t%i",getNumSpecies(pop[0]),getNumReactions(pop[0]));
	}
	
	if (TRACK_NETWORK_PARENTS && (PRINT_EACH_BEST_LINEAGE || PRINT_EACH_ALL_LINEAGE))
	{
		for (i=0; i < max; ++i)
		{
			//printf("\t%i",ids[i]);
			fprintf(LOGFILE,"\t%i",ids[i]);
		}
	}
	
	if (PRINT_EACH_SCRIPT)
	{
		//printf("\n========script=======\n");
		fprintf(LOGFILE,"\n========script=======\n");

		//printNetwork(stdout,pop[0]);
		printNetwork(LOGFILE,pop[0]);

		//printf("\n=====================\n");
		fprintf(LOGFILE,"\n=====================\n");
	}

	printf("\n");
	fprintf(LOGFILE,"\n");

	if (ids && (num > 0))
		free(ids);
	
	if (USER_CALLBACK_FNC)
		return USER_CALLBACK_FNC(iter,pop,popSz);
	
	return 0;
}

void finalCallBackWithLogKepping(int iter, GApopulation pop, int popSz)
{
	unsigned long long * seeds;
	
	//save
	int each_fitness = PRINT_EACH_FITNESS,
		each_script = PRINT_EACH_SCRIPT,
		each_size = PRINT_EACH_SIZE,
		each_best_lineage = PRINT_EACH_BEST_LINEAGE,
		each_all_fitness = PRINT_EACH_ALL_FITNESS,
		each_all_lineage = PRINT_EACH_ALL_LINEAGE;
	
	//cheat
	PRINT_EACH_FITNESS = PRINT_FINAL_FITNESS;
	PRINT_EACH_SCRIPT = PRINT_FINAL_SCRIPT;
	PRINT_EACH_SIZE = PRINT_FINAL_SIZE;
	PRINT_EACH_BEST_LINEAGE = PRINT_FINAL_BEST_LINEAGE;
	PRINT_EACH_ALL_FITNESS = PRINT_FINAL_ALL_FITNESS;
	PRINT_EACH_ALL_LINEAGE = PRINT_FINAL_ALL_LINEAGE;
	
	//printf("\n========final results=======\n");
	fprintf(LOGFILE,"\n========final results=======\n");
	callBackWithLogKeeping(iter,pop,popSz);
	
	//restore
	PRINT_EACH_FITNESS = each_fitness;
	PRINT_EACH_SCRIPT = each_script;
	PRINT_EACH_SIZE = each_size;
	PRINT_EACH_BEST_LINEAGE = each_best_lineage;
	PRINT_EACH_ALL_FITNESS = each_all_fitness;
	PRINT_EACH_ALL_LINEAGE = each_all_lineage;
	
	if (LOGFILE && PRINT_SEEDS)
	{
		seeds = getMTseeds();
		fprintf(LOGFILE,"random number generator seeds: %llf,%llf,%llf,%llf\n",seeds[0],seeds[1],seeds[2],seeds[3]);
	}
}

GApopulation evolveNetworks(int sz0,int sz1,int maxIter, GAFitnessFnc fitness, GACallbackFnc callbackFunc)
{
	GApopulation P;
	
	if (!fitness) return 0;
	
	GAsetFitnessFunction(fitness);
	GAsetCallbackFunction(callbackFunc);
	
	TRACK_NETWORK_PARENTS = GAisLineageTrackingOn();
	
	if (LOGFILE)
	{
		USER_CALLBACK_FNC = callbackFunc;
		callbackFunc = &callBackWithLogKeeping;
	}

	if (!MTrandHasBeenInitialized())
		initMTrand(); 

	P = randomNetworks(sz0);

	GAinit(&deleteNetwork, &cloneNetwork ,GAgetFitnessFunction(), &crossoverNetwork, &mutateNetwork, GAgetSelectionFunction(), callbackFunc);

	P = GArun(P,sz0,sz1,maxIter);
	
	if (LOGFILE)
	{
		USER_CALLBACK_FNC = 0;
		finalCallBackWithLogKepping(maxIter,P,sz1);
	}
	
	if (LOGFILE)
		fclose(LOGFILE);

	return P;
}


/*******************************
   Special fitness function
*******************************/

double compareSteadyStates(GAindividual p, double ** table, int rows, int inputs, int outputs, int corr, double ** res)
{
	int i, j, m, k, g, cols, n, *best;
	double * ss, * iv, closest, temp, sumOfSq, corrcoef, *mXY, *mX, *mY, *mX2, *mY2, oldMaxT = SS_FUNC_MAX_TIME;
	ReactionNetwork * r = (ReactionNetwork*)(p);
	SetFixedSpeciesFunction setFixed;
	
	SS_FUNC_MAX_TIME = 100.0;

	cols = inputs + outputs;
	
	n = getNumSpecies(r);
	
	if (n < cols) return 0.0; // not enough species
	
	setFixed = setFixedSpeciesFunctions[r->type];
	
	for (i=0; i < inputs; ++i)
	{
		setFixed(r->network,i,1);
	}

	sumOfSq = 0.0;
	corrcoef = 0.0;
	iv = (double*) malloc( n * sizeof(double) );
	for (i=0; i < n; ++i)
		iv[i] = r->initialValues[i];

	best = (int*) malloc ( outputs * sizeof(int) );
	mX = (double*)malloc( outputs * sizeof(double) );
	mY = (double*)malloc( outputs * sizeof(double) );
	mXY = (double*)malloc( outputs * sizeof(double) );
	mX2 = (double*)malloc( outputs * sizeof(double) );
	mY2 = (double*)malloc( outputs * sizeof(double) );

	for (i=0; i < outputs; ++i)
	{
		best[i] = -1;
		mXY[i] = mX[i] = mY[i] = mX2[i] = mY2[i] = 0;
	}

	for (m=0; m < rows; ++m)
	{
		for (i=0; i < inputs; ++i)
			r->initialValues[i] = table[m][i];
		
		for (i=inputs; i < n; ++i)
			r->initialValues[i] = 0.0;
		
		ss = networkSteadyState(r);

		if (ss) //error in simulation?
		{
			if (res)
			{
				for (i=0; i < inputs; ++i)
					res[m][i] = ss[i];
			}
			for (i=0; i < outputs; ++i) //for each target output
			{
				if (best[i] < 0)
				{
					closest = -1.0;
					for (j=inputs; j < n; ++j) //find best match
					{
						g = 0;
						for (k=0; k < outputs; ++k)
							if (best[k] == j)
							{
								g = 1;
								break;
							}
						if (g) continue;
						temp = (ss[j] - table[m][inputs+i]);
						if ((closest < 0.0) || ((temp*temp) < closest))
						{
							closest = temp*temp;
							best[i] = j;
						}
					}
				}
				
				j = inputs+i;//best[i];

				if (res)
				{
					res[m][inputs+i] = ss[j];
				}

				temp = (ss[j] - table[m][inputs+i]);
				closest = temp*temp;
				
				sumOfSq += closest;
				mX[i] += ss[j];
				mY[i] += table[m][inputs+i];
				mXY[i] += table[m][inputs+i] * ss[j];
				mX2[i] += ss[j]*ss[j];
				mY2[i] += table[m][inputs+i]*table[m][inputs+i];

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
	
	corrcoef = 0.0;
	for (i=0; i < outputs; ++i)
	{
		mX[i] /= rows;
		mY[i] /= rows;
		mXY[i] /= rows;
		mX2[i] /= rows;
		mY2[i] /= rows;

		if ( (mX2[i] - mX[i]*mX[i]) <= 0.0 || (mY2[i] - mY[i]*mY[i]) <= 0.0 )
		{
			sumOfSq = -1.0;
			break;
		}

		temp = ( (mXY[i] - mX[i]*mY[i])/
				( sqrt(mX2[i] - mX[i]*mX[i])*sqrt(mY2[i] - mY[i]*mY[i])) );   //correlation formula
		corrcoef += (1.0 + temp)/2.0; //between 0 and 1, instead of -1 and 1		
	}

	for (i=0; i < n; ++i)
		r->initialValues[i] = iv[i];
	
	free(iv);
	free(mX);
	free(mY);
	free(mXY);
	free(mX2);
	free(mY2);
	
	//restore the fixed
	for (i=0; i < inputs; ++i)
	{
		setFixed(r->network,i,0);
	}

	SS_FUNC_MAX_TIME = oldMaxT;
	
	if (sumOfSq <= 0.0) return 0.0;
	
	if (corr) return corrcoef;
	return (1.0 / (1.0 + sumOfSq));
}

void enableLogFile(char* filename)
{
	if (LOGFILE)
		fclose(LOGFILE);
			
	LOGFILE = 0;
	
	if (filename)
	{
		LOGFILE = fopen(filename, "w");
	}
}

void disableLogFile()
{
	enableLogFile(0);
}

void configureContinuousLog(int bestNetworkFitness, 
							int bestNetworkScript,
							int bestNetworkSize, 
							int bestNetworkLineage,
							int allFitness,
							int allNetworkLineage )
{
	
	PRINT_EACH_FITNESS = bestNetworkFitness > 0;
	PRINT_EACH_SCRIPT = bestNetworkScript > 0;
	PRINT_EACH_SIZE = bestNetworkSize > 0;
	PRINT_EACH_BEST_LINEAGE = bestNetworkLineage > 0;
	PRINT_EACH_ALL_FITNESS = allFitness > 0;
	PRINT_EACH_ALL_LINEAGE = allNetworkLineage > 0;
	if (PRINT_EACH_ALL_LINEAGE || PRINT_EACH_BEST_LINEAGE)
		TRACK_NETWORK_PARENTS = 1;
}

void configureFinalLog(int bestNetworkFitness, 
							int bestNetworkScript,
							int bestNetworkSize, 
							int bestNetworkLineage,
							int allFitness, 
							int allNetworkLineage,
							int seeds)
{
	PRINT_FINAL_FITNESS = bestNetworkFitness > 0;
	PRINT_FINAL_SCRIPT = bestNetworkScript > 0;
	PRINT_FINAL_SIZE = bestNetworkSize > 0;
	PRINT_FINAL_BEST_LINEAGE = bestNetworkLineage > 0;
	PRINT_FINAL_ALL_FITNESS = allFitness > 0;
	PRINT_FINAL_ALL_LINEAGE = allNetworkLineage > 0;
	PRINT_SEEDS = seeds > 0;
	if (PRINT_FINAL_BEST_LINEAGE || PRINT_FINAL_ALL_LINEAGE)
		TRACK_NETWORK_PARENTS = 1;
}

void setInitialValues( GAindividual x, double * values)
{
	int i,n;
	ReactionNetwork * rnet = (ReactionNetwork*)x;
	
	if (!rnet || !values) return;
	
	n = getNumSpecies(rnet);
	
	for (i=0; i < n; ++i)
		rnet->initialValues[i] = values[i];
}

double * getInitialValues( GAindividual x)
{
	int i,n;
	ReactionNetwork * rnet = (ReactionNetwork*)x;
	
	if (!rnet) return 0;
	
	return rnet->initialValues;
}

void lineageTrackingON()
{
	GAlineageTrackingON();
}

void lineageTrackingOFF()
{
	GAlineageTrackingOFF();
}

int isLineageTrackingOn()
{
	return GAisLineageTrackingOn();
}

int* getOriginalParents(int i, int j)
{
	return GAgetOriginalParents(i,j);
}

int* getImmediateParents(int i, int j)
{
	return GAgetImmediateParents(i,j);
}
