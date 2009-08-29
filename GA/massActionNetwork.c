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
static double MUTATE_REMOVE_REACTION = 0.25;
static double PROB_UNI_UNI = 0.2;
static double PROB_UNI_BI = 0.2;
static double PROB_BI_UNI = 0.2;
static double PROB_BI_BI = 0.2;
static double PROB_NO_REACTANTS = 0.2;
static double PROB_NO_PRODUCTS = 0.2;
static double AVG_RATE_CONSTANT = 1.0;
static int AVG_NUM_SPECIES = 1;
static int AVG_NUM_REACTIONS = 1;

void setParametersForMassActionNetwork(double uni_uni, double uni_bi, double bi_uni, double bi_bi, double no_reactant, double no_product, double avg_rate_constant)
{
	double total = uni_uni + uni_bi + bi_uni + bi_bi + no_reactant + no_product;
	if (avg_rate_constant > 0)
		AVG_RATE_CONSTANT = avg_rate_constant;
	
	if (total <= 0) return;

	PROB_UNI_UNI = uni_uni/total;
	PROB_UNI_BI = uni_bi/total;
	PROB_BI_UNI = bi_uni/total;
	PROB_BI_BI = bi_bi/total;
	PROB_NO_REACTANTS = no_reactant/total;
	PROB_NO_PRODUCTS = no_product/total;
}

void setSizeForMassActionNetwork(int n, int s)
{
	AVG_NUM_SPECIES = n;
	AVG_NUM_REACTIONS = s;
}

void setMutationAndCrossoverRatesForMassActionNetwork(double mutatek, double remove, double add, double crossover)
{
	double total;

	if (mutatek < 0) mutatek = 0;
	if (remove < 0) remove = 0;
	if (add < 0) add = 0;
		
	total = mutatek + remove + add;
	if (total == 0) 
		mutatek = remove = 0.33333;
	else
	{
		mutatek /= total;
		remove /= total;
	}
	if (crossover > 1.0) crossover /= 100.0;
	MUTATE_COEFF_PROB = mutatek;
	MUTATE_REMOVE_REACTION = remove;
	CROSSOVER_PROB = crossover;
}

/*******************************************************
    Clone, delete, mutate, crossover  (required by GA)
*********************************************************/

void deleteMassActionNetwork(GAindividual individual)
{
	MassActionNetwork * net;
	
	if (!individual) return;
	net = (MassActionNetwork*)(individual);
	
	if (net->r1) free(net->r1);  //free each array
	if (net->r2) free(net->r2);
	if (net->p1) free(net->p1);
	if (net->p2) free(net->p2);
	if (net->k) free(net->k);
	if (net->fixed) free(net->fixed);
}

GAindividual cloneMassActionNetwork(GAindividual individual)
{
	int i,m,n;
	MassActionNetwork * net, * net2;
	
	if (!individual) return 0;
	
	net = (MassActionNetwork*)(individual);   //original
	net2 = (MassActionNetwork*) malloc(sizeof(MassActionNetwork)); //soon to be clone
	
	m = net->reactions;    //number of reactions
	n = net->species;    //number of species
	net2->reactions = m;
	net2->species = n;
	
	net2->k = (double*) malloc(m * sizeof(double));   //allocate space
	net2->r1 = (int*) malloc(m * sizeof(int));
	net2->r2 = (int*) malloc(m * sizeof(int));
	net2->p1 = (int*) malloc(m * sizeof(int));
	net2->p2 = (int*) malloc(m * sizeof(int));
	net2->fixed = (int*) malloc(n * sizeof(int));
	
	for (i=0; i < n; ++i)   //copy values
	{
		net2->fixed[i] = net->fixed[i];
	}
	
	for (i=0; i < m; ++i)   //copy values
	{
		net2->k[i] = net->k[i];
		net2->r1[i] = net->r1[i];
		net2->r2[i] = net->r2[i];
		net2->p1[i] = net->p1[i]; 
		net2->p2[i] = net->p2[i];
	}
	
	return (GAindividual)(net2);  //done
}

GAindividual crossoverMassActionNetwork(GAindividual individualA, GAindividual individualB)
{
	int i, i1, i2, n = 0, m = 0; 
	MassActionNetwork * net1, * net2, * net3;
	
	if (mtrand() > CROSSOVER_PROB) return mutateMassActionNetwork(cloneMassActionNetwork(individualA));  //do crossover?
	
	if (!individualA) return mutateMassActionNetwork(cloneMassActionNetwork(individualB));
	if (!individualB) return mutateMassActionNetwork(cloneMassActionNetwork(individualA));
	
	net1 = (MassActionNetwork*)(individualA);  //parents
	net2 = (MassActionNetwork*)(individualB);
	
	if (net1->reactions < 3) return mutateMassActionNetwork(cloneMassActionNetwork(net2));  //if parents are too small
	if (net2->reactions < 3) return mutateMassActionNetwork(cloneMassActionNetwork(net1));
	
	i1 = (int)(mtrand() * (net1->reactions - 1) + 1.0);	//crossover point in net1
	i2 = (int)(mtrand() * (net2->reactions - 2) + 1.0);	//crossover point in net2
	
	n = i1 + net2->reactions - i2;
	
	net3 = newMassActionNetwork(m,n);  //child network
	
	for (i=0; i < i1; ++i)
	{
		net3->k[i] = net1->k[i];
		net3->r1[i] = net1->r1[i];
		net3->r2[i] = net1->r2[i];
		net3->p1[i] = net1->p1[i]; 
		net3->p2[i] = net1->p2[i];
		
		if (net1->r1[i] > m) m = net1->r1[i];
		if (net1->r2[i] > m) m = net1->r2[i];
		if (net1->p1[i] > m) m = net1->p1[i];
		if (net1->p2[i] > m) m = net1->p2[i];
	}
	
	for (i=i2; i < net2->reactions; ++i)
	{
		net3->k[i+i1-i2] = net2->k[i];
		net3->r1[i+i1-i2] = net2->r1[i];
		net3->r2[i+i1-i2] = net2->r2[i];
		net3->p1[i+i1-i2] = net2->p1[i]; 
		net3->p2[i+i1-i2] = net2->p2[i];
		
		if (net2->r1[i] > m) m = net2->r1[i];
		if (net2->r2[i] > m) m = net2->r2[i];
		if (net2->p1[i] > m) m = net2->p1[i];
		if (net2->p2[i] > m) m = net2->p2[i];
	}
	
	net3->species = m + 1;
	
	net3->fixed = (int*) malloc( (m + 1)*sizeof(int) );
	for (i=0; i < (m+1); ++i)
		net3->fixed[i] = 0;
	
	return (GAindividual)(net3);
}

GAindividual mutateMassActionNetwork(GAindividual individual)
{
	int i,j,j2,m,n;
	double r;
	MassActionNetwork * net, * net2;
	
	net = (MassActionNetwork*)individual;

	if (!net) return 0;

	m = net->reactions;
	n = net->species;

	i = (int)(mtrand() * m);  //pick random reaction

	if (mtrand() < MUTATE_COEFF_PROB)   //mutate coefficient
	{
		net->k[i] *= (mtrand() * 2.0);
		return (GAindividual)(net);
	}
	else              //add or remove a new reaction and/or species to the network
	{
		if (mtrand() < (MUTATE_COEFF_PROB/2.0 + MUTATE_REMOVE_REACTION) && m > 2)     //remove a reaction
		{
			net2 = newMassActionNetwork( n, (m - 1) );
			for (j=0,j2=0; j < m; ++j)
			{
				if (j != i)
				{
					net2->k[j2] = net->k[j];
					net2->r1[j2] = net->r1[j];
					net2->r2[j2] = net->r2[j];
					net2->p1[j2] = net->p1[j];
					net2->p2[j2] = net->p2[j];
					++j2;
				}
			}
			deleteMassActionNetwork(net);
			return (GAindividual)(net2);
		}
		else
		{
			net2 = newMassActionNetwork( n, (m + 1) );
			for (j=0; j < m; ++j)
			{
				net2->k[j] = net->k[j];
				net2->r1[j] = net->r1[j];
				net2->r2[j] = net->r2[j];
				net2->p1[j] = net->p1[j];
				net2->p2[j] = net->p2[j];
			}
			net2->k[j] = AVG_RATE_CONSTANT * mtrand();   //reaction rate constant
			
			r = mtrand();

			if (r < PROB_UNI_UNI)
			{
				net2->r1[j] = (int)(net2->species * mtrand());  //first reactant
				net2->p1[j] = (int)(net2->species * mtrand());  //first product
			}
			else
			if (r < (PROB_UNI_BI + PROB_UNI_UNI))
			{
				net2->r1[j] = (int)(net2->species * mtrand());  //first reactant
				net2->p1[j] = (int)(net2->species * mtrand());  //first product
				net2->p2[j] = (int)(net2->species * mtrand());  //second product
			}
			else
			if (r < (PROB_BI_UNI + PROB_UNI_BI + PROB_UNI_UNI))
			{
				net2->r1[j] = (int)(net2->species * mtrand());  //first reactant
				net2->r2[j] = (int)(net2->species * mtrand());  //second reactant
				net2->p1[j] = (int)(net2->species * mtrand());  //second product
			}
			else
			{
				net2->r1[j] = (int)(net2->species * mtrand());  //first reactant
				net2->r2[j] = (int)(net2->species * mtrand());  //second reactant
				net2->p1[j] = (int)(net2->species * mtrand());  //first product
				net2->p2[j] = (int)(net2->species * mtrand());  //second product
			}

			r = mtrand();

			if (r < PROB_NO_REACTANTS)
				net2->r1[j] = net2->r2[j] = -1;
			else
			if (r < PROB_NO_PRODUCTS)
				net2->p1[j] = net2->p2[j] = -1;

			deleteMassActionNetwork(net);
			return (GAindividual)(net2);
		}
	}
    return (net);
}

/*****************************************************
   Functions for simulating and printing
******************************************************/


void ratesForMassActionNetwork(double time,double* u,double* rate,GAindividual p)
{
	int i;
	MassActionNetwork * net;
	
	net = (MassActionNetwork*)(p);

	if (!net) return;
	
	for (i=0; i < net->reactions; ++i)
	{
		rate[i] = net->k[i];
		if (net->r1[i] > -1) rate[i] *= u[ net->r1[i] ];
		if (net->r2[i] > -1) rate[i] *= u[ net->r2[i] ];
	}
}

double * stoichiometryForMassActionNetwork(GAindividual p)
{
	int i,j,n;
	double * N;
	MassActionNetwork * net;
	
	net = (MassActionNetwork*)(p);
	
	n = net->reactions;
	N = (double*) malloc(net->species * n * sizeof(double));
	for (i=0; i < n; ++i)
	{
		for (j=0; j < net->species; ++j)
			getValue(N,n,j,i) = 0;
		
		if ((net->r1[i] >= 0) && ((net->fixed[ net->r1[i] ]) == 0)) 
			getValue(N,n,net->r1[i],i) -= 1;
		
		if ((net->r2[i] >= 0) && ((net->fixed[ net->r2[i] ]) == 0)) 
			getValue(N,n,net->r2[i],i) -= 1;
			
		if ((net->p1[i] >= 0) && ((net->fixed[ net->p1[i] ]) == 0)) 
			getValue(N,n,net->p1[i],i) += 1;
		
		if ((net->p2[i] >= 0) && ((net->fixed[ net->p2[i] ]) == 0)) 
			getValue(N,n,net->p2[i],i) += 1;
	}
	return N;
}


int getNumSpeciesForMassActionNetwork(GAindividual individual)
{
	MassActionNetwork * net = (MassActionNetwork*)(individual);
	if (!net) return 0;
	return (net->species);
}

int getNumReactionsForMassActionNetwork(GAindividual individual)
{
	MassActionNetwork * net = (MassActionNetwork*)(individual);
	if (!net) return 0;
	return (net->reactions);
}

void setFixedSpeciesForMassActionNetwork(GAindividual individual, int i, int value)
{
	MassActionNetwork * net = (MassActionNetwork*)(individual);
	if (i < net->species)
		net->fixed[i] = value;
}


void printMassActionNetwork(FILE *stream, GAindividual individual)
{
	int i,fix;
	MassActionNetwork * net;
	
	if (!individual) return;
	net = (MassActionNetwork*)individual;
	
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
					fprintf(stream, "k%i * s%i",i+1,net->r1[i]+1);
				else
					fprintf(stream, "k%i * s%i",i+1,net->r2[i]+1);
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
	double u;
	MassActionNetwork * net;
	MassActionNetwork ** array;
	
	initMTrand(); /*initialize seeds for MT random number generator*/
	array = (MassActionNetwork**)malloc(num * sizeof(MassActionNetwork*));
	for (i=0; i < num; ++i)
	{
		n = (int)(1 + s * 2.0 * mtrand());
		net = newMassActionNetwork(n,(int)(2 + r * 2.0 * mtrand()));
		
		for (j=0; j < net->reactions; ++j)
		{
			net->k[j] = AVG_RATE_CONSTANT * 2.0 * mtrand();   //reaction rate constant
			net->r1[j] = net->r2[j] = net->p1[j] = net->p2[j] = -1;

			u = mtrand();

			if (u < PROB_UNI_UNI)
			{
				net->r1[j] = (int)(net->species * mtrand());  //first reactant
				net->p1[j] = (int)(net->species * mtrand());  //first product
			}
			else
			if (u < (PROB_UNI_BI + PROB_UNI_UNI))
			{
				net->r1[j] = (int)(net->species * mtrand());  //first reactant
				net->p1[j] = (int)(net->species * mtrand());  //first product
				net->p2[j] = (int)(net->species * mtrand());  //second product
			}
			else
			if (u < (PROB_BI_UNI + PROB_UNI_BI + PROB_UNI_UNI))
			{
				net->r1[j] = (int)(net->species * mtrand());  //first reactant
				net->r2[j] = (int)(net->species * mtrand());  //second reactant
				net->p1[j] = (int)(net->species * mtrand());  //second product
			}
			else
			{
				net->r1[j] = (int)(net->species * mtrand());  //first reactant
				net->r2[j] = (int)(net->species * mtrand());  //second reactant
				net->p1[j] = (int)(net->species * mtrand());  //first product
				net->p2[j] = (int)(net->species * mtrand());  //second product
			}
			
			u = mtrand();

			if (u < PROB_NO_REACTANTS)
				net->r1[j] = net->r2[j] = -1;
			else
			if (u < PROB_NO_PRODUCTS)
				net->p1[j] = net->p2[j] = -1;
		}
		
		array[i] = net;
	}
	
	return (GApopulation)(array);
}

MassActionNetwork * newMassActionNetwork(int s,int r)
{
	int i;
	MassActionNetwork * net;

	net = (MassActionNetwork*)malloc(sizeof(MassActionNetwork));
	net->species = s;
	net->reactions = r;
	
	net->fixed = (int*)malloc(s * sizeof(int));
	for (i=0; i < s; ++i)
		net->fixed[i] = 0; //no fixed species by default
	
	net->k = (double*) malloc(net->reactions * sizeof(double));
	net->r1 = (int*) malloc(net->reactions * sizeof(int));
	net->r2 = (int*) malloc(net->reactions * sizeof(int));
	net->p1 = (int*) malloc(net->reactions * sizeof(int));
	net->p2 = (int*) malloc(net->reactions * sizeof(int));
	
	for (i=0; i < r; ++i)
	{
		net->r1[i] = net->r2[i] = net->p1[i] = net->p2[i] = -1;
		net->k[i] = 1.0;
	}

	return net;
}
