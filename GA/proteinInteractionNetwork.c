/*******************************************************

	Copyright (C) 2009 Deepak Chandran
	see header file

********************************************************/

#include "proteinInteractionNetwork.h"

/************************
    global variables
*************************/

static double PROB_SOURCE = 0.3;
static double VMAX_RANGE = 10.0;
static double KM_RANGE = 10.0;
static int COMPLEX_RANGE = 4;
static int AVG_NUM_PROTEINS = 5;
static int AVG_NUM_REACTIONS = 10;

static double CROSSOVER_PROB = 1.0;
static double MUTATE_KM = 0.2;
static double MUTATE_VMAX = 0.2;
static double MUTATE_COMPLEX = 0.2;
static double ADD_PROTEIN = 0.2;

void setParametersForProteinInteractionNetwork(int sz, double km, double vmax, double prob)
{
	KM_RANGE = km;
	COMPLEX_RANGE = sz;
	VMAX_RANGE = vmax;
	PROB_SOURCE = prob;
}

void setSizeForProteinInteractionNetwork(int n1, int n2)
{
	AVG_NUM_PROTEINS = n1;
	AVG_NUM_REACTIONS = n2;
}


void setMutationAndCrossoverRatesForProteinInteractionNetwork(double km, double vmax, double complex, double add, double remove, double crossover)
{
	double total;

	total = km + vmax + complex + add + remove;
	
	if (crossover > 1.0) crossover /= 100.0;
	CROSSOVER_PROB = crossover;
	MUTATE_KM = km/total;
	MUTATE_VMAX = vmax/total;
	MUTATE_COMPLEX = complex/total;
	ADD_PROTEIN = add/total;
}

/*******************************************************
    Clone, delete, mutate, crossover  (required by GA)
*********************************************************/

void deleteProteinInteractionNetwork(GAindividual individual)
{
	int i;
	ProteinInteractionNetwork * net;
	
	if (!individual) return;
	
	net = (ProteinInteractionNetwork*)(individual);
	
	if (net->complexes)
	{	
		for (i=0; i < net->reactions; ++i)
		{
			if (net->complexes[i].proteins) free (net->complexes[i].proteins);
		}
		free (net->complexes);
	}
	if (net->substrates) free(net->substrates);
	if (net->products) free(net->products);
	if (net->Km) free(net->Km);
	if (net->Vmax) free(net->Vmax);
	if (net->fixed) free (net->fixed);
}

GAindividual cloneProteinInteractionNetwork(GAindividual individual)
{
	int i,j,m,n;
	ProteinInteractionNetwork * net, * net2;
	
	if (!individual) return 0;
	
	net = (ProteinInteractionNetwork*)(individual);   //original
	net2 = (ProteinInteractionNetwork*) malloc(sizeof(ProteinInteractionNetwork)); //soon to be clone
	
	m = net->species;    //number of genes
	n = net->reactions;    //number of complexes
	net2->species = m;
	net2->reactions = n;
	
	net2->complexes = (ProteinComplex*) malloc( n * sizeof (ProteinComplex) );
	net2->substrates = (int*) malloc( n * sizeof(int) );
	net2->products = (int*) malloc( n * sizeof(int) );
	net2->Km = (double*) malloc( n * sizeof(double) );
	net2->Vmax = (double*) malloc( n * sizeof(double) );
	net2->fixed = (int*) malloc( m * sizeof(int) );
	
	for (i=0; i < n; ++i)
	{
		net2->complexes[i].size = net->complexes[i].size;
		net2->complexes[i].proteins = (int*) malloc( net2->complexes[i].size * sizeof(int) );
		for (j=0; j < net2->complexes[i].size; ++j)
			net2->complexes[i].proteins[j] = net->complexes[i].proteins[j];
		
		net2->substrates[i] = net->substrates[i];
		net2->products[i] = net->products[i];
		net2->Km[i] = net->Km[i];
		net2->Vmax[i] = net->Vmax[i];
	}
	
	for (i=0; i < m; ++i)
	{
		net2->fixed[i] = net->fixed[i];
	}
	
	return (GAindividual)(net2);  //done
}

GAindividual crossoverProteinInteractionNetwork(GAindividual individualA, GAindividual individualB)  //crossover between complexes in two networks
{
	int i, j, k, i1, i2, n , m; 
	ProteinInteractionNetwork * net1, * net2, * net3;
	
	if (mtrand() > CROSSOVER_PROB) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(individualA));
	
	if (!individualA) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(individualB));
	if (!individualB) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(individualA));
	
	net1 = (ProteinInteractionNetwork*)(individualA);  //parents
	net2 = (ProteinInteractionNetwork*)(individualB);
	
	if (net1->reactions < 3) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(net2));  //if parents are too small
	if (net2->reactions < 3) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(net1));
		
	i1 = (int)(mtrand() * (net1->reactions - 1) + 1.0);	//crossover point in net1
	i2 = (int)(mtrand() * (net2->reactions - 2) + 1.0);	//crossover point in net2
	
	n = i1 + net2->reactions - i2;
	
	m = 0;
	net3 = newProteinInteractionNetwork(m,n);  //child network
	
	for (i=0; i < i1; ++i)
	{
		net3->complexes[i].size = net1->complexes[i].size;
		net3->complexes[i].proteins = (int*) malloc( net1->complexes[i].size * sizeof(int) );
		for (j=0; j < net1->complexes[i].size; ++j)
		{
			net3->complexes[i].proteins[j] = net1->complexes[i].proteins[j];
			if (net3->complexes[i].proteins[j] > m)
				m = net3->complexes[i].proteins[j];
		}
		
		net3->substrates[i] = net1->substrates[i];
		net3->products[i] = net1->products[i];
		if (net3->substrates[i] > m) m = net3->substrates[i];
		if (net3->products[i] > m) m = net3->products[i];

		net3->Km[i] = net1->Km[i];
	}
	
	for (i=i2; i < net2->reactions; ++i)
	{
		k = i+i1-i2;
		net3->complexes[k].size = net2->complexes[i].size;
		net3->complexes[k].proteins = (int*) malloc( net2->complexes[i].size * sizeof(int) );
		for (j=0; j < net2->complexes[i].size; ++j)
		{
			net3->complexes[k].proteins[j] = net2->complexes[i].proteins[j];
			if (net3->complexes[k].proteins[j] > m)
				m = net3->complexes[k].proteins[j];
		}
		
		net3->substrates[k] = net2->substrates[i];
		net3->products[k] = net2->products[i];
		if (net3->substrates[k] > m) m = net3->substrates[k];
		if (net3->products[k] > m) m = net3->products[k];
		
		net3->Km[k] = net2->Km[i];
	}
	
	net3->species = m + 1;
	net3->fixed = (int*) malloc( (m+1) * sizeof(int) );
	
	for (i=0; i < net3->species; ++i)
	{
		if (i < net1->species)
		{
			net3->fixed[i] = net1->fixed[i];
		}
		else
		if (i < net2->species)
		{
			net3->fixed[i] = net2->fixed[i];
		}
		else
		{
			net3->fixed[i] = 0;
		}
	}
	
	return (GAindividual)(net3);
}

GAindividual mutateProteinInteractionNetwork(GAindividual individual)
{
	int i,j,k,l,m,n;
	double r, * Km, * Vmax;
	int * substrates, * products, * fixed;
	ProteinComplex * complexes;
	ProteinInteractionNetwork * net;
	
	net = (ProteinInteractionNetwork*)individual;
	
	m = net->species;
	n = net->reactions;
	
	i = (int)(mtrand() * m);  //pick random reaction
	
	r = mtrand();

	if (r < MUTATE_KM)   //mutate km
	{
		i = (int)(mtrand() * n);
		net->Km[i] *= (4.0 * mtrand() - 2.0);
		return (GAindividual)(net);
	}
	else
	if (r < (MUTATE_KM+MUTATE_VMAX)) //mutate vmax
	{
		i = (int)(mtrand() * m);
		net->Vmax[i] *= (mtrand() * 2.0);
		return (GAindividual)(net);
	}
	else
	if (r < (MUTATE_KM+MUTATE_VMAX+MUTATE_COMPLEX))
	{
		i = (int)(mtrand() * n);
		j = (int)(mtrand() * net->complexes[i].size);
		net->complexes[i].proteins[j] = (int)(mtrand() * m);
		return (GAindividual)(net);
	}
	else
	if (r < (MUTATE_KM+MUTATE_VMAX+MUTATE_COMPLEX+ADD_PROTEIN))
	{
		++m;
		net->species = m;
		net->reactions = n+1;
		complexes = net->complexes;
		substrates = net->substrates;
		products = net->products;
		Km = net->Km;
		Vmax = net->Vmax;
		fixed = net->fixed;
		
		net->complexes = (ProteinComplex*) malloc( (1+n) * sizeof (ProteinComplex) );
		net->substrates = (int*) malloc( (1+n) * sizeof(int) );
		net->products = (int*) malloc( (1+n) * sizeof(int) );
		net->Km = (double*) malloc( (1+n) * sizeof(double) );
		net->Vmax = (double*) malloc( (1+n) * sizeof(double) );
		net->fixed = (int*) malloc( m * sizeof(int) );
		
		for (j=0; j < n; ++j)
		{
			net->substrates[j] = substrates[j];
			net->products[j] = products[j];
			net->complexes[j] = complexes[j];
			net->Km[j] = Km[j];
			net->Vmax[j] = Vmax[j];
		}
		
		for (j=0; j < (m-1); ++j)
		{
			net->fixed[j] = fixed[j];
		}
		
		net->complexes[n].size = (int)(1 + COMPLEX_RANGE * mtrand());
		net->complexes[n].proteins = (int*) malloc(net->complexes[n].size * sizeof(int));
		net->complexes[n].proteins[0] = (int)(m-1);
		for (i=1; i < net->complexes[n].size; ++i)
		{
			net->complexes[n].proteins[i] = (int)(m * mtrand());
		}
		
		free(substrates);
		free(products);
		free(complexes);
		free(Km);
		free(Vmax);
		free(fixed);

		net->substrates[n] = (int)(mtrand() * m);
		net->products[n] = (int)(mtrand() * m);
		net->Km[n] = (2.0 * mtrand() - 1.0) * KM_RANGE;
		net->Vmax[n] = (2.0 * mtrand() - 1.0) * VMAX_RANGE;
		net->fixed[m-1] = 0;
		
		return (GAindividual)(net);
	}
	else 
	if (m > 2) //remove protein
	{	
		int t = 0, j1 = 0;
		
		i = (int)(mtrand() * m);
		
		--m;
		
		for (j=0; j < n; ++j)
		{
			for (k=0; k < net->complexes[j].size; ++k)
				if (net->complexes[j].proteins[k] == i)
				{
					++t;
					break;
				}
		}
		
		if (n < (t+2)) return (GAindividual)net;
		
		net->species = m;
		complexes = net->complexes;
		substrates = net->substrates;
		products = net->products;
		Km = net->Km;
		Vmax = net->Vmax;
		fixed = net->fixed;
		
		net->reactions = n-t;
		net->complexes = (ProteinComplex*) malloc( (n-t) * sizeof (ProteinComplex) );
		net->substrates = (int*) malloc( (n-t) * sizeof(int) );
		net->products = (int*) malloc( (n-t) * sizeof(int) );
		net->Km = (double*) malloc( (n-t) * sizeof(double) );
		net->Vmax = (double*) malloc( (n-t) * sizeof(double) );
		net->fixed = (int*) malloc( (m) * sizeof(int) );
		
		for (j=0, j1=0; j < (1+m) && j1 < m; ++j)
		{
			if (j != i)
			{
				net->fixed[j1] = fixed[j];
				++j1;
			}
		}
		
		for (j=0, j1 = 0; j < n; ++j)
		{
			l = 1;
			for (k=0; k < complexes[j].size; ++k)
				if (complexes[j].proteins[k] == i)
				{
					l = 0;
					break;
				}
			if (l)
			{
				net->complexes[j1].size = complexes[j].size;
				net->complexes[j1].proteins = complexes[j].proteins;
				net->Km[j1] = Km[j];
				net->Vmax[j1] = Vmax[j];
				net->substrates[j1] = substrates[j];
				net->products[j1] = products[j];
				++j1;
			}
			else
			{
				free(complexes[j].proteins);
			}
		}
		
		free(substrates);
		free(products);
		free(complexes);
		free(Km);
		free(Vmax);
		free(fixed);
		
		return (GAindividual)(net);
	}
	
    return (GAindividual)(net);
}

/*****************************************************
   Functions for simulating and printing
******************************************************/

int getNumSpeciesForProteinInteractionNetwork(GAindividual individual)
{
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);
	return (net->species);
}

int getNumReactionsForProteinInteractionNetwork(GAindividual individual)
{
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);
	return (net->reactions);
}

void setFixedSpeciesForProteinInteractionNetwork(GAindividual individual, int i, int value)
{
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);
	if (i < net->species)
		net->fixed[i] = value;
}

void ratesForProteinInteractionNetwork(double time,double* u,double* rate,GAindividual p)
{
	int j,k;
	double r,a,b,c,d;
	ProteinInteractionNetwork * net;
	
	net = (ProteinInteractionNetwork*)(p);
	
	for (j=0; j < net->reactions; ++j)
	{
		r = net->Vmax[j];
		a = net->Vmax[j];
		b = 1.0;
		for (k=0; k < net->complexes[j].size; ++k)  //complex
		{
			r *= u[ net->complexes[j].proteins[k] ];
			b *= u[ net->complexes[j].proteins[k] ];
		}
		
		if (net->substrates[j] >= 0)
		{
			rate[j] = r * u [ net->substrates[j] ]/(net->Km[j] + u [ net->substrates[j] ]);
			c  = u [ net->substrates[j] ]/(net->Km[j] + u [ net->substrates[j] ]);
		}
		else
			rate[j] = r;

		r = rate[j];
	}
}


double * stoichiometryForProteinInteractionNetwork(GAindividual p)
{
	int i,j;
	double *N;
	ProteinInteractionNetwork * net;
	
	net = (ProteinInteractionNetwork*)(p);
	N = (double*) malloc (net->species * net->reactions * sizeof(double));
	
	for (i=0; i < net->species; ++i)
		for (j=0; j < net->reactions; ++j)
			getValue(N,net->reactions,i,j) = 0.0;

	for (j=0; j < net->reactions; ++j)
	{
		if (net->substrates[j] >= 0)
			getValue(N,net->reactions, net->substrates[j] , j ) -= 1.0;
		
		if (net->products[j] >= 0)
			getValue(N,net->reactions, net->products[j] , j ) += 1.0;	
	}
	return N;
}

void printProteinInteractionNetwork(GAindividual p)
{
	int i,j,k,fix;
	ProteinInteractionNetwork * net;
	
	net = (ProteinInteractionNetwork*)(p);

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
		printf("const s%i",fix);
		for (i=0; i < net->species; ++i)
		{
			if (net->fixed[i])			
				printf(", s%i",i+1);			
		}
		printf("\n");
	}
	
	for (j=0; j < net->reactions; ++j)
	{
		if (net->substrates[j] >= 0)
			printf("s%i -> ",net->substrates[j]+1);
		else
			printf("$X -> ");

		if (net->products[j] >= 0)
			printf("s%i; ",net->products[j]+1);
		else
			printf("$X; ");

		printf("vmax%i",j);

		for (k=0; k < net->complexes[j].size; ++k)  //complex
			printf("*s%i",net->complexes[j].proteins[k]);
		
		if (net->substrates[j] >= 0)
			printf("/(km%i + s%i)",j,net->substrates[j]+1);

		printf(";\n");
	}
}

/***********************
  GA related functions
***********************/

GApopulation randomProteinInteractionNetworks(int num)
{
	int g,i,j,k,n,m;
	double r;
	ProteinInteractionNetwork * net, **array;
	
	g = AVG_NUM_PROTEINS;
	r = AVG_NUM_REACTIONS;
	
	initMTrand(); /*initialize seeds for MT random number generator*/
	
	array = (ProteinInteractionNetwork**) malloc( num * sizeof(ProteinInteractionNetwork*) );
	
	for (k=0; k < num; ++k)
	{
		m = (int)(2 + g * 1.5 * (0.25 + mtrand()));
		n = (int)(2 + r * 1.5 * (0.25 + mtrand()));
		net = newProteinInteractionNetwork(m,n);
		
		for (i=0; i < n; ++i)
		{
			net->complexes[i].size = (int)(1 + COMPLEX_RANGE * mtrand());
			net->complexes[i].proteins = (int*) malloc(net->complexes[i].size * sizeof(int));
			for (j=0; j < net->complexes[i].size; ++j)
			{
				net->complexes[i].proteins[j] = (int)(m * mtrand());
			}

			net->substrates[i] = (int)(mtrand() * m);
			net->products[i] = (int)(mtrand() * m);

			if (mtrand() < PROB_SOURCE) //source or sink
				if (mtrand() < 0.5)
					net->substrates[i] = -1;
				else
					net->products[i] = -1;

			net->Km[i] = mtrand() * KM_RANGE;
			net->Vmax[i] = mtrand() * VMAX_RANGE;
		}
		array[k] = net;
	}
	
	return (GApopulation)(array);
}

ProteinInteractionNetwork * newProteinInteractionNetwork(int m,int n)
{
	int i;
	ProteinInteractionNetwork * net;
	
	net = (ProteinInteractionNetwork*) malloc(sizeof(ProteinInteractionNetwork));
	net->species = m;    //number of genes
	net->reactions = n;    //number of complexes
	
	net->complexes = (ProteinComplex*) malloc( n * sizeof (ProteinComplex) );
	net->substrates = (int*) malloc( n * sizeof(int) );
	net->products = (int*) malloc( n * sizeof(int) );
	net->Km = (double*) malloc( n * sizeof(double) );	
	net->Vmax = (double*) malloc( n * sizeof(double) );

	for (i=0; i < n; ++i)
	{
		net->complexes[i].size = 0;
		net->complexes[i].proteins = 0;		
		net->substrates[i] = 0;
		net->products[i] = 0;
		net->Km[i] = 0.0;
		net->Vmax[i] = 0.0;
	}
	
	net->fixed = 0;

	if (m > 0)
	{
		net->fixed = (int*) malloc( m * sizeof(int) );
		for (i=0; i < m; ++i)
		{
			net->fixed[i] = 0;
		}
	}
	return net;
}

/*
int main()  //just for testing
{
	int i,j,n = 10;
	GApopulation pop = randomProteinInteractionNetworks(n,4,4);
	
	printf("generated networks\n");
	
	for (j=0; j < 50; ++j)
		for (i=0; i < n; ++i)
		{
			printProteinInteractionNetwork(pop[i]);
			pop[i] = mutateProteinInteractionNetwork(pop[i]);
			printf("\n");
			printProteinInteractionNetwork(pop[i]);
			printf("\n\n");
		}
	
	printf("mutation() test = ok\n\n");
	
	for (i=0; i < 50; ++i)
	{
		GAindividual net = crossoverProteinInteractionNetworks(pop[ (int)(mtrand()*n) ],pop[ (int)(mtrand()*n) ]);
		printProteinInteractionNetwork(net);
		net = mutateProteinInteractionNetwork(net);
		deleteProteinInteractionNetwork(net);
	}
	
	printf("crossover() test = ok\n\n");
	
	printf("deleting...\n");
	
	for (i=0; i < n; ++i)
	{
		deleteProteinInteractionNetwork(pop[i]);
	}
	
	//testing simulation functions using ring oscillator
	
	ProteinInteractionNetwork * net = newProteinInteractionNetwork(3,3);
	
	net->Vmax[0] = net->Vmax[1] = net->Vmax[2] = 2.0;
	net->degradation[0] = net->degradation[1] = net->degradation[2] = 0.1;
	
	net->complexes[0].size = 3;
	net->complexes[0].proteins = malloc(3*sizeof(int));
	net->complexes[0].proteins[0] = net->complexes[0].proteins[1] = net->complexes[0].proteins[2] = 0;
	net->Km[0] = -1.0;
	net->targetGene[0] = 1;
	
	net->complexes[1].size = 3;
	net->complexes[1].proteins = malloc(3*sizeof(int));
	net->complexes[1].proteins[0] = net->complexes[1].proteins[1] = net->complexes[1].proteins[2] = 1;
	net->Km[1] = -1.0;
	net->targetGene[1] = 2;
	
	net->complexes[2].size = 3;
	net->complexes[2].proteins = malloc(3*sizeof(int));
	net->complexes[2].proteins[0] = net->complexes[2].proteins[1] = net->complexes[2].proteins[2] = 2;
	net->Km[2] = -1.0;
	net->targetGene[2] = 0;
	
	double x0[] = { 1.0, 1.0, 10.0 };
	
	double * N = stoichiometryForProteinInteractionNetwork((GAindividual)net);
	int sz;
	
	double * y = SSA(3, 6, N, &(SSAfunction), x0, 0,500,100000,&sz,net); //simulate
	//double * y = ODEsim(3, x0, &(ODEfunction), 0, 500, 0.1, net);
	
	free(N);
	if (y)
	{
		writeToFile("temp.tab",y,sz,4);
		free(y);
	}
	
	deleteProteinInteractionNetwork((GAindividual)net);
	
}
*/
