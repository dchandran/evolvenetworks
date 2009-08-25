/*******************************************************

	Copyright (C) 2009 Deepak Chandran
	see header file

********************************************************/

#include "geneRegulationNetwork.h"

/************************
    global variables
*************************/

static double VMAX_RANGE = 10.0;
static double KA_RANGE = 10.0;
static double TF_RANGE = 4.0;
static double DEG_RANGE = 2.0;
static int AVG_NUM_GENES = 5;
static int AVG_NUM_REGULATIONS = 2;

static double CROSSOVER_PROB = 1.0;
static double MUTATE_KA = 0.2;
static double MUTATE_DEGRADE = 0.2;
static double MUTATE_COMPLEX = 0.2;
static double ADD_GENE = 0.2;

void setParametersForGeneRegulationNetwork(double ka, double tfs, double vmax, double deg)
{
	KA_RANGE = ka;
	TF_RANGE = tfs;
	VMAX_RANGE = vmax;
	DEG_RANGE = deg;
}

void setSizeForGeneRegulationNetwork(int n1, int n2)
{
	AVG_NUM_GENES = n1;
	AVG_NUM_REGULATIONS = n2;
}


void setMutationAndCrossoverRatesForGeneRegulationNetwork(double ka, double degrade, double complex, double add, double remove, double crossover)
{
	double total;

	total = ka + degrade + complex + add + remove;
	
	if (crossover > 1.0) crossover /= 100.0;
	CROSSOVER_PROB = crossover;
	MUTATE_KA = ka/total;
	MUTATE_DEGRADE = degrade/total;
	MUTATE_COMPLEX = complex/total;
	ADD_GENE = add/total;
}

/*******************************************************
    Clone, delete, mutate, crossover  (required by GA)
*********************************************************/

void deleteGeneRegulationNetwork(void * individual)
{
	int i;
	GeneRegulationNetwork * net;
	
	if (!individual) return;
	
	net = (GeneRegulationNetwork*)(individual);
	
	if ((*net).complexes)
	{	
		for (i=0; i < (*net).numComplexes; ++i)
		{
			if ((*net).complexes[i].TFs) free ((*net).complexes[i].TFs);
		}
		free ((*net).complexes);
	}
	if ((*net).targetGene) free((*net).targetGene);
	if ((*net).Ka) free((*net).Ka);
	if ((*net).degradation) free((*net).degradation);
}

void* cloneGeneRegulationNetwork(void * individual)
{
	int i,j,m,n;
	GeneRegulationNetwork * net, * net2;
	
	if (!individual) return 0;
	
	net = (GeneRegulationNetwork*)(individual);   //original
	net2 = malloc(sizeof(GeneRegulationNetwork)); //soon to be clone
	
	m = (*net).species;    //number of genes
	n = (*net).numComplexes;    //number of complexes
	(*net2).species = m;
	(*net2).numComplexes = n;
	
	(*net2).complexes = malloc( n * sizeof (complex) );
	(*net2).targetGene = malloc( n * sizeof(int) );
	(*net2).Ka = malloc( n * sizeof(double) );
	(*net2).degradation = malloc( m * sizeof(double) );
	(*net2).Vmax = malloc( m * sizeof(double) );
	
	for (i=0; i < n; ++i)
	{
		(*net2).complexes[i].size = (*net).complexes[i].size;
		(*net2).complexes[i].TFs = malloc( (*net2).complexes[i].size * sizeof(int) );
		for (j=0; j < (*net2).complexes[i].size; ++j)
			(*net2).complexes[i].TFs[j] = (*net).complexes[i].TFs[j];
		
		(*net2).targetGene[i] = (*net).targetGene[i];
		(*net2).Ka[i] = (*net).Ka[i];
	}
	
	for (i=0; i < m; ++i)
	{
		(*net2).degradation[i] = (*net).degradation[i];
		(*net2).Vmax[i] = (*net).Vmax[i];
	}
	
	return (void*)(net2);  //done
}

void* crossoverGeneRegulationNetwork(void * individualA, void * individualB)  //crossover between complexes in two networks
{
	int i, j, k, i1, i2, n , m; 
	GeneRegulationNetwork * net1, * net2, * net3;
	
	if (mtrand() > CROSSOVER_PROB) return mutateGeneRegulationNetwork(cloneGeneRegulationNetwork(individualA));
	
	if (!individualA) return mutateGeneRegulationNetwork(cloneGeneRegulationNetwork(individualB));
	if (!individualB) return mutateGeneRegulationNetwork(cloneGeneRegulationNetwork(individualA));
	
	net1 = (GeneRegulationNetwork*)(individualA);  //parents
	net2 = (GeneRegulationNetwork*)(individualB);
	
	if ((*net1).numComplexes < 3) return mutateGeneRegulationNetwork(cloneGeneRegulationNetwork(net2));  //if parents are too small
	if ((*net2).numComplexes < 3) return mutateGeneRegulationNetwork(cloneGeneRegulationNetwork(net1));
		
	i1 = (int)(mtrand() * ((*net1).numComplexes - 1) + 1.0);	//crossover point in net1
	i2 = (int)(mtrand() * ((*net2).numComplexes - 2) + 1.0);	//crossover point in net2
	
	n = i1 + (*net2).numComplexes - i2;
	
	net3 = newGeneRegulationNetwork(m,n);  //child network
	
	for (i=0; i < i1; ++i)
	{
		(*net3).complexes[i].size = (*net1).complexes[i].size;
		(*net3).complexes[i].TFs = malloc( (*net1).complexes[i].size * sizeof(int) );
		for (j=0; j < (*net1).complexes[i].size; ++j)
		{
			(*net3).complexes[i].TFs[j] = (*net1).complexes[i].TFs[j];
			if ((*net3).complexes[i].TFs[j] > m)
				m = (*net3).complexes[i].TFs[j];
		}
		
		(*net3).targetGene[i] = (*net1).targetGene[i];
		if ((*net3).targetGene[i] > m)
			m = (*net3).targetGene[i];
		
		(*net3).Ka[i] = (*net1).Ka[i];
	}
	
	for (i=i2; i < (*net2).numComplexes; ++i)
	{
		k = i+i1-i2;
		(*net3).complexes[k].size = (*net2).complexes[i].size;
		(*net3).complexes[k].TFs = malloc( (*net2).complexes[i].size * sizeof(int) );
		for (j=0; j < (*net2).complexes[i].size; ++j)
		{
			(*net3).complexes[k].TFs[j] = (*net2).complexes[i].TFs[j];
			if ((*net3).complexes[k].TFs[j] > m)
				m = (*net3).complexes[k].TFs[j];
		}
		
		(*net3).targetGene[k] = (*net2).targetGene[i];
		if ((*net3).targetGene[k] > m)
			m = (*net3).targetGene[k];
		
		(*net3).Ka[k] = (*net2).Ka[i];
	}
	
	(*net3).species = m + 1;
	(*net3).degradation = malloc( (m+1) * sizeof(double) );
	(*net3).Vmax = malloc( (m+1) * sizeof(double) );
	
	for (i=0; i < (*net3).species; ++i)
	{
		if (i < (*net1).species)
		{
			(*net3).degradation[i] = (*net1).degradation[i];
			(*net3).Vmax[i] = (*net1).Vmax[i];
		}
		else
		if (i < (*net2).species)
		{
			(*net3).degradation[i] = (*net2).degradation[i];
			(*net3).Vmax[i] = (*net2).Vmax[i];
		}
		else
		{
			(*net3).degradation[i] = mtrand() * DEG_RANGE;
			(*net3).Vmax[i] = mtrand() * VMAX_RANGE;
		}
	}	
	
	return (void*)(net3);
}

void* mutateGeneRegulationNetwork(void * individual)
{
	int i,j,k,l,m,n;
	double r, * Ka, * degradation, * Vmax;
	int * targetGene;
	complex * complexes;
	GeneRegulationNetwork * net;
	
	net = (GeneRegulationNetwork*)individual;
	
	m = (*net).species;
	n = (*net).numComplexes;
	
	i = (int)(mtrand() * m);  //pick random reaction
	
	r = mtrand();

	if (r < MUTATE_KA)   //mutate ka
	{
		i = (int)(mtrand() * n);
		(*net).Ka[i] *= (4.0 * mtrand() - 2.0);
		return (void*)(net);
	}
	else
	if (r < (MUTATE_KA+MUTATE_DEGRADE)) //mutate degradation
	{
		i = (int)(mtrand() * m);
		if (mtrand() > 0.5)
			(*net).degradation[i] *= (mtrand() * 2.0);
		else
			(*net).Vmax[i] *= (mtrand() * 2.0);
		return (void*)(net);
	}
	else
	if (r < (MUTATE_KA+MUTATE_DEGRADE+MUTATE_COMPLEX))
	{
		i = (int)(mtrand() * n);
		j = (int)(mtrand() * (*net).complexes[i].size);
		(*net).complexes[i].TFs[j] = (int)(mtrand() * m);
		return (void*)(net);
	}
	else
	if (r < (MUTATE_KA+MUTATE_DEGRADE+MUTATE_COMPLEX+ADD_GENE))
	{
		++m;
		(*net).species = m;
		(*net).numComplexes = n+1;
		complexes = (*net).complexes;
		targetGene = (*net).targetGene;
		Ka = (*net).Ka;
		Vmax = (*net).Vmax;
		degradation = (*net).degradation;
		
		(*net).complexes = malloc( (1+n) * sizeof (complex) );
		(*net).targetGene = malloc( (1+n) * sizeof(int) );
		(*net).Ka = malloc( (1+n) * sizeof(double) );
		(*net).degradation = malloc( m * sizeof(double) );
		(*net).Vmax = malloc( m * sizeof(double) );
		
		for (j=0; j < n; ++j)
		{
			(*net).targetGene[j] = targetGene[j];
			(*net).complexes[j] = complexes[j];
			(*net).Ka[j] = Ka[j];
		}
		
		for (j=0; j < (m-1); ++j)
		{
			(*net).degradation[j] = degradation[j];
			(*net).Vmax[j] = Vmax[j];
		}
		
		(*net).complexes[n].size = (int)(1 + TF_RANGE * mtrand());
		(*net).complexes[n].TFs = malloc((*net).complexes[n].size * sizeof(int));
		(*net).complexes[n].TFs[0] = (int)(m-1);
		for (i=1; i < (*net).complexes[n].size; ++i)
		{
			(*net).complexes[n].TFs[i] = (int)(m * mtrand());
		}
		
		free(targetGene);
		free(complexes);
		free(Ka);
		free(degradation);
		free(Vmax);

		(*net).targetGene[n] = (int)(mtrand() * m);
		(*net).Ka[n] = (2.0 * mtrand() - 1.0) * KA_RANGE;
		(*net).degradation[m-1] = mtrand() * DEG_RANGE;
		(*net).Vmax[m-1] = mtrand() * VMAX_RANGE;
		
		return (void*)(net);
	}
	else 
	if (m > 2) //remove gene
	{	
		int t = 0, j1 = 0;
		
		i = (int)(mtrand() * m);
		
		--m;
		(*net).species = m;
		
		for (j=0; j < n; ++j)
		{
			for (k=0; k < (*net).complexes[j].size; ++k)
				if ((*net).complexes[j].TFs[k] == i)
				{
					++t;
					break;
				}
		}
		
		if (n < (t+2)) return (void*)net;
		
		complexes = (*net).complexes;
		targetGene = (*net).targetGene;
		Ka = (*net).Ka;
		degradation = (*net).degradation;
		Vmax = (*net).Vmax;
		
		(*net).numComplexes = n-t;
		(*net).complexes = malloc( (n-t) * sizeof (complex) );
		(*net).targetGene = malloc( (n-t) * sizeof(int) );
		(*net).Ka = malloc( (n-t) * sizeof(double) );
		(*net).degradation = malloc( (m) * sizeof(double) );
		(*net).Vmax = malloc( (m) * sizeof(double) );
		
		for (j=0, j1=0; j < (1+m) && j1 < m; ++j)
		{
			if (j != i)
			{
				(*net).degradation[j1] = degradation[j];
				(*net).Vmax[j1] = Vmax[j];
				++j1;
			}
		}
		
		for (j=0, j1 = 0; j < n; ++j)
		{
			l = 1;
			for (k=0; k < complexes[j].size; ++k)
				if (complexes[j].TFs[k] == i)
				{
					l = 0;
					break;
				}
			if (l)
			{
				(*net).complexes[j1].size = complexes[j].size;
				(*net).complexes[j1].TFs = complexes[j].TFs;
				(*net).Ka[j1] = Ka[j];
				(*net).targetGene[j1] = targetGene[j];
				++j1;
			}
			else
			{
				free(complexes[j].TFs);
			}
		}
		
		free(targetGene);
		free(complexes);
		free(Ka);
		free(degradation);
		free(Vmax);
		
		return (void*)(net);
	}
	
    return (void*)(net);
}

/*****************************************************
   Functions for simulating and printing
******************************************************/

void ratesForGeneRegulationNetwork(double time,double* u,double* rate,void * p)
{
	int i,j,k;
	double prod, num, denom;
	GeneRegulationNetwork * net;
	
	net = (GeneRegulationNetwork*)(p);
	
	for (i=0; i < (*net).species; ++i)
	{
		num = denom = 0;
		for (j=0; j < (*net).numComplexes; ++j)
		{
			if ((*net).targetGene[j] == i)
			{
				prod = 1.0;
				for (k=0; k < (*net).complexes[j].size; ++k)  //complex
					prod *= u[ (*net).complexes[j].TFs[k] ];
				
				if ((*net).Ka[j] > 0)
				{
					num += (*net).Ka[j] * prod;
					denom += (*net).Ka[j] * prod;
				}
				else
				{
					denom += -1.0 * (*net).Ka[j] * prod;
				}
			}
		}
		if (num == 0) num = 1.0;
		rate[i] = (*net).Vmax[i] *  num/(1+denom);
		rate[ i + (*net).species ] = (*net).degradation[i] * u[i];
	}
}


double * stoichiometryForGeneRegulationNetwork(void * p)
{
	int i,j,m,n;
	double prod, num, denom, *N;
	GeneRegulationNetwork * net;
	
	net = (GeneRegulationNetwork*)(p);
	
	m = (*net).species;
	n = 2 * m;
	N = malloc(m * n * sizeof(double));
	
	for (i=0; i < m; ++i)
		for (j=0; j < n; ++j)
			getValue(N,n,i,j) = 0.0;
			
	for (i=0; i < m; ++i)
	{
		getValue(N,n,i,i) = 1.0;
		getValue(N,n,i,i+m) = -1.0;
	}
	return N;
}

void printGeneRegulationNetwork(void * individual)
{
	int i,j,k,p;
	double prod, num, denom;
	GeneRegulationNetwork * net;
	
	if (!individual) return;
	
	net = (GeneRegulationNetwork*)(individual);
	
	for (i=0; i < (*net).species; ++i)
	{
		printf("$X -> s%i; %lf * (",i+1, (*net).Vmax[i]);
		num = denom = 0;
		p = 0;
		for (j=0; j < (*net).numComplexes; ++j)
		{
			if ((*net).targetGene[j] == i && (*net).Ka[j] > 0)
			{
				if (p > 0)
					printf("+");
				++p;
				printf("%lf",(*net).Ka[j]);
				for (k=0; k < (*net).complexes[j].size; ++k)
					printf("*s%i",(*net).complexes[j].TFs[k]+1);
			}
		}
		
		if (p < 1) 
			printf("1.0");
		
		printf(")/(1.0");
		
		p = 0;
		for (j=0; j < (*net).numComplexes; ++j)
		{
			if ((*net).targetGene[j] == i)
			{
				printf("+");
				++p;
				if ((*net).Ka[j] != 0)
				{
					if ((*net).Ka[j] < 0)
						printf("%lf",-(*net).Ka[j]);
					else
						printf("%lf",(*net).Ka[j]);
					for (k=0; k < (*net).complexes[j].size; ++k)
						printf("*s%i",(*net).complexes[j].TFs[k]+1);
				}
			}
		}
		
		printf(");\n");
		printf("s%i -> $X; %lf*s%i;\n",i+1,(*net).degradation[i],i+1);
	}
}

/***********************
  GA related functions
***********************/

GApopulation randomGeneRegulationNetworks(int num)
{
	int g,i,j,k,n,m;
	double r;
	GeneRegulationNetwork * net, **array;
	
	g = AVG_NUM_GENES;
	r = AVG_NUM_REGULATIONS;
	
	initMTrand(); /*initialize seeds for MT random number generator*/
	
	array = malloc( num * sizeof(GeneRegulationNetwork*) );
	
	for (k=0; k < num; ++k)
	{
		m = (int)(2 + g * 1.5 * (0.25 + mtrand()));
		n = (int)(2 + r * 1.5 * (0.25 + mtrand()));
		net = newGeneRegulationNetwork(m,n);
		
		for (i=0; i < n; ++i)
		{
			(*net).complexes[i].size = (int)(1 + TF_RANGE * mtrand());
			(*net).complexes[i].TFs = malloc((*net).complexes[i].size * sizeof(int));
			for (j=0; j < (*net).complexes[i].size; ++j)
			{
				(*net).complexes[i].TFs[j] = (int)(m * mtrand());
			}
			(*net).targetGene[i] = (int)(mtrand() * m);
			(*net).Ka[i] = (2.0 * mtrand() - 1.0) * KA_RANGE;
		}
		for (i=0; i < m; ++i)
		{
			(*net).degradation[i] = mtrand() * DEG_RANGE;
			(*net).Vmax[i] = mtrand() * VMAX_RANGE;
		}
		array[k] = net;
	}
	
	return (GApopulation)(array);
}

GeneRegulationNetwork * newGeneRegulationNetwork(int m,int n)
{
	int i;
	GeneRegulationNetwork * net;
	
	net = malloc(sizeof(GeneRegulationNetwork));
	(*net).species = m;    //number of genes
	(*net).numComplexes = n;    //number of complexes
	
	(*net).complexes = malloc( n * sizeof (complex) );
	(*net).targetGene = malloc( n * sizeof(int) );
	(*net).Ka = malloc( n * sizeof(double) );
	
	(*net).Vmax = malloc( m * sizeof(double) );
	(*net).degradation = malloc( m * sizeof(double) );
	
	for (i=0; i < n; ++i)
	{
		(*net).complexes[i].size = 0;
		(*net).complexes[i].TFs = 0;		
		(*net).targetGene[i] = 0;
		(*net).Ka[i] = 0.0;
	}
	for (i=0; i < m; ++i)
	{
		(*net).degradation[i] = 0.0;
		(*net).Vmax[i] = 0.0;
	}
	return net;
}

/*
int main()  //just for testing
{
	int i,j,n = 10;
	GApopulation pop = randomGeneRegulationNetworks(n,4,4);
	
	printf("generated networks\n");
	
	for (j=0; j < 50; ++j)
		for (i=0; i < n; ++i)
		{
			printGeneRegulationNetwork(pop[i]);
			pop[i] = mutateGeneRegulationNetwork(pop[i]);
			printf("\n");
			printGeneRegulationNetwork(pop[i]);
			printf("\n\n");
		}
	
	printf("mutation() test = ok\n\n");
	
	for (i=0; i < 50; ++i)
	{
		void * net = crossoverGeneRegulationNetworks(pop[ (int)(mtrand()*n) ],pop[ (int)(mtrand()*n) ]);
		printGeneRegulationNetwork(net);
		net = mutateGeneRegulationNetwork(net);
		deleteGeneRegulationNetwork(net);
	}
	
	printf("crossover() test = ok\n\n");
	
	printf("deleting...\n");
	
	for (i=0; i < n; ++i)
	{
		deleteGeneRegulationNetwork(pop[i]);
	}
	
	//testing simulation functions using ring oscillator
	
	GeneRegulationNetwork * net = newGeneRegulationNetwork(3,3);
	
	(*net).Vmax[0] = (*net).Vmax[1] = (*net).Vmax[2] = 2.0;
	(*net).degradation[0] = (*net).degradation[1] = (*net).degradation[2] = 0.1;
	
	(*net).complexes[0].size = 3;
	(*net).complexes[0].TFs = malloc(3*sizeof(int));
	(*net).complexes[0].TFs[0] = (*net).complexes[0].TFs[1] = (*net).complexes[0].TFs[2] = 0;
	(*net).Ka[0] = -1.0;
	(*net).targetGene[0] = 1;
	
	(*net).complexes[1].size = 3;
	(*net).complexes[1].TFs = malloc(3*sizeof(int));
	(*net).complexes[1].TFs[0] = (*net).complexes[1].TFs[1] = (*net).complexes[1].TFs[2] = 1;
	(*net).Ka[1] = -1.0;
	(*net).targetGene[1] = 2;
	
	(*net).complexes[2].size = 3;
	(*net).complexes[2].TFs = malloc(3*sizeof(int));
	(*net).complexes[2].TFs[0] = (*net).complexes[2].TFs[1] = (*net).complexes[2].TFs[2] = 2;
	(*net).Ka[2] = -1.0;
	(*net).targetGene[2] = 0;
	
	double x0[] = { 1.0, 1.0, 10.0 };
	
	double * N = stoichiometryForGeneRegulationNetwork((void*)net);
	int sz;
	
	double * y = SSA(3, 6, N, &(SSAfunction), x0, 0,500,100000,&sz,net); //simulate
	//double * y = ODEsim(3, x0, &(ODEfunction), 0, 500, 0.1, net);
	
	free(N);
	if (y)
	{
		writeToFile("temp.tab",y,sz,4);
		free(y);
	}
	
	deleteGeneRegulationNetwork((void*)net);
	
}
*/
