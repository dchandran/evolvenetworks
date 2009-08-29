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
static double DEG_RANGE = 2.0;
static int TF_RANGE = 4;
static int AVG_NUM_GENES = 5;
static int AVG_NUM_REGULATIONS = 2;

static double CROSSOVER_PROB = 1.0;
static double MUTATE_KA = 0.2;
static double MUTATE_PRODUCTION = 0.2;
static double MUTATE_COMPLEX = 0.2;
static double ADD_GENE = 0.2;

void setParametersForGeneRegulationNetwork(int tfs, double ka, double vmax, double deg)
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
	MUTATE_PRODUCTION = degrade/total;
	MUTATE_COMPLEX = complex/total;
	ADD_GENE = add/total;
}

/*******************************************************
    Clone, delete, mutate, crossover  (required by GA)
*********************************************************/

void deleteGeneRegulationNetwork(GAindividual individual)
{
	int i;
	GeneRegulationNetwork * net;
	
	if (!individual) return;
	
	net = (GeneRegulationNetwork*)(individual);
	
	if (net->complexes)
	{	
		for (i=0; i < net->numComplexes; ++i)
		{
			if (net->complexes[i].TFs) free (net->complexes[i].TFs);
		}
		free (net->complexes);
	}
	if (net->targetGene) free(net->targetGene);
	if (net->Ka) free(net->Ka);
	if (net->degradation) free(net->degradation);
	if (net->fixed)
		free (net->fixed);
}

GAindividual cloneGeneRegulationNetwork(GAindividual individual)
{
	int i,j,m,n;
	GeneRegulationNetwork * net, * net2;
	
	if (!individual) return 0;
	
	net = (GeneRegulationNetwork*)(individual);   //original
	net2 = (GeneRegulationNetwork*) malloc(sizeof(GeneRegulationNetwork)); //soon to be clone
	
	m = net->species;    //number of genes
	n = net->numComplexes;    //number of complexes
	net2->species = m;
	net2->numComplexes = n;
	
	net2->complexes = (TFComplex*) malloc( n * sizeof (TFComplex) );
	net2->targetGene = (int*) malloc( n * sizeof(int) );
	net2->Ka = (double*) malloc( n * sizeof(double) );
	net2->degradation = (double*) malloc( m * sizeof(double) );
	net2->Vmax = (double*) malloc( m * sizeof(double) );
	net2->fixed = (int*) malloc( m * sizeof(int) );
	
	for (i=0; i < n; ++i)
	{
		net2->complexes[i].size = net->complexes[i].size;
		net2->complexes[i].TFs = (int*) malloc( net2->complexes[i].size * sizeof(int) );
		for (j=0; j < net2->complexes[i].size; ++j)
			net2->complexes[i].TFs[j] = net->complexes[i].TFs[j];
		
		net2->targetGene[i] = net->targetGene[i];
		net2->Ka[i] = net->Ka[i];
	}
	
	for (i=0; i < m; ++i)
	{
		net2->fixed[i] = net->fixed[i];
		net2->degradation[i] = net->degradation[i];
		net2->Vmax[i] = net->Vmax[i];
	}
	
	return (GAindividual)(net2);  //done
}

GAindividual crossoverGeneRegulationNetwork(GAindividual individualA, GAindividual individualB)  //crossover between complexes in two networks
{
	int i, j, k, i1, i2, n , m; 
	GeneRegulationNetwork * net1, * net2, * net3;
	
	if (mtrand() > CROSSOVER_PROB) return mutateGeneRegulationNetwork(cloneGeneRegulationNetwork(individualA));
	
	if (!individualA) return mutateGeneRegulationNetwork(cloneGeneRegulationNetwork(individualB));
	if (!individualB) return mutateGeneRegulationNetwork(cloneGeneRegulationNetwork(individualA));
	
	net1 = (GeneRegulationNetwork*)(individualA);  //parents
	net2 = (GeneRegulationNetwork*)(individualB);
	
	if (net1->numComplexes < 3) return mutateGeneRegulationNetwork(cloneGeneRegulationNetwork(net2));  //if parents are too small
	if (net2->numComplexes < 3) return mutateGeneRegulationNetwork(cloneGeneRegulationNetwork(net1));
		
	i1 = (int)(mtrand() * (net1->numComplexes - 1) + 1.0);	//crossover point in net1
	i2 = (int)(mtrand() * (net2->numComplexes - 2) + 1.0);	//crossover point in net2
	
	n = i1 + net2->numComplexes - i2;
	
	m = 0;
	net3 = newGeneRegulationNetwork(m,n);  //child network
	
	for (i=0; i < i1; ++i)
	{
		net3->complexes[i].size = net1->complexes[i].size;
		net3->complexes[i].TFs = (int*) malloc( net1->complexes[i].size * sizeof(int) );
		for (j=0; j < net1->complexes[i].size; ++j)
		{
			net3->complexes[i].TFs[j] = net1->complexes[i].TFs[j];
			if (net3->complexes[i].TFs[j] > m)
				m = net3->complexes[i].TFs[j];
		}
		
		net3->targetGene[i] = net1->targetGene[i];
		if (net3->targetGene[i] > m)
			m = net3->targetGene[i];
		
		net3->Ka[i] = net1->Ka[i];
	}
	
	for (i=i2; i < net2->numComplexes; ++i)
	{
		k = i+i1-i2;
		net3->complexes[k].size = net2->complexes[i].size;
		net3->complexes[k].TFs = (int*) malloc( net2->complexes[i].size * sizeof(int) );
		for (j=0; j < net2->complexes[i].size; ++j)
		{
			net3->complexes[k].TFs[j] = net2->complexes[i].TFs[j];
			if (net3->complexes[k].TFs[j] > m)
				m = net3->complexes[k].TFs[j];
		}
		
		net3->targetGene[k] = net2->targetGene[i];
		if (net3->targetGene[k] > m)
			m = net3->targetGene[k];
		
		net3->Ka[k] = net2->Ka[i];
	}
	
	net3->species = m + 1;
	
	net3->degradation = (double*) malloc( (m+1) * sizeof(double) );
	net3->Vmax = (double*) malloc( (m+1) * sizeof(double) );
	net3->fixed = (int*) malloc( (m+1) * sizeof(int) );
	
	for (i=0; i < net3->species; ++i)
	{
		if (i < net1->species)
		{
			net3->degradation[i] = net1->degradation[i];
			net3->Vmax[i] = net1->Vmax[i];
			net3->fixed[i] = net1->fixed[i];
		}
		else
		if (i < net2->species)
		{
			net3->degradation[i] = net2->degradation[i];
			net3->Vmax[i] = net2->Vmax[i];
			net3->fixed[i] = net2->fixed[i];
		}
		else
		{
			net3->degradation[i] = mtrand() * DEG_RANGE;
			net3->Vmax[i] = mtrand() * VMAX_RANGE;
			net3->fixed[i] = 0;
		}
	}	
	
	return (GAindividual)(net3);
}

GAindividual mutateGeneRegulationNetwork(GAindividual individual)
{
	int i,j,k,l,m,n;
	double r, * Ka, * degradation, * Vmax;
	int * targetGene, * fixed;
	TFComplex * complexes;
	GeneRegulationNetwork * net;
	
	net = (GeneRegulationNetwork*)individual;
	
	m = net->species;
	n = net->numComplexes;
	
	i = (int)(mtrand() * m);  //pick random reaction
	
	r = mtrand();

	if (r < MUTATE_KA)   //mutate ka
	{
		i = (int)(mtrand() * n);
		net->Ka[i] *= (4.0 * mtrand() - 2.0);
		return (GAindividual)(net);
	}
	else
	if (r < (MUTATE_KA+MUTATE_PRODUCTION)) //mutate degradation and vmax
	{
		i = (int)(mtrand() * m);
		if (mtrand() > 0.5)
			net->degradation[i] *= (mtrand() * 2.0);
		else
			net->Vmax[i] *= (mtrand() * 2.0);
		return (GAindividual)(net);
	}
	else
	if (r < (MUTATE_KA+MUTATE_PRODUCTION+MUTATE_COMPLEX))
	{
		i = (int)(mtrand() * n);
		j = (int)(mtrand() * net->complexes[i].size);
		net->complexes[i].TFs[j] = (int)(mtrand() * m);
		return (GAindividual)(net);
	}
	else
	if (r < (MUTATE_KA+MUTATE_PRODUCTION+MUTATE_COMPLEX+ADD_GENE))
	{
		++m;
		net->species = m;
		net->numComplexes = n+1;
		complexes = net->complexes;
		targetGene = net->targetGene;
		Ka = net->Ka;
		Vmax = net->Vmax;
		degradation = net->degradation;
		fixed = net->fixed;
		
		net->complexes = (TFComplex*) malloc( (1+n) * sizeof (TFComplex) );
		net->targetGene = (int*) malloc( (1+n) * sizeof(int) );
		net->Ka = (double*) malloc( (1+n) * sizeof(double) );
		net->degradation = (double*) malloc( m * sizeof(double) );
		net->Vmax = (double*) malloc( m * sizeof(double) );
		net->fixed = (int*) malloc( m * sizeof(int) );
		
		for (j=0; j < n; ++j)
		{
			net->targetGene[j] = targetGene[j];
			net->complexes[j] = complexes[j];
			net->Ka[j] = Ka[j];
		}
		
		for (j=0; j < (m-1); ++j)
		{
			net->degradation[j] = degradation[j];
			net->Vmax[j] = Vmax[j];
			net->fixed[j] = fixed[j];
		}
		
		net->complexes[n].size = (int)(1 + TF_RANGE * mtrand());
		net->complexes[n].TFs = (int*) malloc(net->complexes[n].size * sizeof(int));
		net->complexes[n].TFs[0] = (int)(m-1);
		for (i=1; i < net->complexes[n].size; ++i)
		{
			net->complexes[n].TFs[i] = (int)(m * mtrand());
		}
		
		free(targetGene);
		free(complexes);
		free(Ka);
		free(degradation);
		free(Vmax);
		free(fixed);

		net->targetGene[n] = (int)(mtrand() * m);
		net->Ka[n] = (2.0 * mtrand() - 1.0) * KA_RANGE;
		net->degradation[m-1] = mtrand() * DEG_RANGE;
		net->Vmax[m-1] = mtrand() * VMAX_RANGE;
		net->fixed[m-1] = 0;
		
		return (GAindividual)(net);
	}
	else 
	if (m > 2) //remove gene
	{	
		int t = 0, j1 = 0;
		
		i = (int)(mtrand() * m);
		
		--m;
		
		for (j=0; j < n; ++j)
		{
			for (k=0; k < net->complexes[j].size; ++k)
				if (net->complexes[j].TFs[k] == i)
				{
					++t;
					break;
				}
		}
		
		if (n < (t+2)) return (GAindividual)net;
		
		net->species = m;
		complexes = net->complexes;
		targetGene = net->targetGene;
		Ka = net->Ka;
		degradation = net->degradation;
		Vmax = net->Vmax;
		fixed = net->fixed;
		
		net->numComplexes = n-t;
		net->complexes = (TFComplex*) malloc( (n-t) * sizeof (TFComplex) );
		net->targetGene = (int*) malloc( (n-t) * sizeof(int) );
		net->Ka = (double*) malloc( (n-t) * sizeof(double) );
		net->degradation = (double*) malloc( (m) * sizeof(double) );
		net->Vmax = (double*) malloc( (m) * sizeof(double) );
		net->fixed = (int*) malloc( (m) * sizeof(int) );
		
		for (j=0, j1=0; j < (1+m) && j1 < m; ++j)
		{
			if (j != i)
			{
				net->fixed[j1] = fixed[j];
				net->degradation[j1] = degradation[j];
				net->Vmax[j1] = Vmax[j];
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
				net->complexes[j1].size = complexes[j].size;
				net->complexes[j1].TFs = complexes[j].TFs;
				net->Ka[j1] = Ka[j];
				net->targetGene[j1] = targetGene[j];
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
		free(fixed);
		
		return (GAindividual)(net);
	}
	
    return (GAindividual)(net);
}

/*****************************************************
   Functions for simulating and printing
******************************************************/

int getNumSpeciesForGeneRegulationNetwork(GAindividual individual)
{
	GeneRegulationNetwork * net = (GeneRegulationNetwork*)(individual);
	return (net->species);
}

int getNumReactionsForGeneRegulationNetwork(GAindividual individual)
{
	GeneRegulationNetwork * net = (GeneRegulationNetwork*)(individual);
	return (2 * net->species);
}

void setFixedSpeciesForGeneRegulationNetwork(GAindividual individual, int i, int value)
{
	GeneRegulationNetwork * net = (GeneRegulationNetwork*)(individual);
	if (i < net->species)
		net->fixed[i] = value;
}

void ratesForGeneRegulationNetwork(double time,double* u,double* rate,GAindividual p)
{
	int i,j,k;
	double prod, num, denom;
	GeneRegulationNetwork * net;
	
	net = (GeneRegulationNetwork*)(p);
	
	for (i=0; i < net->species; ++i)
	{
		num = denom = 0;
		for (j=0; j < net->numComplexes; ++j)
		{
			if (net->targetGene[j] == i)
			{
				prod = 1.0;
				for (k=0; k < net->complexes[j].size; ++k)  //complex
					prod *= u[ net->complexes[j].TFs[k] ];
				
				if (net->Ka[j] > 0)
				{
					num += net->Ka[j] * prod;
					denom += net->Ka[j] * prod;
				}
				else
				{
					denom += -1.0 * net->Ka[j] * prod;
				}
			}
		}
		if (num == 0) num = 1.0;
		rate[i] = net->Vmax[i] *  num/(1+denom);
		rate[ i + net->species ] = net->degradation[i] * u[i];
	}
}


double * stoichiometryForGeneRegulationNetwork(GAindividual p)
{
	int i,j,m,n;
	double *N;
	GeneRegulationNetwork * net;
	
	net = (GeneRegulationNetwork*)(p);
	
	m = net->species;
	n = 2 * m;
	N = (double*) malloc(m * n * sizeof(double));
	
	for (i=0; i < m; ++i)
		for (j=0; j < n; ++j)
			getValue(N,n,i,j) = 0.0;
			
	for (i=0; i < m; ++i)
		if (net->fixed[i] == 0)
		{
			getValue(N,n,i,i) = 1.0;
			getValue(N,n,i,i+m) = -1.0;
		}
	return N;
}

void printGeneRegulationNetwork(FILE *stream, GAindividual individual)
{
	int i,j,k,p,fix;
	double num, denom;
	GeneRegulationNetwork * net;
	
	if (!individual) return;
	
	net = (GeneRegulationNetwork*)(individual);
	
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
		fprintf(stream,"const s%i",fix);
		for (i=0; i < net->species; ++i)
		{
			if (net->fixed[i])			
				fprintf(stream,", s%i",i+1);			
		}
		fprintf(stream,"\n");
	}
	
	for (i=0; i < net->species; ++i)
	{
		fprintf(stream,"$X -> s%i; vmax%i * (",i+1,i+1);
		num = denom = 0;
		p = 0;
		for (j=0; j < net->numComplexes; ++j)
		{
			if (net->targetGene[j] == i && net->Ka[j] > 0)
			{
				if (p > 0)
					fprintf(stream,"+");
				++p;
				fprintf(stream,"Ka%i",j+1);
				for (k=0; k < net->complexes[j].size; ++k)
					fprintf(stream,"*s%i",net->complexes[j].TFs[k]+1);
			}
		}
		
		if (p < 1) 
			fprintf(stream,"1.0");
		
		fprintf(stream,")/(1.0");
		
		p = 0;
		for (j=0; j < net->numComplexes; ++j)
		{
			if (net->targetGene[j] == i)
			{
				fprintf(stream,"+");
				++p;
				if (net->Ka[j] != 0)
				{
					fprintf(stream,"Ka%i",j+1);
					for (k=0; k < net->complexes[j].size; ++k)
						fprintf(stream,"*s%i",net->complexes[j].TFs[k]+1);
				}
			}
		}
		
		fprintf(stream,");\n");
		fprintf(stream,"s%i -> $X; deg%i*s%i;\n",i+1,i+1,i+1);
	}
	
	fprintf(stream,"\n");
	
	for (i=0; i < net->species; ++i)
	{
		fprintf(stream,"vmax%i = %lf;\n",i+1, net->Vmax[i]);
		fprintf(stream,"deg%i = %lf;\n",i+1,net->degradation[i]);
	}
	for (j=0; j < net->numComplexes; ++j)
	{
		if (net->Ka[j] != 0)
		{
			if (net->Ka[j] < 0)
				fprintf(stream,"Ka%i = %lf;\n",j+1,-net->Ka[j]);
			else
				fprintf(stream,"Ka%i = %lf;\n",j+1,net->Ka[j]);
		}
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
	
	array = (GeneRegulationNetwork**) malloc( num * sizeof(GeneRegulationNetwork*) );
	
	for (k=0; k < num; ++k)
	{
		m = (int)(2 + g * 1.5 * (0.25 + mtrand()));
		n = (int)(2 + r * 1.5 * (0.25 + mtrand()));
		net = newGeneRegulationNetwork(m,n);
		
		for (i=0; i < n; ++i)
		{
			net->complexes[i].size = (int)(1 + TF_RANGE * mtrand());
			net->complexes[i].TFs = (int*) malloc(net->complexes[i].size * sizeof(int));
			for (j=0; j < net->complexes[i].size; ++j)
			{
				net->complexes[i].TFs[j] = (int)(m * mtrand());
			}
			net->targetGene[i] = (int)(mtrand() * m);
			net->Ka[i] = (2.0 * mtrand() - 1.0) * KA_RANGE;
		}
		for (i=0; i < m; ++i)
		{
			net->degradation[i] = mtrand() * DEG_RANGE;
			net->Vmax[i] = mtrand() * VMAX_RANGE;
		}
		array[k] = net;
	}
	
	return (GApopulation)(array);
}

GeneRegulationNetwork * newGeneRegulationNetwork(int m,int n)
{
	int i;
	GeneRegulationNetwork * net;
	
	net = (GeneRegulationNetwork*) malloc(sizeof(GeneRegulationNetwork));
	net->species = m;    //number of genes
	net->numComplexes = n;    //number of complexes
	
	net->complexes = (TFComplex*) malloc( n * sizeof (TFComplex) );
	net->targetGene = (int*) malloc( n * sizeof(int) );
	net->Ka = (double*) malloc( n * sizeof(double) );
	
	net->Vmax = 0;
	net->degradation = 0;
	net->fixed = 0;

	if (m > 0)
	{
		net->Vmax = (double*) malloc( m * sizeof(double) );
		net->degradation = (double*) malloc( m * sizeof(double) );
		net->fixed = (int*) malloc( m * sizeof(int) );
		
		for (i=0; i < n; ++i)
		{
			net->complexes[i].size = 0;
			net->complexes[i].TFs = 0;		
			net->targetGene[i] = 0;
			net->Ka[i] = 0.0;
		}
		for (i=0; i < m; ++i)
		{
			net->degradation[i] = 0.0;
			net->Vmax[i] = 0.0;
			net->fixed[i] = 0;
		}
	}
	return net;
}

