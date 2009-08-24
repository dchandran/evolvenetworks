/*******************************************************

	Copyright (C) 2009 Deepak Chandran
	see header file

********************************************************/

#include "proteinInteractionNetwork.h"

/*************
global parameters
**************/

static double AVG_TOTAL = 10.0;
static double KM_RANGE = 10.0;
static double VMAX_RANGE = 10.0;
static double AVG_NUM_REGULATIONS = 0.2;
static double MUTATE_REWIRE = 0.2;
static double MUTATE_CHANGE_PARAM = 0.5;
static double MUTATE_TOTAL_CONC = 0.15;
static double CROSSOVER_PROB = 1.0;
static int AVG_NUM_SPECIES = 10;

void setParametersForProteinInteractionNetwork(double ka, double vmax, double total)
{
	KM_RANGE = ka;
	VMAX_RANGE = vmax;
	AVG_TOTAL = total;
}

void setSizeForProteinInteractionNetwork(int s, int r)
{
	AVG_NUM_SPECIES = s;
	AVG_NUM_REGULATIONS = (double)r/(double)s;
}

void setMutationAndCrossoverRatesForProteinInteractionNetwork(double a, double b, double c, double d, double e)
{
	double total = a+b+c+d;
	MUTATE_REWIRE = a/total;
	MUTATE_CHANGE_PARAM = b/total;
	MUTATE_TOTAL_CONC = c/total;
	
	if (e > 1.0) e /= 100.0;
	CROSSOVER_PROB = e;
}

/********************************************************
    Clone, delete, mutate, crossover  (required by GA)
*********************************************************/

void deleteProteinInteractionNetwork(void * individual)
{
	if (!individual) return;
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);
	
	int i;
	if ((*net).regulators)
	{
		for (i=0; i < (*net).species; ++i)
		{
			if ((*net).regulators[i].proteins) free ((*net).regulators[i].proteins);
			if ((*net).regulators[i].Km) free ((*net).regulators[i].Km);
			if ((*net).regulators[i].Vmax) free ((*net).regulators[i].Vmax);
		}
		free ((*net).regulators);
	
	if ((*net).totals)
		free ((*net).totals);
	}
}

void* cloneProteinInteractionNetwork(void * individual)
{
	if (!individual) return 0;
	
	int i,j,m,n;
	
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);   //original
	ProteinInteractionNetwork * net2 = malloc(sizeof(ProteinInteractionNetwork)); //soon to be clone
	
	n = (*net).species;    //number of species
	(*net2).species = n;
	
	(*net2).regulators = malloc(n * sizeof(Regulators));   //allocate space
	(*net2).totals = malloc(n * sizeof(double));
	
	for (i=0; i < n; ++i)   //copy regulators
	{
		(*net2).totals[i] = (*net).totals[i];
		m = (*net).regulators[i].size;
		(*net2).regulators[i].size = m;
		(*net2).regulators[i].proteins = malloc(m * sizeof(int));
		(*net2).regulators[i].Km = malloc(m * sizeof(double));
		(*net2).regulators[i].Vmax = malloc(m * sizeof(double));
		for (j=0; j < m; ++j)  //copy values for each regulator
		{
			(*net2).regulators[i].proteins[j] = (*net).regulators[i].proteins[j];
			(*net2).regulators[i].Km[j] = (*net).regulators[i].Km[j];
			(*net2).regulators[i].Vmax[j] = (*net).regulators[i].Vmax[j];
		}
	}
	
	return (void*)(net2);  //done
}

void* crossoverProteinInteractionNetwork(void * individualA, void * individualB)
{
	int i, j, k, i1, i2, n = 0, m = 0;
	
	if (mtrand() > CROSSOVER_PROB) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(individualA)); 
	
	if (!individualA) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(individualB));
	if (!individualB) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(individualA));
	
	ProteinInteractionNetwork * net1 = (ProteinInteractionNetwork*)(individualA);  //parents
	ProteinInteractionNetwork * net2 = (ProteinInteractionNetwork*)(individualB);
	
	if ((*net1).species < 3) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(net2));  //if parents are too small
	if ((*net2).species < 3) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(net1));
	
	i1 = (int)(mtrand() * ((*net1).species - 1) + 1.0);	//crossover point in net1
	i2 = (int)(mtrand() * ((*net2).species - 2) + 1.0);	//crossover point in net2
	
	n = i1 + (*net2).species - i2;
	
	ProteinInteractionNetwork * net3 = malloc(sizeof(ProteinInteractionNetwork));  //child network
	(*net3).species = n;
	(*net3).regulators = malloc(n * sizeof(Regulators));
	(*net3).totals = malloc(n * sizeof(double));
	
	for (i=0; i < i1; ++i) //copy all regulators from net1
	{
		n = (*net1).regulators[i].size;
		(*net3).regulators[i].size = n;
		(*net3).totals[i] = (*net1).totals[i];
		(*net3).regulators[i].proteins = malloc(n * sizeof(int));
		(*net3).regulators[i].Km = malloc(n * sizeof(double));
		(*net3).regulators[i].Vmax = malloc(n * sizeof(double));
		for (j=0; j < n; ++j)
		{
			(*net3).regulators[i].proteins[j] = (*net1).regulators[i].proteins[j];
			if ((*net3).regulators[i].proteins[j] >= (*net3).species)
				(*net3).regulators[i].proteins[j] = (int)(mtrand() * (*net3).species);
			(*net3).regulators[i].Km[j] = (*net1).regulators[i].Km[j];
			(*net3).regulators[i].Vmax[j] = (*net1).regulators[i].Vmax[j];
		}
	}
	
	for (i=i2; i < (*net2).species; ++i)
	{
		k = i+i1-i2;
		n = (*net2).regulators[i].size;
		(*net3).regulators[k].size = n;
		(*net3).totals[k] = (*net2).totals[i];
		(*net3).regulators[k].proteins = malloc(n * sizeof(int));
		(*net3).regulators[k].Km = malloc(n * sizeof(double));
		(*net3).regulators[k].Vmax = malloc(n * sizeof(double));
		for (j=0; j < n; ++j)
		{
			(*net3).regulators[k].proteins[j] = (*net2).regulators[i].proteins[j];
			if ((*net3).regulators[k].proteins[j] >= (*net3).species)
				(*net3).regulators[k].proteins[j] = (int)(mtrand() * (*net3).species);
			(*net3).regulators[k].Km[j] = (*net2).regulators[i].Km[j];
			(*net3).regulators[k].Vmax[j] = (*net2).regulators[i].Vmax[j];
		}
	}
	
	return (void*)(net3);
}

void* mutateProteinInteractionNetwork(void * individual)
{
	int i,j,j2,k,m,n;
	
	double r = mtrand();
	
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)individual;

	n = (*net).species;

	i = (int)(mtrand() * n);  //pick random protein

	

	if (r < MUTATE_REWIRE)   //mutate one of the regulators
	{
		j = (int)(mtrand() * (*net).regulators[i].size);
		if (mtrand() < MUTATE_REWIRE) //rewire
			(*net).regulators[i].proteins[j] = (int)(mtrand() * (*net).species);
	}
	else
	if (r < (MUTATE_REWIRE+MUTATE_CHANGE_PARAM)) //change parameter
	{
		j = (int)(mtrand() * (*net).regulators[i].size);
		
		if (mtrand() < 0.5)
			(*net).regulators[i].Km[j]  *= (mtrand() * 2.0);
		else
			(*net).regulators[i].Vmax[j]  *= (mtrand() * 2.0);
		
		return (void*)(net);
	}
	else
	if (r < (MUTATE_REWIRE+MUTATE_CHANGE_PARAM+MUTATE_TOTAL_CONC))  //change total concentrations (conservation rule)
	{
		(*net).totals[i] *= (2.0 * mtrand());
		return (void*)(net);
	}
	else              //add or remove a new protein to the network
	{
		j = mtrand();
		if (j < 0.5 && n > 2)     //remove a protein
		{
			ProteinInteractionNetwork * net2 = malloc(sizeof(ProteinInteractionNetwork));
			(*net2).species = n-1;
			(*net2).regulators = malloc( (n-1) * sizeof(Regulators) );
			(*net2).totals = malloc( (n-1) * sizeof(double) );
			
			for (j=0,j2=0; j < n; ++j) //copy all proteins
			{
				if (j != i) //except protein i
				{
					(*net2).totals[j2] = (*net).totals[j];
					m = (*net).regulators[j].size;
					(*net2).regulators[j2].size = m;
					(*net2).regulators[j2].proteins = malloc(m * sizeof(int));
					(*net2).regulators[j2].Km = malloc(m * sizeof(double));
					(*net2).regulators[j2].Vmax = malloc(m * sizeof(double));
					for (k=0; k < m; ++k)  //copy values for each regulator
					{
						(*net2).regulators[j2].proteins[k] = (*net).regulators[j].proteins[k];
						(*net2).regulators[j2].Km[k] = (*net).regulators[j].Km[k];
						(*net2).regulators[j2].Vmax[k] = (*net).regulators[j].Vmax[k];
					}
					++j2;
				}
			}
			deleteNetwork(net);
			
			return (void*)(net2);
		}
		else	//add a protein
		{
			ProteinInteractionNetwork * net2 = malloc(sizeof(ProteinInteractionNetwork));
			(*net2).species = n+1;
			(*net2).regulators = malloc( (n+1) * sizeof(Regulators) );
			(*net2).totals = malloc( (n+1) * sizeof(double) );
			
			for (j=0; j < n; ++j) //copy all proteins
			{
				(*net2).totals[j] = (*net).totals[j];
				m = (*net).regulators[j].size;
				(*net2).regulators[j].size = m;
				(*net2).regulators[j].proteins = malloc(m * sizeof(int));
				(*net2).regulators[j].Km = malloc(m * sizeof(double));
				(*net2).regulators[j].Vmax = malloc(m * sizeof(double));
				for (k=0; k < m; ++k)  //copy values for each regulator
				{
					(*net2).regulators[j].proteins[k] = (*net).regulators[j].proteins[k];
					(*net2).regulators[j].Km[k] = (*net).regulators[j].Km[k];
					(*net2).regulators[j].Vmax[k] = (*net).regulators[j].Vmax[k];
				}
			}
			
			//the new protein
			(*net2).totals[n] = (2.0 * mtrand()) * AVG_TOTAL;
			m = (*net2).regulators[n].size = (int)(mtrand() * n * AVG_NUM_REGULATIONS);
			(*net2).regulators[n].proteins = malloc(m * sizeof(int));
			(*net2).regulators[n].Km = malloc(m * sizeof(double));
			(*net2).regulators[n].Vmax = malloc(m * sizeof(double));
			for (j=0; j < m; ++j)  //random values for the new protein
			{
				(*net2).regulators[n].proteins[j] = (int)(mtrand() * (*net2).species);
				(*net2).regulators[n].Km[j] = mtrand() * KM_RANGE;
				(*net2).regulators[n].Vmax[j] = mtrand() * VMAX_RANGE;
				
				if (mtrand() < 0.5)
					(*net2).regulators[n].Vmax[j] *= -1.0;  //negative regulator
			}
			
			deleteNetwork(net);
			return (void*)(net2);
		}
	}

    return (net);
}

/*****************************************************
   Functions for simulating and printing
******************************************************/

void ratesForProteinInteractionNetwork(double time,double* u,double* rate,void * individual)
{
	int i,j,n,p;
	double km,vmax,tot,f;
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);
	n = (*net).species;
	int forward = 0, backward = 0;
	for (i=0; i < n; ++i)
	{
		tot = (*net).totals[i]; //total concentration of u[i]
		forward = backward = 0;
		rate[2*i] = 0;
		rate[2*i+1] = 0;
		for (j=0; j < (*net).regulators[i].size; ++j)
		{
			vmax = (*net).regulators[i].Vmax[j]; //vmax for this regulation
			km = (*net).regulators[i].Km[j];  //km for this regulation
			p = (*net).regulators[i].proteins[j]; //index of regulating protein
			if (vmax < 0) 
			{
				forward = 1;
				rate[2*i] += -vmax * u[p] * u[i] / (km + u[i]);
			}
			else
			{
				backward = 1;
				rate[2*i+1] += vmax * u[p] * (tot - u[i]) / (km + (tot - u[i]) );
			}
		}
		if (!forward)
			rate[2*i] += VMAX_RANGE/2.0 * u[i] / (km + u[i]);
		
		if (!backward)
			rate[2*i+1] += VMAX_RANGE/2.0 * (tot - u[i]) / (km + (tot - u[i]));
	}
}


double * stoichiometryForProteinInteractionNetwork(void * p)
{
	int i,j,n;
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(p);
	n = (*net).species;
	double * N = malloc(n * 2 * n * sizeof(double));
	for (i=0; i < (2*n*n); ++i)
	{
		N[i] = 0.0;
	}
	for (i=0; i < n; ++i)
	{
		getValue(N,(2*n),i,i*2) = -1.0;
		getValue(N,(2*n),i,i*2+1) = 1.0;
	}
	return N;
}

void printProteinInteractionNetwork(void * individual)
{
	int i,j,n,p;
	double km,vmax,tot,f;
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);
	n = (*net).species;
	
	for (i=0; i < n; ++i)
	{
		printf("$X -> s%i; ", i+1);
		f = 0;
		for (j=0; j < (*net).regulators[i].size; ++j)
		{
			tot = (*net).totals[i]; //total concentration of u[i]
			vmax = (*net).regulators[i].Vmax[j]; //vmax for this regulation
			km = (*net).regulators[i].Km[j];  //km for this regulation
			p = (*net).regulators[i].proteins[j]; //index of regulating protein
			if (vmax > 0)
				if (!f)
				{
					printf("%lf * s%i * (%lf - s%i) / (%lf + (%lf - s%i)) ",vmax,i+1,tot,p+1,km,tot,i+1);
					f = 1;
				}
				else
					printf(" + %lf * s%i * (%lf - s%i) / (%lf + (%lf - s%i)) ",vmax,i+1,tot,p+1,km,tot,i+1);
		}
		
		for (j=0; j < (*net).regulators[i].size; ++j)
		{
			tot = (*net).totals[i]; //total concentration of u[i]
			vmax = (*net).regulators[i].Vmax[j]; //vmax for this regulation
			km = (*net).regulators[i].Km[j];  //km for this regulation
			p = (*net).regulators[i].proteins[j]; //index of regulating protein
			if (vmax < 0) 
			{
				printf(" - %lf * s%i * s%i / (%lf + s%i)",-vmax,p+1,i+1,km,i+1);
				if (!f) f = 1;
			}
		}
		
		if (!f) 
			printf("0;\n");
		else
			printf(";\n");
	}
}

/***********************
  GA related functions
***********************/

GApopulation randomProteinInteractionNetworks(int num)
{
	int s = AVG_NUM_SPECIES;
	int i,j,k,n,m;
	initMTrand(); /*initialize seeds for MT random number generator*/
	ProteinInteractionNetwork ** array = malloc(num * sizeof(ProteinInteractionNetwork*));
	for (i=0; i < num; ++i)
	{
		n = (int)(1 + s * 2.0 * mtrand());
		ProteinInteractionNetwork * net = malloc(sizeof(ProteinInteractionNetwork)); //new network
	
		(*net).species = n;    //number of proteins
		
		(*net).regulators = malloc(n * sizeof(Regulators));   //allocate space
		(*net).totals = malloc(n * sizeof(double));
		
		for (j=0; j < n; ++j)   //random regulators for each protein
		{
			(*net).totals[j] = 2.0*mtrand()*AVG_TOTAL;
			m = (int)(2 + mtrand() * n * AVG_NUM_REGULATIONS);
			(*net).regulators[j].size = m;
			(*net).regulators[j].proteins = malloc(m * sizeof(int));
			(*net).regulators[j].Km = malloc(m * sizeof(double));
			(*net).regulators[j].Vmax = malloc(m * sizeof(double));
			for (k=0; k < m; ++k)  //random values for the new protein
			{
				(*net).regulators[j].proteins[k] = (int)(mtrand() * (*net).species);
				(*net).regulators[j].Km[k] = mtrand() * KM_RANGE;
				(*net).regulators[j].Vmax[k] = mtrand() * VMAX_RANGE;
				
				if (mtrand() < 0.5 || k == (m-1))
					(*net).regulators[j].Vmax[k] *= -1.0;  //negative regulator
			}
		}
		
		//printNetwork(net);
		
		array[i] = net;
	}
	
	return (GApopulation)(array);
}

/*
int main()  //just for testing
{
	int i,j;
	
	int n = 100;
	printf("random start...\n");
	GApopulation pop = randomProteinInteractionNetworks(n,8);
	
	printf("print initial nets...\n");
	for (i=0; i < n; ++i)
	{
		printProteinInteractionNetwork(pop[i]);
		printf("\n");
	}
	
	printf("crossover and mutate many times...\n");
	for (i=0; i < 100; ++i)
	{
		void * net = crossoverProteinInteractionNetwork(pop[0],pop[1]);
		printNetwork(net);
		net = mutateProteinInteractionNetwork(net);
		printf("\n");
		deleteNetwork(net);
	}
	
	for (i=0; i < n; ++i)
	{
		deleteProteinInteractionNetwork(pop[i]);
	}
	
	free(pop);
	
	initMTrand();
	ProteinInteractionNetwork p;
	int n = 3,m;
	p.species = n;
	p.regulators = malloc(n * sizeof(Regulators));   //allocate space
	p.totals = malloc(n * sizeof(double));
		
	for (j=0; j < n; ++j)
	{
		p.totals[j] = 10.0;
		m = 1;
		p.regulators[j].size = m;
		p.regulators[j].proteins = malloc(m * sizeof(int));
		p.regulators[j].Km = malloc(m * sizeof(double));
		p.regulators[j].Vmax = malloc(m * sizeof(double));
		
		p.regulators[j].proteins[0] = j + 1;
		p.regulators[j].Km[0] = 0.1;
		p.regulators[j].Vmax[0] = 1.0;
		if (j + 1 >= n)
		{
			p.regulators[j].proteins[0] = 0;
			p.regulators[j].Vmax[0] = -1.0;
		}		
	}
	
	double init[] = {0,0,0};
	int sz = 1000;
	double * N = stoichiometryForProteinInteractionNetwork(&p);
	//double * y = ODEsim2(n,2*n,N,&ratesForProteinInteractionNetwork,init,0,100,0.1,&p);	
	double * y = SSA(n,2*n,N,&ratesForProteinInteractionNetwork,init,0,100,100000,&sz,&p);
	
	free(N);
	writeToFile("data.tab",y,sz,4);
	free(y);
	return 0;
}


*/
