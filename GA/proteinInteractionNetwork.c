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
	double total;

	total = a+b+c+d;
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
	int i;
	ProteinInteractionNetwork * net;
	
	if (!individual) return;
	net = (ProteinInteractionNetwork*)(individual);
	
	if (net->species < 1) return;

	if (net->regulators)
	{
		for (i=0; i < net->species; ++i)
		{
			if (net->regulators[i].size > 0)
			{
				if (net->regulators[i].proteins) free (net->regulators[i].proteins);
				if (net->regulators[i].Km) free (net->regulators[i].Km);
				if (net->regulators[i].Vmax) free (net->regulators[i].Vmax);
			}
		}
		free (net->regulators);
	
	if (net->totals)
		free (net->totals);
	
	if (net->fixed)
		free (net->fixed);
		
	}
}

void* cloneProteinInteractionNetwork(void * individual)
{
	int i,j,m,n;
	ProteinInteractionNetwork * net, * net2;
	
	if (!individual) return 0;
	
	net = (ProteinInteractionNetwork*)(individual);   //original
	net2 = malloc(sizeof(ProteinInteractionNetwork)); //soon to be clone
	
	n = net->species;    //number of species
	net2->species = n;
	
	net2->regulators = malloc(n * sizeof(Regulators));   //allocate space
	net2->totals = malloc(n * sizeof(double));
	net2->fixed = malloc(n * sizeof(int));
	
	for (i=0; i < n; ++i)   //copy regulators
	{
		net2->fixed[i] = net->fixed[i];
		net2->totals[i] = net->totals[i];
		m = net->regulators[i].size;
		net2->regulators[i].size = m;
		net2->regulators[i].proteins = malloc(m * sizeof(int));
		net2->regulators[i].Km = malloc(m * sizeof(double));
		net2->regulators[i].Vmax = malloc(m * sizeof(double));
		for (j=0; j < m; ++j)  //copy values for each regulator
		{
			net2->regulators[i].proteins[j] = net->regulators[i].proteins[j];
			net2->regulators[i].Km[j] = net->regulators[i].Km[j];
			net2->regulators[i].Vmax[j] = net->regulators[i].Vmax[j];
		}
	}
	
	return (void*)(net2);  //done
}

void* crossoverProteinInteractionNetwork(void * individualA, void * individualB)
{
	int i, j, k, i1, i2, n;
	ProteinInteractionNetwork * net1, * net2, * net3;

	if (mtrand() > CROSSOVER_PROB) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(individualA)); 
	
	if (!individualA) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(individualB));
	if (!individualB) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(individualA));
	
	net1 = (ProteinInteractionNetwork*)(individualA);  //parents
	net2 = (ProteinInteractionNetwork*)(individualB);
	
	if (net1->species < 3) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(net2));  //if parents are too small
	if (net2->species < 3) return mutateProteinInteractionNetwork(cloneProteinInteractionNetwork(net1));
	
	i1 = (int)(mtrand() * (net1->species - 1) + 1.0);	//crossover point in net1
	i2 = (int)(mtrand() * (net2->species - 2) + 1.0);	//crossover point in net2
	
	n = i1 + net2->species - i2;
	
	net3 = malloc(sizeof(ProteinInteractionNetwork));  //child network
	net3->species = n;
	net3->regulators = malloc(n * sizeof(Regulators));
	net3->totals = malloc(n * sizeof(double));
	net3->fixed = malloc(n * sizeof(int));
	
	for (i=0; i < i1; ++i) //copy all regulators from net1
	{
		n = net1->regulators[i].size;
		net3->regulators[i].size = n;
		net3->totals[i] = net1->totals[i];
		net3->fixed[i] = net1->fixed[i];
		net3->regulators[i].proteins = malloc(n * sizeof(int));
		net3->regulators[i].Km = malloc(n * sizeof(double));
		net3->regulators[i].Vmax = malloc(n * sizeof(double));
		for (j=0; j < n; ++j)
		{
			net3->regulators[i].proteins[j] = net1->regulators[i].proteins[j];
			if (net3->regulators[i].proteins[j] >= net3->species)
				net3->regulators[i].proteins[j] = (int)(mtrand() * net3->species);
			net3->regulators[i].Km[j] = net1->regulators[i].Km[j];
			net3->regulators[i].Vmax[j] = net1->regulators[i].Vmax[j];
		}
	}
	
	for (i=i2; i < net2->species; ++i)
	{
		k = i+i1-i2;
		n = net2->regulators[i].size;
		net3->regulators[k].size = n;
		net3->totals[k] = net2->totals[i];
		net3->fixed[k] = net2->fixed[i];
		net3->regulators[k].proteins = malloc(n * sizeof(int));
		net3->regulators[k].Km = malloc(n * sizeof(double));
		net3->regulators[k].Vmax = malloc(n * sizeof(double));
		for (j=0; j < n; ++j)
		{
			net3->regulators[k].proteins[j] = net2->regulators[i].proteins[j];
			if (net3->regulators[k].proteins[j] >= net3->species)
				net3->regulators[k].proteins[j] = (int)(mtrand() * net3->species);
			net3->regulators[k].Km[j] = net2->regulators[i].Km[j];
			net3->regulators[k].Vmax[j] = net2->regulators[i].Vmax[j];
		}
	}
	return (void*)(net3);
}

void* mutateProteinInteractionNetwork(void * individual)
{
	int i,j,j2,k,m,n,n2, * fixed;
	double r, * totals;
	ProteinInteractionNetwork * net;
	Regulators * regulators;

	net = (ProteinInteractionNetwork*)individual;

	n = net->species;

	i = (int)(mtrand() * n);  //pick random protein
	r = mtrand();

	if (r < MUTATE_REWIRE)   //mutate one of the regulators
	{
		j = (int)(mtrand() * net->regulators[i].size);
		if (mtrand() < MUTATE_REWIRE) //rewire
			net->regulators[i].proteins[j] = (int)(mtrand() * net->species);
	}
	else
	if (r < (MUTATE_REWIRE+MUTATE_CHANGE_PARAM)) //change parameter
	{
		j = (int)(mtrand() * net->regulators[i].size);
		
		if (mtrand() < 0.5)
			net->regulators[i].Km[j]  *= (mtrand() * 2.0);
		else
			net->regulators[i].Vmax[j]  *= (mtrand() * 2.0);
		
		return (void*)(net);
	}
	if (r < (MUTATE_REWIRE+MUTATE_CHANGE_PARAM+MUTATE_TOTAL_CONC))  //change total concentrations (conservation rule)
	{
		net->totals[i] *= (2.0 * mtrand());
		return (void*)(net);
	}
	else              //add or remove a new protein to the network
	{
		if (mtrand() < 0.5 && n > 2)     //remove a protein
		{
			n2 = n-1;
			net->species = n2;
			regulators = net->regulators;
			totals = net->totals;
			fixed = net->fixed;

			net->regulators = malloc( n2 * sizeof(Regulators) );
			net->totals = malloc( n2 * sizeof(double) );
			net->fixed = malloc( n2 * sizeof(int) );
			
			for (j=0,j2=0; j < n; ++j) //copy all proteins
			{
				m = regulators[j].size;

				if (j != i) //except protein i
				{
					net->totals[j2] = totals[j];
					net->fixed[j2] = fixed[j];
					net->regulators[j2].size = m;
					net->regulators[j2].proteins = malloc(m * sizeof(int));
					net->regulators[j2].Km = malloc(m * sizeof(double));
					net->regulators[j2].Vmax = malloc(m * sizeof(double));
					for (k=0; k < m; ++k)  //copy values for each regulator
					{
						net->regulators[j2].proteins[k] = regulators[j].proteins[k];
						net->regulators[j2].Km[k] = regulators[j].Km[k];
						net->regulators[j2].Vmax[k] = regulators[j].Vmax[k];
					}
					++j2;
				}

				if (m != 0)
				{
					free(regulators[j].Km);
					free(regulators[j].Vmax);
					free(regulators[j].proteins);
				}
			}

			if (n != 0)
			{
				free(fixed);
				free(totals);
				free(regulators);
			}

			return (void*)(net);
		}
		else	//add a protein
		{
			n2 = n+1;
			regulators = net->regulators;
			totals = net->totals;
			fixed = net->fixed;

			net->species = n2;
			net->regulators = malloc( n2 * sizeof(Regulators) );
			net->totals = malloc( n2 * sizeof(double) );
			net->fixed = malloc( n2 * sizeof(int) );
			
			for (j=0; j < n; ++j) //copy all proteins
			{
				net->totals[j] = totals[j];
				net->fixed[j] = fixed[j];
				m = regulators[j].size;
				net->regulators[j].size = m;
				net->regulators[j].proteins = malloc(m * sizeof(int));
				net->regulators[j].Vmax = malloc(m * sizeof(double));
				net->regulators[j].Km = malloc(m * sizeof(double));
				for (k=0; k < m; ++k)  //copy values for each regulator
				{
					net->regulators[j].proteins[k] = regulators[j].proteins[k];
					net->regulators[j].Km[k] = regulators[j].Km[k];
					net->regulators[j].Vmax[k] = regulators[j].Vmax[k];
				}

				if (m != 0)
				{
					free(regulators[j].Km);
					free(regulators[j].Vmax);
					free(regulators[j].proteins);
				}
			}

			if (n != 0)
			{
				free(fixed);
				free(totals);
				free(regulators);
			}

			//the new protein
			net->totals[n] = (2.0 * mtrand()) * AVG_TOTAL;
			m = (int)(mtrand() * n * AVG_NUM_REGULATIONS);
			net->regulators[n].size = m;
			net->regulators[n].proteins = malloc(m * sizeof(int));
			net->regulators[n].Km = malloc(m * sizeof(double));
			net->regulators[n].Vmax = malloc(m * sizeof(double));
			for (j=0; j < m; ++j)  //random values for the new protein
			{
				net->regulators[n].proteins[j] = (int)(mtrand() * net->species);
				net->regulators[n].Km[j] = mtrand() * KM_RANGE;
				net->regulators[n].Vmax[j] = mtrand() * VMAX_RANGE;
				
				if (mtrand() < 0.5)
					net->regulators[n].Vmax[j] *= -1.0;  //negative regulator
			}

			return (void*)(net);
		}
	}

    return (net);
}

/*****************************************************
   Functions for simulating and printing
******************************************************/


int getNumSpeciesForProteinInteractionNetwork(void * individual)
{
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);
	return (net->species);
}

int getNumReactionsForProteinInteractionNetwork(void * individual)
{
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);
	return (2 * net->species);
}

void setFixedSpeciesForProteinInteractionNetwork(void * individual, int i, int value)
{
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);
	if (i < net->species)
		net->fixed[i] = value;
}

void ratesForProteinInteractionNetwork(double time,double* u,double* rate,void * individual)
{
	int i,j,n,p,forward, backward;
	double km,vmax,tot;
	ProteinInteractionNetwork * net;
	
	net = (ProteinInteractionNetwork*)(individual);
	n = net->species;
	
	for (i=0; i < n; ++i)
	{
		tot = net->totals[i]; //total concentration of u[i]
		forward = backward = 0;
		rate[2*i] = 0;
		rate[2*i+1] = 0;
		for (j=0; j < net->regulators[i].size; ++j)
		{
			vmax = net->regulators[i].Vmax[j]; //vmax for this regulation
			km = net->regulators[i].Km[j];  //km for this regulation
			p = net->regulators[i].proteins[j]; //index of regulating protein
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
			rate[2*i] += VMAX_RANGE/2.0 * u[i] / (tot/2.0 + u[i]);
		
		if (!backward)
			rate[2*i+1] += VMAX_RANGE/2.0 * (tot - u[i]) / (tot/2.0 + (tot - u[i]));
	}
}


double * stoichiometryForProteinInteractionNetwork(void * p)
{
	int i,n;
	double * N;
	ProteinInteractionNetwork * net;
	
	net = (ProteinInteractionNetwork*)(p);
	n = net->species;
	N = malloc(n * 2 * n * sizeof(double));
	for (i=0; i < (2*n*n); ++i)
	{
		N[i] = 0.0;
	}
	for (i=0; i < n; ++i)
		if (net->fixed[i] == 0)
		{
			getValue(N,(2*n),i,i*2) = -1.0;
			getValue(N,(2*n),i,i*2+1) = 1.0;
		}
	return N;
}

void printProteinInteractionNetwork(void * individual)
{
	int i,j,n,p,fix;
	double km,vmax,tot,f;
	ProteinInteractionNetwork * net;
	
	net = (ProteinInteractionNetwork*)(individual);
	n = net->species;
	
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
	
	for (i=0; i < n; ++i)
	{
		printf("$X -> s%i; ", i+1);
		f = 0;
		for (j=0; j < net->regulators[i].size; ++j)
		{
			p = net->regulators[i].proteins[j]; //index of regulating protein
			vmax = net->regulators[i].Vmax[j]; //vmax for this regulation

			if (vmax > 0)
				if (!f)
				{
					printf("vmax%i%i * s%i * (tot%i - s%i) / (km%i%i + (tot%i - s%i)) ",i+1,j+1,i+1,i+1,p+1,i+1,j+1,i+1,i+1);
					f = 1;
				}
				else
					printf(" + %lf * s%i * (tot%i - s%i) / (km%i%i + (tot%i - s%i)) ",i+1,j+1,i+1,i+1,p+1,i+1,j+1,i+1,i+1);
		}
		
		for (j=0; j < net->regulators[i].size; ++j)
		{
			vmax = net->regulators[i].Vmax[j]; //vmax for this regulation
			p = net->regulators[i].proteins[j]; //index of regulating protein
			if (vmax < 0) 
			{
				printf(" - vmax%i%i * s%i * s%i / (km%i%i + s%i)",i+1,j+1,p+1,i+1,i+1,j+1,i+1);
				if (!f) f = 1;
			}
		}
		
		if (!f) 
			printf("0;\n");
		else
			printf(";\n");
	}
	
	printf("\n");
	
	for (i=0; i < n; ++i)
	{
		tot = net->totals[i]; //total concentration of u[i]
		printf("tot%i = %lf;\n",i,tot);
		
		for (j=0; j < net->regulators[i].size; ++j)
		{
			vmax = net->regulators[i].Vmax[j]; //vmax for this regulation
			km = net->regulators[i].Km[j];  //km for this regulation
			if (vmax != 0)
				printf("vmax%i%i = %lf;\n ",vmax,i+1,j+1,vmax);
			if (km != 0)
				printf("km%i%i = %lf;\n ",km,i+1,j+1,km);	
		}
	}
}

/***********************
  GA related functions
***********************/

GApopulation randomProteinInteractionNetworks(int num)
{
	int s = AVG_NUM_SPECIES;
	int i,j,k,n,m;
	ProteinInteractionNetwork * net;
	ProteinInteractionNetwork ** array;
	
	initMTrand(); /*initialize seeds for MT random number generator*/
	
	array = malloc(num * sizeof(ProteinInteractionNetwork*));
	for (i=0; i < num; ++i)
	{
		n = (int)(1 + s * 2.0 * mtrand());
		net = malloc(sizeof(ProteinInteractionNetwork)); //new network
	
		net->species = n;    //number of proteins
		net->fixed = malloc(n * sizeof(int));
		
		net->regulators = malloc(n * sizeof(Regulators));   //allocate space
		net->totals = malloc(n * sizeof(double));
		
		for (j=0; j < n; ++j)   //random regulators for each protein
		{
			net->fixed[j] = 0; //no fixed species by default
			net->totals[j] = 2.0*mtrand()*AVG_TOTAL;
			m = (int)(2 + mtrand() * n * AVG_NUM_REGULATIONS);
			net->regulators[j].size = m;
			net->regulators[j].proteins = malloc(m * sizeof(int));
			net->regulators[j].Km = malloc(m * sizeof(double));
			net->regulators[j].Vmax = malloc(m * sizeof(double));
			for (k=0; k < m; ++k)  //random values for the new protein
			{
				net->regulators[j].proteins[k] = (int)(mtrand() * net->species);
				net->regulators[j].Km[k] = mtrand() * KM_RANGE;
				net->regulators[j].Vmax[k] = mtrand() * VMAX_RANGE;
				
				if (mtrand() < 0.5 || k == (m-1))
					net->regulators[j].Vmax[k] *= -1.0;  //negative regulator
			}
		}
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
		deleteProteinInteractionNetwork(net);
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
