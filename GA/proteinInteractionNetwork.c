/*******************************************************

	Copyright (C) 2009 Deepak Chandran
	see header file

********************************************************/

#include "proteinInteractionNetwork.h"

/*************
global parameters
**************/

static double MIN_TOTAL = 0.001;
static double MAX_TOTAL = 10.0;
static double KM_MIN = 0.0001;
static double KM_MAX = 10.0;
static double VMAX_MIN = 0.0001;
static double VMAX_MAX = 10.0;
static double MIN_NUM_REGULATIONS = 0.1;
static double MAX_NUM_REGULATIONS = 0.5;
static double MUTATE_REWIRE = 0.2;
static double MUTATE_CHANGE_PARAM = 0.5;
static double MUTATE_TOTAL_CONC = 0.15;
static double CROSSOVER_PROB = 1.0;
static int MIN_NUM_SPECIES = 4;
static int MAX_NUM_SPECIES = 16;

void setRateConstantsForProteinInteractionNetwork(double min_ka, double max_ka, double min_vmax, double max_vmax, double min_total, double max_total)
{
	double d;
	if (min_ka > 0.0) KM_MIN = min_ka;
	if (max_ka > 0.0) KM_MAX = min_ka;
	if (KM_MIN > KM_MAX)
	{
		d = KM_MIN;
		KM_MIN = KM_MAX;
		KM_MAX = d;
	}
	if (min_ka > 0.0) KM_MIN = min_ka;
	if (max_ka > 0.0) KM_MAX = min_ka;
	if (KM_MIN > KM_MAX)
	{
		d = KM_MIN;
		KM_MIN = KM_MAX;
		KM_MAX = d;
	}
void setRateConstantsForProteinInteractionNetwork(double ka, double vmax, double total)
{
	KM_RANGE = ka;
	VMAX_RANGE = vmax;
	AVG_TOTAL = total;
}

void setSizeForProteinInteractionNetwork(int s0, int s1, int r0, int r1)
{
	double d;
	MIN_NUM_SPECIES = s0;
	MAX_NUM_SPECIES = s1;
	if (MIN_NUM_SPECIES < 2) MIN_NUM_SPECIES = 2;
	if (MAX_NUM_SPECIES < MIN_NUM_SPECIES)
	{
		d = MIN_NUM_SPECIES;
		MIN_NUM_SPECIES = MAX_NUM_SPECIES;
		MAX_NUM_SPECIES = d;
	}
	
	if (r0 < 1) r0 = 1;
	if (r1 < r0) r1 = r0;
	MIN_NUM_REGULATIONS = (double)r0/((double)(s0 + s1)/2.0);
	MAX_NUM_REGULATIONS = (double)r1/((double)(s0 + s1)/2.0);
}

void setMutationRatesForProteinInteractionNetwork(double a, double b, double c, double d)
{
	double total;

	total = a+b+c+d;
	MUTATE_REWIRE = a/total;
	MUTATE_CHANGE_PARAM = b/total;
	MUTATE_TOTAL_CONC = c/total;
}

void setCrossoverRateForProteinInteractionNetwork(double e)
{
	CROSSOVER_PROB = e;
}

/********************************************************
    Clone, delete, mutate, crossover  (required by GA)
*********************************************************/

void deleteProteinInteractionNetwork(GAindividual individual)
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

GAindividual cloneProteinInteractionNetwork(GAindividual individual)
{
	int i,j,m,n;
	ProteinInteractionNetwork * net, * net2;
	
	if (!individual) return 0;
	
	net = (ProteinInteractionNetwork*)(individual);   //original
	net2 = (ProteinInteractionNetwork*) malloc(sizeof(ProteinInteractionNetwork)); //soon to be clone
	
	n = net->species;    //number of species
	net2->species = n;
	
	net2->regulators = (Regulators*) malloc(n * sizeof(Regulators));   //allocate space
	net2->totals = (double*) malloc(n * sizeof(double));
	net2->fixed = (int*) malloc(n * sizeof(int));
	
	for (i=0; i < n; ++i)   //copy regulators
	{
		net2->fixed[i] = net->fixed[i];
		net2->totals[i] = net->totals[i];
		m = net->regulators[i].size;
		net2->regulators[i].size = m;
		net2->regulators[i].proteins = (int*) malloc(m * sizeof(int));
		net2->regulators[i].Km = (double*) malloc(m * sizeof(double));
		net2->regulators[i].Vmax = (double*) malloc(m * sizeof(double));
		for (j=0; j < m; ++j)  //copy values for each regulator
		{
			net2->regulators[i].proteins[j] = net->regulators[i].proteins[j];
			net2->regulators[i].Km[j] = net->regulators[i].Km[j];
			net2->regulators[i].Vmax[j] = net->regulators[i].Vmax[j];
		}
	}
	
	return (GAindividual)(net2);  //done
}

GAindividual crossoverProteinInteractionNetwork(GAindividual individualA, GAindividual individualB)
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
	
	net3 = (ProteinInteractionNetwork*) malloc(sizeof(ProteinInteractionNetwork));  //child network
	net3->species = n;
	net3->regulators = (Regulators*) malloc(n * sizeof(Regulators));
	net3->totals = (double*) malloc(n * sizeof(double));
	net3->fixed = (int*) malloc(n * sizeof(int));
	
	for (i=0; i < i1; ++i) //copy all regulators from net1
	{
		n = net1->regulators[i].size;
		net3->regulators[i].size = n;
		net3->totals[i] = net1->totals[i];
		net3->fixed[i] = net1->fixed[i];
		net3->regulators[i].proteins = (int*) malloc(n * sizeof(int));
		net3->regulators[i].Km = (double*) malloc(n * sizeof(double));
		net3->regulators[i].Vmax = (double*) malloc(n * sizeof(double));
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
		net3->regulators[k].proteins = (int*) malloc(n * sizeof(int));
		net3->regulators[k].Km = (double*) malloc(n * sizeof(double));
		net3->regulators[k].Vmax = (double*) malloc(n * sizeof(double));
		for (j=0; j < n; ++j)
		{
			net3->regulators[k].proteins[j] = net2->regulators[i].proteins[j];
			if (net3->regulators[k].proteins[j] >= net3->species)
				net3->regulators[k].proteins[j] = (int)(mtrand() * net3->species);
			net3->regulators[k].Km[j] = net2->regulators[i].Km[j];
			net3->regulators[k].Vmax[j] = net2->regulators[i].Vmax[j];
		}
	}
	return (GAindividual)(net3);
}

GAindividual mutateProteinInteractionNetwork(GAindividual individual)
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
		
		return (GAindividual)(net);
	}
	if (r < (MUTATE_REWIRE+MUTATE_CHANGE_PARAM+MUTATE_TOTAL_CONC))  //change total concentrations (conservation rule)
	{
		net->totals[i] *= (2.0 * mtrand());
		return (GAindividual)(net);
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

			net->regulators = (Regulators*)malloc( n2 * sizeof(Regulators) );
			net->totals = (double*)malloc( n2 * sizeof(double) );
			net->fixed = (int*)malloc( n2 * sizeof(int) );
			
			for (j=0,j2=0; j < n; ++j) //copy all proteins
			{
				m = regulators[j].size;

				if (j != i) //except protein i
				{
					net->totals[j2] = totals[j];
					net->fixed[j2] = fixed[j];
					net->regulators[j2].size = m;
					net->regulators[j2].proteins = (int*)malloc(m * sizeof(int));
					net->regulators[j2].Km = (double*)malloc(m * sizeof(double));
					net->regulators[j2].Vmax = (double*)malloc(m * sizeof(double));
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

			return (GAindividual)(net);
		}
		else	//add a protein
		{
			n2 = n+1;
			regulators = net->regulators;
			totals = net->totals;
			fixed = net->fixed;

			net->species = n2;
			net->regulators = (Regulators*) malloc( n2 * sizeof(Regulators) );
			net->totals = (double*) malloc( n2 * sizeof(double) );
			net->fixed = (int*) malloc( n2 * sizeof(int) );
			
			for (j=0; j < n; ++j) //copy all proteins
			{
				net->totals[j] = totals[j];
				net->fixed[j] = fixed[j];
				m = regulators[j].size;
				net->regulators[j].size = m;
				net->regulators[j].proteins = (int*) malloc(m * sizeof(int));
				net->regulators[j].Vmax = (double*) malloc(m * sizeof(double));
				net->regulators[j].Km = (double*) malloc(m * sizeof(double));
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
			net->regulators[n].proteins = (int*) malloc(m * sizeof(int));
			net->regulators[n].Km = (double*) malloc(m * sizeof(double));
			net->regulators[n].Vmax = (double*) malloc(m * sizeof(double));
			for (j=0; j < m; ++j)  //random values for the new protein
			{
				net->regulators[n].proteins[j] = (int)(mtrand() * net->species);
				net->regulators[n].Km[j] = mtrand() * KM_RANGE;
				net->regulators[n].Vmax[j] = mtrand() * VMAX_RANGE;
				
				if (mtrand() < 0.5)
					net->regulators[n].Vmax[j] *= -1.0;  //negative regulator
			}

			return (GAindividual)(net);
		}
	}

    return (net);
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
	return (2 * net->species);
}

void setFixedSpeciesForProteinInteractionNetwork(GAindividual individual, int i, int value)
{
	ProteinInteractionNetwork * net = (ProteinInteractionNetwork*)(individual);
	if (i < net->species)
		net->fixed[i] = value;
}

void ratesForProteinInteractionNetwork(double time,double* u,double* rate,GAindividual individual)
{
	int i,j,n,p,forward, backward;
	double km,vmax,tot;
	ProteinInteractionNetwork * net;
	
	net = (ProteinInteractionNetwork*)(individual);
	n = net->species;
	
	//this is cheating!...but safeguards against bad initial conditions
	//should be a one-time issue, generally
	for (i=0; i < n; ++i)
	{
		if (u[i] > net->totals[i])
			net->totals[i] += u[i];
	}
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


double * stoichiometryForProteinInteractionNetwork(GAindividual p)
{
	int i,n;
	double * N;
	ProteinInteractionNetwork * net;
	
	net = (ProteinInteractionNetwork*)(p);
	n = net->species;
	N = (double*) malloc(n * 2 * n * sizeof(double));
	for (i=0; i < (2*n*n); ++i)
	{
		N[i] = 0.0;
	}
	for (i=0; i < n; ++i)
	{
		if (net->fixed[i] == 0)
		{
			getValue(N,(2*n),i,i*2) = -1.0;
			getValue(N,(2*n),i,i*2+1) = 1.0;
		}
	}
	return N;
}

void printProteinInteractionNetwork(FILE* stream, GAindividual individual)
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
		fprintf(stream,"const s%i",fix);
		for (i=0; i < net->species; ++i)
		{
			if (net->fixed[i])			
				fprintf(stream,", s%i",i+1);			
		}
		fprintf(stream,"\n");
	}
	
	for (i=0; i < n; ++i)
	{
		fprintf(stream,"$X -> s%i; ", i+1);
		f = 0;
		for (j=0; j < net->regulators[i].size; ++j)
		{
			p = net->regulators[i].proteins[j]; //index of regulating protein
			vmax = net->regulators[i].Vmax[j]; //vmax for this regulation

			if (vmax > 0)
				if (!f)
				{
					fprintf(stream,"vmax%i%i * s%i * (tot%i - s%i) / (km%i%i + (tot%i - s%i)) ",i+1,j+1,i+1,i+1,p+1,i+1,j+1,i+1,i+1);
					f = 1;
				}
				else
					fprintf(stream," + %lf * s%i * (tot%i - s%i) / (km%i%i + (tot%i - s%i)) ",i+1,j+1,i+1,i+1,p+1,i+1,j+1,i+1,i+1);
		}
		
		for (j=0; j < net->regulators[i].size; ++j)
		{
			vmax = net->regulators[i].Vmax[j]; //vmax for this regulation
			p = net->regulators[i].proteins[j]; //index of regulating protein
			if (vmax < 0) 
			{
				fprintf(stream," - vmax%i%i * s%i * s%i / (km%i%i + s%i)",i+1,j+1,p+1,i+1,i+1,j+1,i+1);
				if (!f) f = 1;
			}
		}
		
		if (!f) 
			fprintf(stream,"0;\n");
		else
			fprintf(stream,";\n");
	}
	
	fprintf(stream,"\n");
	
	for (i=0; i < n; ++i)
	{
		tot = net->totals[i]; //total concentration of u[i]
		fprintf(stream,"tot%i = %lf;\n",i,tot);
		
		for (j=0; j < net->regulators[i].size; ++j)
		{
			vmax = net->regulators[i].Vmax[j]; //vmax for this regulation
			km = net->regulators[i].Km[j];  //km for this regulation
			if (vmax != 0)
				fprintf(stream,"vmax%i%i = %lf;\n ",i+1,j+1,vmax);
			if (km != 0)
				fprintf(stream,"km%i%i = %lf;\n ",i+1,j+1,km);	
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
	
	array = (ProteinInteractionNetwork**) malloc(num * sizeof(ProteinInteractionNetwork*));
	for (i=0; i < num; ++i)
	{
		n = (int)(1 + s * 2.0 * mtrand());
		net = (ProteinInteractionNetwork*) malloc(sizeof(ProteinInteractionNetwork)); //new network
	
		net->species = n;    //number of proteins
		net->fixed = (int*) malloc(n * sizeof(int));
		
		net->regulators = (Regulators*) malloc(n * sizeof(Regulators));   //allocate space
		net->totals = (double*) malloc(n * sizeof(double));
		
		for (j=0; j < n; ++j)   //random regulators for each protein
		{
			net->fixed[j] = 0; //no fixed species by default
			net->totals[j] = 2.0*mtrand()*AVG_TOTAL;
			m = (int)(2 + mtrand() * n * AVG_NUM_REGULATIONS);
			net->regulators[j].size = m;
			net->regulators[j].proteins = (int*) malloc(m * sizeof(int));
			net->regulators[j].Km = (double*) malloc(m * sizeof(double));
			net->regulators[j].Vmax = (double*) malloc(m * sizeof(double));
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
