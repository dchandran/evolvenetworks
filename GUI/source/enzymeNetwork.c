/*******************************************************

	Copyright (C) 2009 Deepak Chandran
	see header file

********************************************************/

#include "enzymeNetwork.h"

/************************
    global variables
*************************/

static double MUTATE_KEQ_PROB = 0.1;
static double KEQ_LN_RANGE = 4.0;

static double MUTATE_ALPHA_PROB = 0.1;
static double ALPHA_LN_RANGE = 2.0;

static double MUTATE_H_PROB = 0.1;
static double H_RANGE = 2.0;

static double MUTATE_S_HALF_PROB = 0.1;
static double S_HALF_RANGE = 20.0;

static double MUTATE_P_HALF_PROB = 0.1;
static double P_HALF_RANGE = 20.0;

static double MUTATE_ENZYME_PROB = 0.1;
static double CROSSOVER_PROB = 0.2;

void setDistributionOfEnzymeNetwork(double uni_uni, double uni_bi, double bi_uni, double bi_bi, double no_reactant, double no_product)
{
	setDistributionOfMassActionNetwork(uni_uni, uni_bi, bi_uni, bi_bi, no_reactant, no_product);
}

void setRateConstantsForEnzymeNetwork(double max_kcat, double max_log_keq, double max_alpha, double max_h, double max_s_half, double max_p_half)
{
	setRateConstantForMassActionNetwork(max_kcat);
	if (max_log_keq > 0)
		KEQ_LN_RANGE = max_log_keq;
	if (max_alpha > 0)
		ALPHA_LN_RANGE = max_alpha;
	if (max_s_half > 0)
		S_HALF_RANGE = max_s_half;
	if (max_p_half > 0)
		P_HALF_RANGE = max_p_half;
	if (max_h > 2)
		H_RANGE = max_h-2.0;
}

void setSizeForEnzymeNetwork(int n, int s)
{
	setSizeForMassActionNetwork(n,s);
}

void setMutationRatesForEnzymeNetwork(double enzyme, 
									  double k_cat, 
									  double k_eq, 
									  double alpha,
									  double h, 
									  double s_half, 
									  double p_half, 
									  double remove, 
									  double add)
{
	double total;

	if (k_cat < 0) k_cat = 0;
	if (k_eq < 0) k_eq = 0;
	if (alpha < 0) alpha = 0;
	if (s_half < 0) s_half = 0;
	if (p_half < 0) p_half = 0;
	if (h < 0) h = 0;
	if (remove < 0) remove = 0;
	if (add < 0) add = 0;
	if (enzyme < 0) enzyme = 0;

	total = k_cat + h + k_eq + s_half + p_half + remove + add + enzyme;

	if (total == 0) return;
	
	alpha /= total;
	k_cat /= total;
	k_eq /= total; 
	s_half /= total; 
	p_half /= total; 
	h /= total; 
	remove /= total; 
	add /= total; 
	enzyme /= total;

	setMutationRatesForMassActionNetwork(k_cat, remove, add);

	MUTATE_KEQ_PROB = k_eq;
	MUTATE_H_PROB = h;
	MUTATE_ALPHA_PROB = alpha;
	MUTATE_S_HALF_PROB = s_half;
	MUTATE_P_HALF_PROB = p_half;
	MUTATE_ENZYME_PROB = enzyme;
}

void setCrossoverRateForEnzymeNetwork(double crossover)
{
	CROSSOVER_PROB = crossover;
}

/********************************************************
    Clone, delete, mutate, crossover  (required by GA)
*********************************************************/

void deleteEnzymeNetwork(GAindividual individual)
{
	EnzymeNetwork * net;
	
	if (!individual) return;
	
	net = (EnzymeNetwork*)(individual);
	
	if (net->massActionNetwork)
	{
		deleteMassActionNetwork(net->massActionNetwork);
	}

	if (net->enzymes)
		free (net->enzymes);
	
	if (net->Keq)
		free(net->Keq);

	if (net->h)
		free(net->h);

	if (net->alpha)
		free(net->alpha);
	
	if (net->P_half)
		free(net->P_half);

	if (net->S_half)
		free(net->S_half);

	free(net);
}

GAindividual cloneEnzymeNetwork(GAindividual individual)
{
	int i,n;
	EnzymeNetwork * net, * net2;
	
	if (!individual) return 0;
	
	net = (EnzymeNetwork*)(individual);   //original
	net2 = (EnzymeNetwork*) malloc(sizeof(EnzymeNetwork)); //soon to be clone
	
	net2->massActionNetwork = (MassActionNetwork*) cloneMassActionNetwork(net->massActionNetwork);

	n = net->massActionNetwork->reactions;
	net2->enzymes = (int*) malloc (n * sizeof(int));
	net2->Keq = (double*) malloc (n * sizeof(double));
	net2->h = (double*) malloc (n * sizeof(double));
	net2->alpha = (double*) malloc (n * sizeof(double));
	net2->P_half = (double*) malloc (n * sizeof(double));
	net2->S_half = (double*) malloc (n * sizeof(double));

	for (i=0; i < n; ++i)
	{
		net2->enzymes[i] = net->enzymes[i];
		net2->Keq[i] = net->Keq[i];
		net2->h[i] = net->h[i];
		net2->alpha[i] = net->alpha[i];
		net2->P_half[i] = net->P_half[i];
		net2->S_half[i] = net->S_half[i];
	}

	return (GAindividual)(net2);  //done
}

GAindividual crossoverEnzymeNetwork(GAindividual individualA, GAindividual individualB)  //crossover between complexes in two networks
{
	int i,j, n;
	EnzymeNetwork * net1, * net2, * net3;
	
	if (mtrand() > CROSSOVER_PROB) return mutateEnzymeNetwork(cloneEnzymeNetwork(individualA));
	
	if (!individualA) return mutateEnzymeNetwork(cloneEnzymeNetwork(individualB));
	if (!individualB) return mutateEnzymeNetwork(cloneEnzymeNetwork(individualA));
	
	net1 = (EnzymeNetwork*)(individualA);  //parents
	net2 = (EnzymeNetwork*)(individualB);
	
	if (net1->massActionNetwork->reactions < 3) return mutateEnzymeNetwork(cloneEnzymeNetwork(net2));  //if parents are too small
	if (net2->massActionNetwork->reactions < 3) return mutateEnzymeNetwork(cloneEnzymeNetwork(net1));
	
	net3 =  (EnzymeNetwork*) malloc (sizeof(EnzymeNetwork));

	net3->massActionNetwork = crossoverMassActionNetwork(net1->massActionNetwork,net2->massActionNetwork);

	n = net3->massActionNetwork->reactions;

	net3->enzymes = (int*)malloc(n * sizeof(int));
	net3->Keq = (double*)malloc(n * sizeof(double));
	net3->h = (double*)malloc(n * sizeof(double));
	net3->alpha = (double*)malloc(n * sizeof(double));
	net3->P_half = (double*)malloc(n * sizeof(double));
	net3->S_half = (double*)malloc(n * sizeof(double));

	//get some set of enzymes from one parent and some from the other

	for (i=0; i < n; ++i)
	{
		net3->enzymes[i] = (int)(mtrand() * net3->massActionNetwork->species);
		net3->Keq[i] = mtrand() * pow(2, (KEQ_LN_RANGE * (2.0 * mtrand() - 1)));
		net3->h[i] = 2.0 + H_RANGE * mtrand();
		net3->alpha[i] = mtrand() * pow(2, (ALPHA_LN_RANGE * (2.0 * mtrand() - 1)));
		net3->S_half[i] = S_HALF_RANGE * mtrand();
		net3->P_half[i] = P_HALF_RANGE * mtrand();
	}

	j = (int)(mtrand() * net3->massActionNetwork->reactions);
	for (i=0; i < j && i < net1->massActionNetwork->reactions && i < net3->massActionNetwork->reactions; ++i)
	{
		net3->enzymes[i] = net1->enzymes[i];
		net3->Keq[i] = net1->Keq[i];
		net3->h[i] = net1->h[i];
		net3->alpha[i] = net1->alpha[i];
		net3->S_half[i] = net1->S_half[i];
		net3->P_half[i] = net1->P_half[i];

		if (net3->enzymes[i] >= net3->massActionNetwork->species)
			net3->enzymes[i] = (int)(mtrand() * net3->massActionNetwork->species);
	}
	for (i=0; (i+j) < net3->massActionNetwork->reactions && i < net2->massActionNetwork->reactions; ++i)
	{
		net3->enzymes[i+j] = net2->enzymes[i];
		net3->Keq[i+j] = net2->Keq[i];
		net3->h[i+j] = net2->h[i];
		net3->alpha[i+j] = net2->alpha[i];
		net3->S_half[i+j] = net2->S_half[i];
		net3->P_half[i+j] = net2->P_half[i];
		if (net3->enzymes[i+j] >= net3->massActionNetwork->species)
			net3->enzymes[i+j] = (int)(mtrand() * net3->massActionNetwork->species);
	}
	return (GAindividual)(net3);
}

GAindividual mutateEnzymeNetwork(GAindividual individual)
{
	int i,j,n;
	double r;
	EnzymeNetwork * net, *net2;
	MassActionNetwork * mnet;

	if (!individual) 
		return individual;
	
	net = (EnzymeNetwork*)individual;
	mnet = net->massActionNetwork;

	if (!mnet) 
		return individual;
	
	r = mtrand();

	n = 0;

	for (i=0; i < mnet->reactions; ++i)
	{
		if (((mnet->reactant1[i] == -1 && mnet->reactant2[i] > -1) ||  //uni-uni
			 (mnet->reactant1[i] > -1 && mnet->reactant2[i] == -1))
			&&
			((mnet->product1[i] == -1 && mnet->product2[i] > -1) ||
			 (mnet->product1[i] > -1 && mnet->product2[i] == -1))
			)
			++n;
	}

	if (n > 0)
	{
		j = (int)(mtrand() * n);

		for (i=0; i < mnet->reactions; ++i)
		{
			if (((mnet->reactant1[i] == -1 && mnet->reactant2[i] > -1) ||  //uni-uni
				 (mnet->reactant1[i] > -1 && mnet->reactant2[i] == -1))
				&&
				((mnet->product1[i] == -1 && mnet->product2[i] > -1) ||
				 (mnet->product1[i] > -1 && mnet->product2[i] == -1))
				)
				--j;
			if (j < 0)
				break;
		}

		if (i >= mnet->reactions)
			return individual;

		//i = (int)(mtrand() * mnet->reactions);

		if (r < MUTATE_KEQ_PROB) //mutate km
		{
			net->Keq[i] *= 2.0 * mtrand();
			return (GAindividual)(net);
		}
		if (r < (MUTATE_KEQ_PROB+MUTATE_S_HALF_PROB))   //mutate s_half
		{
			net->S_half[i] *= 2.0 * mtrand();
			return (GAindividual)(net);
		}
		if (r < (MUTATE_KEQ_PROB+MUTATE_S_HALF_PROB+MUTATE_P_HALF_PROB))   //mutate s_half
		{
			net->P_half[i] *= 2.0 * mtrand();
			return (GAindividual)(net);
		}
		if (r < (MUTATE_KEQ_PROB+MUTATE_S_HALF_PROB+MUTATE_P_HALF_PROB+MUTATE_ENZYME_PROB))   //mutate enzyme
		{
			net->enzymes[i] = (int)(mtrand() * net->massActionNetwork->species);
			return (GAindividual)(net);
		}
		if (r < (MUTATE_KEQ_PROB+MUTATE_S_HALF_PROB+MUTATE_P_HALF_PROB+MUTATE_ENZYME_PROB+MUTATE_H_PROB))   //mutate h
		{
			net->h[i] *= 2.0 * mtrand();
			return (GAindividual)(net);
		}
		if (r < (MUTATE_H_PROB+MUTATE_KEQ_PROB+MUTATE_S_HALF_PROB+MUTATE_P_HALF_PROB+MUTATE_ENZYME_PROB+MUTATE_ALPHA_PROB)) //mutate alpha
		{
			net->alpha[i] *= 2.0 * mtrand();
			return (GAindividual)(net);
		}
	}
	
	n = net->massActionNetwork->reactions;

	mnet = (MassActionNetwork*)mutateMassActionNetwork(net->massActionNetwork);

	//mnet = net->massActionNetwork;
	
	if (n == mnet->reactions)
	{
		net->massActionNetwork = mnet;
		net2 = net;
	}
	else
	{
		net->massActionNetwork = 0;

		net2 = (EnzymeNetwork*)malloc(sizeof(EnzymeNetwork));

		net2->massActionNetwork = mnet;
		net2->enzymes = (int*)malloc(mnet->reactions * sizeof(int));
		net2->Keq = (double*)malloc(mnet->reactions * sizeof(double));
		net2->h = (double*)malloc(mnet->reactions * sizeof(double));
		net2->alpha = (double*)malloc(mnet->reactions * sizeof(double));
		net2->P_half = (double*)malloc(mnet->reactions * sizeof(double));
		net2->S_half = (double*)malloc(mnet->reactions * sizeof(double));

		for (i=0; i < mnet->reactions; ++i)
		{
			net2->enzymes[i] = (int)(mtrand() * mnet->species);
			net2->Keq[i] = mtrand() * pow(2, (KEQ_LN_RANGE * (2.0 * mtrand() - 1)));
			net2->h[i] = 2.0 + (mtrand() * H_RANGE);
			net2->alpha[i] = mtrand() * pow(2, (ALPHA_LN_RANGE * (2.0 * mtrand() - 1)));
			net2->P_half[i] = (mtrand() * P_HALF_RANGE);
			net2->S_half[i] = (mtrand() * S_HALF_RANGE);
		}

		for (i=0; i < n && i < mnet->reactions; ++i)
		{
			net2->enzymes[i] = net->enzymes[i];
			net2->Keq[i] = net->Keq[i];
			net2->h[i] = net->h[i];
			net2->alpha[i] = net->alpha[i];
			net2->P_half[i] = net->P_half[i];
			net2->S_half[i] = net->S_half[i];
		}

		deleteEnzymeNetwork(net);
	}

	return (GAindividual)(net2);
}

/*****************************************************
   Functions for simulating and printing
******************************************************/

int getNumSpeciesForEnzymeNetwork(GAindividual individual)
{
	EnzymeNetwork * net = (EnzymeNetwork*)(individual);
	if (!net) return 0;
	return (getNumSpeciesForMassActionNetwork(net->massActionNetwork));
}

int getNumReactionsForEnzymeNetwork(GAindividual individual)
{
	EnzymeNetwork * net = (EnzymeNetwork*)(individual);
	if (!net) return 0;
	return (getNumReactionsForMassActionNetwork(net->massActionNetwork));
}

void setFixedSpeciesForEnzymeNetwork(GAindividual individual, int i, int value)
{
	EnzymeNetwork * net = (EnzymeNetwork*)(individual);
	if (!net) return;
	setFixedSpeciesForMassActionNetwork(net->massActionNetwork,i,value);
}

void ratesForEnzymeNetwork(double time,double* u,double* rate,GAindividual individual)
{
	int i,uni;
	double s1,s2,p1,p2,e;
	EnzymeNetwork * enet;
	MassActionNetwork * net;
	
	enet = (EnzymeNetwork*)(individual);
	
	if (!enet) return;

	net = enet->massActionNetwork;

	if (!net) return;
	
	for (i=0; i < net->reactions; ++i)
	{
		s1 = s2 = p1 = p2 = -1.0;
		if (net->reactant1[i] > -1)
			s1 = u[ net->reactant1[i] ];
		if (net->reactant2[i] > -1)
			s2 = u[ net->reactant2[i] ];
		if (net->product1[i] > -1)
			p1 = u[ net->product1[i] ];
		if (net->product2[i] > -1)
			p2 = u[ net->product2[i] ];
		uni =(((net->reactant1[i] == -1 && net->reactant2[i] > -1) ||  //uni-uni
			(net->reactant1[i] > -1 && net->reactant2[i] == -1))
			&&
			((net->product1[i] == -1 && net->product2[i] > -1) ||
			(net->product1[i] > -1 && net->product2[i] == -1)));
		
		rate[i] = net->k[i];
		if (uni && enet->enzymes[i] > -1)
		{
			e = u[ enet->enzymes[i] ];
			if (enet->h[i] < 2.0) 
				enet->h[i] = 2.0;
			rate[i] *= s1/enet->S_half[i] * 
				(1 - (p1/s1)/enet->Keq[i]) * pow((s1/enet->S_half[i] + p1/enet->P_half[i]),enet->h[i]-1) *
				1.0/ ((1 + e)/(1 + enet->alpha[i]*e) + pow((s1/enet->S_half[i] + p1/enet->P_half[i]),enet->h[i]));
		}
		else
		{
			if (net->reactant1[i] > -1)
			{
				rate[i] *= u[ net->reactant1[i] ];
			}
			if (net->reactant2[i] > -1) 
			{
				rate[i] *= u[ net->reactant2[i] ];
			}
		}
	}
}

double * stoichiometryForEnzymeNetwork(GAindividual individual)
{
	EnzymeNetwork * net = (EnzymeNetwork*)(individual);
	if (!net) return 0;
	return (stoichiometryForMassActionNetwork(net->massActionNetwork));
}

void printEnzymeNetwork(FILE * stream,GAindividual individual)
{
	int i,fix, r, p;
	MassActionNetwork * net;
	EnzymeNetwork * enet;
	
	if (!individual) return;
	enet = (EnzymeNetwork*)individual;
	net = enet->massActionNetwork;

	if (!net) return;
	
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
		if (net->reactant1[i] > -1 && net->reactant2[i] > -1)
			fprintf(stream, "s%i + s%i -> ",net->reactant1[i]+1,net->reactant2[i]+1);
		else
		{
			if (net->reactant1[i] > -1 || net->reactant2[i] > -1)
			{
				if (net->reactant1[i] > -1)
					fprintf(stream, "s%i ->",net->reactant1[i]+1);
				else
					fprintf(stream, "s%i ->",net->reactant2[i]+1);
			}
			else
			{
				fprintf(stream, "$x -> ");
			}
		}

		if (net->product1[i] > -1 && net->product2[i] > -1)
			fprintf(stream, "s%i + s%i; ",net->product1[i]+1,net->product2[i]+1);
		else
		{
			if (net->product1[i] > -1 || net->product2[i] > -1)
			{
				if (net->product1[i] > -1)
					fprintf(stream, "s%i; ",net->product1[i]+1);
				else
					fprintf(stream, "s%i; ",net->product2[i]+1);
			}
			else
			{
				fprintf(stream, "$x; ");
			}
		}
		
		//rate
		if (net->reactant1[i] > -1 && net->reactant2[i] > -1)
			fprintf(stream, "k%i * s%i * s%i ",i+1,net->reactant1[i]+1,net->reactant2[i]+1);
		else
		{
			if (net->reactant1[i] > -1 || net->reactant2[i] > -1)
			{
				if (net->reactant1[i] > -1)
				{
					r = net->reactant1[i]+1;
				}
				else
				{
					r = net->reactant2[i]+1;
				}
				if ((net->product1[i] == -1 && net->product2[i] > -1) ||  //uni-uni
					(net->product1[i] > -1 && net->product2[i] == -1))
				{
					if (net->product1[i] == -1)
						p = net->product2[i]+1;
					else
						p = net->product1[i]+1;
					fprintf(stream, "k%i * s%i/shalf%i * (1 - (s%i/s%i)/keq%i) * (s%i/shalf%i + s%i/phalf%i)^(h%i-1) / ((s%i/shalf%i + s%i/phalf%i)^(h%i) + (1+s%i)/(1+alpha%i*s%i))",
						i+1,r,i+1,r,p,i+1,r,i+1,p,i+1,enet->h[i],r,i+1,p,i+1,enet->h[i],enet->enzymes[i]+1,i+1,enet->enzymes[i]+1);
				}
				else
				{
					fprintf(stream, "k%i * s%i",i+1,r);
				}
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
		if (((net->reactant1[i] == -1 && net->reactant2[i] > -1) ||  //uni-uni
			(net->reactant1[i] > -1 && net->reactant2[i] == -1))
			&&
			((net->product1[i] == -1 && net->product2[i] > -1) ||
			(net->product1[i] > -1 && net->product2[i] == -1)))
		{
			fprintf(stream, "keq%i = %lf;\n",i+1,enet->Keq[i]);
			fprintf(stream, "h%i = %lf;\n",i+1,enet->h[i]);
			fprintf(stream, "alpha%i = %lf;\n",i+1,enet->alpha[i]);
			fprintf(stream, "shalf%i = %lf;\n",i+1,enet->S_half[i]);
			fprintf(stream, "phalf%i = %lf;\n",i+1,enet->P_half[i]);
		}
	}
}

/***********************
  GA related functions
***********************/

GApopulation randomEnzymeNetworks(int num)
{
	int i,j;
	MassActionNetwork * mnet;
	EnzymeNetwork * enet;
	GApopulation pop = randomMassActionNetworks(num);

	for (i=0; i < num; ++i)
	{
		mnet = (MassActionNetwork*)pop[i];
		enet = (EnzymeNetwork*) malloc(sizeof(EnzymeNetwork));
		enet->enzymes = (int*) malloc( mnet->reactions * sizeof(int) );
		enet->Keq = (double*) malloc( mnet->reactions * sizeof(double) );
		enet->h = (double*) malloc( mnet->reactions * sizeof(double) );
		enet->alpha = (double*) malloc( mnet->reactions * sizeof(double) );
		enet->S_half = (double*) malloc( mnet->reactions * sizeof(double) );
		enet->P_half = (double*) malloc( mnet->reactions * sizeof(double) );
		for (j=0; j < mnet->reactions; ++j)
		{
			enet->enzymes[j] = (int)(mtrand() * mnet->species);
			enet->Keq[j] = mtrand() * pow(2, (KEQ_LN_RANGE * (2.0 * mtrand() - 1)));
			enet->alpha[j] = mtrand() * pow(2, (ALPHA_LN_RANGE * (2.0 * mtrand() - 1)));
			enet->h[j] = H_RANGE * mtrand();
			enet->S_half[j] = S_HALF_RANGE * mtrand();
			enet->P_half[j] = P_HALF_RANGE * mtrand();
		}
		enet->massActionNetwork = mnet;
		pop[i] = enet;
	}
	
	return pop;
}

EnzymeNetwork * newEnzymeNetwork(int m,int n)
{
	int i;
	EnzymeNetwork * net;

	net = (EnzymeNetwork*) malloc(sizeof(EnzymeNetwork));
	net->massActionNetwork = newMassActionNetwork(m,n);
	net->enzymes = (int*) malloc( n * sizeof(int) );
	net->Keq = (double*) malloc( n * sizeof(double) );
	net->h = (double*) malloc( n * sizeof(double) );
	net->alpha = (double*) malloc( n * sizeof(double) );
	net->S_half = (double*) malloc( n * sizeof(double) );
	net->P_half = (double*) malloc( n * sizeof(double) );

	for (i=0; i < n; ++i)
	{
		net->enzymes[i] = 0;
		net->Keq[i] = 1.0;
		net->h[i] = 1.0;
		net->alpha[i] = 1.0;
		net->S_half[i] = 1.0;
		net->P_half[i] = 1.0;
	}
	
	return net;
}
