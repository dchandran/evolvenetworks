#include <iostream>
#include <vector>
#include <time.h>
#include "GASimpleGA.h"
#include "GASStateGA.h"
#include "GA1DArrayGenome.h"
#include "sbml_sim.h"
extern "C"
{
	#include "mtrand.h"
	#include "opt.h"
}

using namespace std;
typedef GA1DArrayGenome<float> RealGenome;

vector< vector<double> > actual;
double end_time = 20;
double dt = 0.1;

void initializeGenome(GAGenome & x)
{
	RealGenome & g = (RealGenome &)x;
	for (int i=0; i < g.size(); ++i)
		g.gene(i,0) = mtrand() * pow(10.0, 2.0*mtrand());
}

float EuclideanDistance(const GAGenome & c1, const GAGenome & c2)
{
  const RealGenome & a = (RealGenome &)c1;
  const RealGenome & b = (RealGenome &)c2;

  float x=0.0;
  for(int i=0; i < b.length() && i < a.length(); ++i)
	  x += (a.gene(i) - b.gene(i))*(a.gene(i) - b.gene(i));

  return (float)sqrt(x);
}
/*
double diff(int n, double * p)
{
	vector<double> params(n,0);
	for (int i=0; i < n; ++i) params[i] = p[i];

	model->setParameters(params);

	vector< vector<double> > res = model->simulate(end_time, dt);
	
	if (res.size() < 1)
		return 100.0;

	double sumsq = 0.0, d = 0.0;
	
	for (int i=1; i < res.size(); ++i)
	{
		for (int j=0; j < res[i].size(); ++j)
		{
			d = res[i][j] - actual[i][j];
			sumsq += d*d;
		}
	}

	return (float)(sumsq/ (res[0].size()));// * (res.size()-1)));
}
*/
float Objective1(GAGenome & x)
{
	RealGenome & g = (RealGenome &)x;
	
	vector<double> params(g.size(),0);
	for (int i=0; i < g.size(); ++i) params[i] = g.gene(i);

	SBML_sim * model = (SBML_sim*)(g.geneticAlgorithm()->userData());
	
	model->reset();
	model->setParameters(params);

	vector< vector<double> > res = model->simulate(end_time, dt);
	
	if (res.size() < 1)
		return 100.0;

	double sumsq = 0.0, d = 0.0;
	
	for (int i=1; i < res.size(); ++i)
	{
		for (int j=0; j < res[i].size(); ++j)
		{
			d = res[i][j] - actual[i][j];
			sumsq += d*d;
		}
	}

	return (float)(sumsq/ (res[0].size()));// * (res.size()-1)));
}

int main()
{
	string sbml("feedback.xml");
	SBML_sim sim(sbml);
	FILE * file1 = fopen("out1.tab","w");
	FILE * file2 = fopen("out2.tab","w");
	FILE * file3 = fopen("out3.tab","w");
	FILE * file4 = fopen("pop.tab","w");
	
	initMTrand();
	
	sim.reset();
	actual = sim.simulate(end_time, dt);
	int n = actual.size();
	if (n > 0)
	{
		int m = actual[0].size();
		for (int i=0; i < m; ++i)
		{
			for (int j=0; j < n; ++j)
				fprintf(file1, "%lf\t",actual[j][i]);
			fprintf(file1, "\n");
		}
	}
	
	vector<double> params = sim.getParameterValues();
	vector<string> paramNames = sim.getParameterNames();
	int numParams = params.size();
	double * x0 = (double*)malloc(numParams * sizeof(double));
	
	for (int i=0; i < numParams; ++i)
		x0[i] = mtrand() * 100.0;
	double fopt = 0.0;
	
	for (int i=0; i < numParams; ++i)
		cout << paramNames[i] << "\t";
	cout << endl << endl;
	
	for (int i=0; i < numParams; ++i)
		cout << params[i] << "\t";

	cout << endl << endl;
	
//	NelderMeadSimplexMethod(sim.getVariableNames().size(), diff , x0 , 10 , &fopt, 10000, 1E-3);
	
	RealGenome genome( numParams, &Objective1);
	genome.initializer(&initializeGenome);
	GASteadyStateGA ga(genome);
	ga.userData(&sim);

	//model_data * u = (model_data*)makeModelData();
	//ga.userData(u);
	
	time_t seconds;
	seconds = time(NULL);


	GASharing dist(EuclideanDistance);
	ga.scaling(dist);
	ga.pReplacement(1.0);
	ga.minimize();
	ga.populationSize(1000);
	ga.nGenerations(100);
	ga.pMutation(0.2);
	ga.pCrossover(0.9);
	GAPopulation pop;
	ga.evolve();
	
	pop = ga.population();
	pop.order(GAPopulation::HIGH_IS_BEST);
	pop.sort(gaTrue);
	
	for (int i=0; i < pop.size(); ++i)
	{
		RealGenome & g = (RealGenome &)(pop.individual(i));
		for (int j=0; j < g.size(); ++j)
			fprintf(file4, "%lf\t", g.gene(j));
		fprintf(file4,"%lf\n", g.score());
	}
	
	RealGenome & g = (RealGenome &)(pop.individual(0));
	cout << g.score() << endl;
	
	sim.reset();	
	sim.setParameters(params);	
	actual = sim.simulate(end_time, dt);
	n = actual.size();
	
	if (n > 0)
	{
		int m = actual[0].size();
		for (int i=0; i < m; ++i)
		{
			for (int j=0; j < n; ++j)
				fprintf(file2, "%lf\t",actual[j][i]);
			fprintf(file2, "\n");
		}
	}
	
	g = (RealGenome &)(pop.individual(pop.size()-1));
	cout << g.score() << endl;
	
	sim.reset();
	sim.setParameters(params);
	actual = sim.simulate(end_time, dt);
	n = actual.size();
	
	if (n > 0)
	{
		int m = actual[0].size();
		for (int i=0; i < m; ++i)
		{
			for (int j=0; j < n; ++j)
				fprintf(file3, "%lf\t",actual[j][i]);
			fprintf(file3, "\n");
		}
	}
	
	printf("time = %i seconds\n", (int)(time(NULL)-seconds));
	
	fclose(file1);
	fclose(file2);
	fclose(file3);
	fclose(file4);
	delete x0;
	
	return 0;
}


