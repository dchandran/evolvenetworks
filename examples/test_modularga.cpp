#include <iostream>
#include <vector>
#include <time.h>
#include "GAmodular.h"
#include "sbml_sim.h"

using namespace std;

vector< vector<double> > actual;
SBML_sim * model;
double end_time = 20;
double dt = 0.1;


double objective(Genome * x)
{	
	model->setParameters(x->params);

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

	return (sumsq/ (res[0].size())); //* (res.size()-1)));

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
	
	actual = sim.simulate(end_time, dt);
	model = &sim;	
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
	
	for (int i=0; i < numParams; ++i)
		cout << paramNames[i] << "\t";
	cout << endl << endl;
	
	for (int i=0; i < numParams; ++i)
		cout << params[i] << "\t";

	cout << endl << endl;

	ModularGA ga;
	ga.setObjective( objective , ModularGA::Minimize );	
	ga.setPopulationSize(1000);
	ga.setGenerations(100);
	ga.setParameterSize(params.size());
	
	time_t seconds;
	seconds = time(NULL);


	ga.evolve();
	
	vector<Genome*> pop = ga.population();
	
	for (int i=0; i < pop.size(); ++i)
	{
		Genome * g = pop[i];
		for (int j=0; j < g->params.size(); ++j)
			fprintf(file4, "%lf\t", g->params[j]);
		fprintf(file4,"%lf\n", g->score);
	}
	
	Genome * g = ga.best();
	cout << g->score << endl;
		
	sim.setParameters(g->params);	
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
	
	g = ga.worst();
	cout << g->score << endl;
		
	sim.setParameters(g->params);
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
	return 0;
}
