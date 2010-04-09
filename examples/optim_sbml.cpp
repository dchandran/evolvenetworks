#include <iostream>
#include <vector>
#include "sbml_sim.h"
extern "C"
{
	#include "opt.h"
}

using namespace std;

vector< vector<double> > actual;
SBML_sim * model;
double end_time = 20;
double dt = 0.1;

double diff(int n, double * p)
{
	vector<double> params(n,0);
	for (int i=0; i < n; ++i) params[i] = p[i];

	model->setParameters(params);

	vector< vector<double> > res = model->simulate(end_time, dt);

	double sumsq = 0.0, d = 0.0;
	
	for (int i=1; i < 2; ++i)
	{
		for (int j=1; j < res[i].size(); ++j)
		{
			d = res[i][j] - actual[i][j];
			sumsq += d*d;
		}
	}

	return (sumsq/ (res[0].size() * (res.size()-1)));
}

int main()
{
	string sbml("feedback.xml");
	SBML_sim sim(sbml);
	FILE * file1 = fopen("out1.tab","w");
	FILE * file2 = fopen("out2.tab","w");
	
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
	int numParams = params.size();
	double * x0 = (double*)malloc(numParams * sizeof(double));
	for (int i=0; i < numParams; ++i)
		x0[i] = 100.0;
	double fopt = 0.0;
	
	for (int i=0; i < numParams; ++i)
		cout << params[i] << " ";

	cout << endl << endl;
	
	NelderMeadSimplexMethod(sim.getVariableNames().size(), diff , x0 , 10 , &fopt, 10000, 1E-3);
	
	for (int i=0; i < numParams; ++i)
		params[i] = x0[i];
	
	for (int i=0; i < numParams; ++i)
		cout << params[i] << " ";

	cout << endl << endl;
		
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
	
	fclose(file1);
	fclose(file2);
	
	delete x0;
	
	return 0;
}


