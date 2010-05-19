#include <iostream>
#include <vector>
#include <time.h>
#include "sbml_sim.h"
extern "C"
{
	#include "mtrand.h"
	#include "opt.h"
}

using namespace std;
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
	
	time_t seconds;
	seconds = time(NULL);
	
//	NelderMeadSimplexMethod(sim.getVariableNames().size(), diff , x0 , 10 , &fopt, 10000, 1E-3);

	vector <vector<double> > results = sim.optimize(actual,100);

	for (int i=0; i < results.size(); ++i)
	{
		for (int j=0; j < results[i].size(); ++j)
			fprintf(file4, "%lf\t", results[i][j]);
		fprintf(file4,"\n");
	}
	
	sim.reset();	
	params = results[ results.size() - 1];
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
	
	sim.reset();
	params = results[0];
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


