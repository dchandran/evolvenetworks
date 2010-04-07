#include <vector>
#include "sbml_sim.h"
extern "C"
{
	#include "opt.h"
}

using namespace std;

vector< vector<double> > actual;
SBML_sim * model;
double time = 20, dt = 0.1;

double diff(int n, double * p)
{
	vector<double> params(n,0);
	for (int i=0; i < n; ++i) params[i] = p[i];
	
	vector< vector<double> > res = model->simulate(time, dt);
	
	double sumsq = 0.0, d = 0.0;
	
	for (int i=1; i < res.size(); ++i)
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
	
	initMTrand();
	
	actual = sim.simulate(time, dt);
	model = &sim;
	
	int n = res.size();
	if (n > 0)
	{
		int m = res[0].size();
		for (int i=0; i < m; ++i)
		{
			for (int j=0; j < n; ++j)
				printf("%lf\t",res[j][i]);
			printf("\n");
		}
	}
	
	NelderMeadSimplexMethod(
	
	return 0;
}


