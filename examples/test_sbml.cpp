#include <iostream>
#include <vector>
#include <time.h>
#include "sbml_sim.h"

using namespace std;

int main()
{
	string sbml("test.xml");
	SBML_sim sim(sbml);
	sim.reset();
	vector< vector<double> > output = sim.simulate(100.0, 0.1);
	int n = output.size();

	//print to file
	
	FILE * file = fopen("out.tab","w");
	int m = output[0].size();
	for (int i=0; i < m; ++i)
	{
		for (int j=0; j < n; ++j)
			fprintf(file, "%lf\t",output[j][i]);
		fprintf(file, "\n");
	}
	
	fclose(file);
	return 0;
}


