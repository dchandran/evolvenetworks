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

int main()
{

/************************ LOAD SBML MODEL ***********************************/
	string sbml("feedback.xml");
	SBML_sim sim(sbml);

/************************** OUTPUT FILES ************************************/
	FILE * file1 = fopen("out1.tab","w");  //actual output
	FILE * file2 = fopen("out2.tab","w");  //best optimized output
	FILE * file3 = fopen("out3.tab","w");  //worst optimized output
	FILE * file4 = fopen("soln.tab","w");   //range of optimized values
	
	initMTrand(); //random number generator
	
/*********** GET ACTUAL OUTPUT OF THE ORIGINAL MODEL ************************/

	double end_time = 20.0, dt = 0.1;
	
	vector< vector<double> > actual = sim.simulate(end_time, dt);

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
	
/***********************  OPTIMIZATION ***************************************/

	time_t seconds;    //record time
	seconds = time(NULL);	
	
	vector< vector<double> > soln = sim.optimize(actual,100);  //optimize for 100 iterations
	
	printf("time = %i seconds\n", (int)(time(NULL)-seconds));  //print time
	
/*********************  PRINT OPTIMIZED PARAMETERS ***************************/

	for (int i=0; i < soln.size(); ++i)
	{
		for (int j=0; j < soln[i].size(); ++j)
			fprintf(file4, "%lf\t", soln[i][j]);
		fprintf(file4, "\n");
	}

/*********************  SIMULATE USING BEST PARAMETERS **********************/

	sim.setParameters(soln[0]);	  //index 0 = best
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

/*********************  SIMULATE USING WORST PARAMETERS *********************/

	sim.setParameters(soln[ soln.size()-1 ]);	//worst
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

//close files
	fclose(file1);
	fclose(file2);
	fclose(file3);
	fclose(file4);
	
	return 0;
}


