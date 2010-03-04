#include <math.h>
#include <stdlib.h>
#include <stdlio.h>
#include "antimony_api.h"
#include "DesignFeedforwardGRN.h"

/*number of inputs*/
static int NUM_INPUTS = 2;

/*number of grid lines for each dimension*/
static int * GRID_SIZE = 0;

/*target function that maps inputs to outputs*/
static FFN_IO_FUNC TARGET_FUNC = 0;
/*! 
* \brief set the number of inputs in the system
*/
void setNumInputs(int n)
{
	NUM_INPUTS = n;
}

/*! 
* \brief set the intput->output mapping function
*/
void setTargetFunction(FFN_IO_FUNC f)
{
	TARGET_FUNC = f;
}

/*! 
* \brief set the intput->output mapping function
*/
void setGridSizes(int * grid_sz)
{
	GRID_SIZE = grid_sz;
}

/*! 
* \brief generates the graphviz + parameters output file
* \param FILE* file for graphviz output
* \param FILE* file for list of parameter values
*/
void generateGraphFile(FILE * graphFile, FILE * paramsFile)
{
	int i,j,k,n,a,deft,on;
	double * x;
	char * green = "
	
	if (!TARGET_FUNC || !GRID_SIZE) return;
	x = (double*)malloc(NUM_INPUTS * sizeof(double)); //inputs
	
	fprintf(graphFile, "digraph G {\n");
	
	//input 0
	n = GRID_SIZE[0]; //num dividers in input 1	
	for (i=0; i < NUM_INPUTS; ++i)	x[i] = 0.0;  //initialize
	deft = TARGET_FUNC(x);  //gene is on by default?
	
	for (i=1; i < n; ++i)
	{
		x[0] = i;
		k = TARGET_FUNC(x);
		if (k != on)
		{
			fprintf(graphFile, "   x0 -> g$i; [ arrowhead=\"dot\" ]", ++a);
			if (k > 0)
				fprintf(graphFile, "   g$i -> y; [ arrowhead=\"dot\" color=\"forestgreen\"]", a);
			
			else			
				fprintf(graphFile, "   g$i -> y; [ arrowhead=\"tee\" color=\"firebrick1\"]", a);
		}
		on = k;
	}
	//end input 0
	
	//done with input 1. generating the rest is a systematic approach
	for (i=1; i < NUM_INPUTS; ++i)
	{
		
	}
	
	fprintf("}\n");
	free(x);
}


