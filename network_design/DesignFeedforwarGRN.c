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
void generateGraphFile(FILE * graphFile)
{
	int i,j,k,n,neg;
	double * input;
	
	if (!TARGET_FUNC || !GRID_SIZE) return;
	
	fprintf(graphFile, "digraph G {\n");
	
	//make the initial set of modules (input 1)
	
	n = GRID_SIZE[0]; //num dividers in input 1
	
	input = (double*)malloc(NUM_INPUTS * sizeof (double));
	
	neg = TARGET_FUNC(
	
	//done with input 1. generating the rest is a systematic approach
	for (i=0; i < NUM_INPUTS; ++i)
	{
		
	}
	
	fprintf("}\n");
}


