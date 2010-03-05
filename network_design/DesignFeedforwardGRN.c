#include <math.h>
#include <stdlib.h>
#include <stdio.h>
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
	//NUM_INPUTS = n;
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
	int i, i0,i1,i2,i3,j,k,n0,n1,n2,a,b,c0,c1,on;
	double * x;
	
	if (!TARGET_FUNC || !GRID_SIZE) return;
	
	a = b = c = 0;
	x = (double*)malloc(NUM_INPUTS * sizeof(double)); //inputs
	
	fprintf(graphFile, "digraph G {\n\n");
	n0 = GRID_SIZE[0]; //num dividers in input 1
	for (i=0; i < NUM_INPUTS; ++i)	x[i] = 0.0;  //initialize

	for (i=1; i < NUM_INPUTS; ++i)
	{		
	
		n1 = GRID_SIZE[i];
		for (i1=0; i1 < n1; ++i1) //input 1 gradually increasing
		{
			x[i] = i; //set input 1
			
			if (c > 0 && b > 0)
			{
				fprintf(graphFile, "   h%i -> g%i [ arrowhead=\"dot\" , color=\"forestgreen\" ]\n", b, c);  //turn off previous "on button"
			}
			
			//input 0 at each level of input 1
			x[0] = 0;
			on = TARGET_FUNC(x);  //is the output on when input 0 == 0?
			fprintf(graphFile, "   g%i -> Y [ arrowhead=\"dot\" , color=\"forestgreen\" ]\n", ++a); //the "on button" for this module
			c1 = 0;
			c0 = a; //save on button
			if (!on)
			{
				fprintf(graphFile, "   g%i -> Y [ arrowhead=\"tee\" , color=\"firebrick1\" ]\n", ++a); //repressor
				c1 = a;
			}
	
			for (i0=1; i0 < n0; ++i0)
			{
				x[0] = i0;
				k = TARGET_FUNC(x);
				if (k != on)
				{
					fprintf(graphFile, "   x0 -> g%i [ arrowhead=\"dot\" , color=\"forestgreen\" ]\n", ++a);
					fprintf(graphFile, "   g%i -> g%i [ arrowhead=\"tee\" , color=\"firebrick1\" ]\n", a, a-1);
				}
				on = k;
			}
			//end input 0
			
			if (i1 > 0)
			{
				fprintf(graphFile, "   x1 -> h%i [ arrowhead=\"dot\" , color=\"forestgreen\" ]\n", ++b);
				fprintf(graphFile, "   h%i -> g%i [ arrowhead=\"tee\" , color=\"firebrick1\" ]\n", b, ++a); //lift the repressor
				fprintf(graphFile, "   g%i -> g%i [ arrowhead=\"tee\" , color=\"firebrick1\" ]\n", a, c0); //repress "on button"
				fprintf(graphFile, "   g%i -> g%i [ arrowhead=\"tee\" , color=\"firebrick1\" ]\n", a, c0); //repress "on button"
			}
		}
	}
	
	fprintf(graphFile,"\n}\n");
	free(x);
}


