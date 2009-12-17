#include <stdlib.h>
#include <stdio.h>

#define bool int

#include "antimony_api.h"

void convertToBlock(const char * name, const char * antimonyFile, const char * file1, const char * file2)
{
	double ** N;
	int i,j,numInputs;
	unsigned long m,n,numSpecies, numParams;
	char ** rates, ** species, ** params, ** paramValues, ** speciesValues;
	char * mainModule = "__main";
	FILE * fp1, * fp2;

	if (loadFile(antimonyFile) == -1)
	{
		printf("%s\n",getLastError());
		return;
	}

	fp1 = fopen(file1,"a");
	fp2 = fopen(file2,"a");

	N = getStoichiometryMatrix (mainModule);
	m = getStoichiometryMatrixNumRows (mainModule);
	n = getStoichiometryMatrixNumColumns (mainModule);
	rates = getReactionRates (mainModule);

	numSpecies = getNumSymbolsOfType (mainModule, 1);
	species = getSymbolNamesOfType (mainModule, 1);
	speciesValues = getSymbolEquationsOfType (mainModule, 1);

	numParams = getNumSymbolsOfType (mainModule, 16);
	params = getSymbolNamesOfType (mainModule, 16);
	paramValues = getSymbolEquationsOfType (mainModule, 16);

	fprintf(fp1,"void %s_init(Block * block) {\n",name);

	for (i=0; i < numParams; ++i)
		fprintf(fp1,"\tblock->params[%i] = %s;\n",i,paramValues[i]);

    for (i=0; i < numSpecies; ++i)
		fprintf(fp1,"\tblock->initVals[%i] = %s;\n",i,speciesValues[i]);

	fprintf(fp1,"}\n\nvoid %s_rates(double t, double * y, double * r, Block * block) {\n",name);

	for (i=0; i < numParams; ++i)
		fprintf(fp1,"\tdouble %s = block->params[%i];\n",params[i],i);

	for (i=0; i < numSpecies; ++i)
		fprintf(fp1,"\tdouble %s = y[%i];\n",species[i],i);

	for (i=0; i < n; ++i)
		fprintf(fp1,"\tr[%i] = %s;\n",i,rates[i]);

	fprintf(fp1,"}\n\nvoid %s_stoic(Matrix * m, Block * block, int i) {\n",name);

	for (i=0; i < m; ++i)
	{
		for (j=0; j < n; ++j)
			if (N[i][j] != 0)
				fprintf(fp1,"\tvalueAt(*m, block->externals[%i], i+%i) += %lf;\n",i,j,N[i][j]);
		if (i < (m-1))
			fprintf(fp1,"\n");
	}

	fprintf(fp1,"}\n\n");

	fprintf(fp1,"double %s_paramLowerBound[] = {",name);
	for (i=0; i < numParams; ++i)
		if (i < (numParams-1))
			fprintf(fp1,"0.0,",name);
		else
			fprintf(fp1,"0.0};\n\n",name);

	fprintf(fp1,"double %s_paramUpperBound[] = {",name);
	for (i=0; i < numParams; ++i)
		if (i < (numParams-1))
			fprintf(fp1,"1.0E30,",name);
		else
			fprintf(fp1,"1.0E30};\n\n",name);


	numInputs = 0;
	fprintf(fp2,"    {\"%s\",&%s_stoic,&%s_rates,&%s_init,%i,%i,%i,%i,1,%s_paramLowerBound,%s_paramUpperBound},\n",name,name,name,name,n,numSpecies,numInputs,numParams,name,name);

	fclose(fp1);
	fclose(fp2);

	freeAll();
}

int main(int args, char ** argv)
{
	convertToBlock(argv[1],argv[2],argv[3],argv[4]);
}
