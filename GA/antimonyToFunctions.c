#include <stdlib.h>
#include <stdio.h>

#define bool int

#include "antimony_api.h"

static int stringInList(const char * a, char ** list, int len) // a == b
{
	int i=0,j=0,k=0;
	char * b;
	if (!a || !list || len < 1) return 0;

	for (j=0; j < len; ++j)
	{
		b = list[j];
		i = 0;
		while (a[i] && b[i] && a[i]==b[i]) ++i;
		if (!a[i] && !b[i] && i > 0) return 1;
	}
	return 0;
}

void convertToBlock(const char * antimonyFile, const char * file1, const char * file2)
{
	double ** N;
	int i,j,k,l,p;
	unsigned long m,n,numSpecies, numParams, numModules, numInterfaces;
	char ** rates, ** speciesNames, ** params, ** paramValues, ** speciesValues, ** interfaceNames;
	char * module, * name;
	FILE * fp1, * fp2;

	if (loadFile(antimonyFile) == -1)
	{
		printf("%s\n",getLastError());
		return;
	}

	fp1 = fopen(file1,"w");
	fp2 = fopen(file2,"w");

	fprintf(fp2,"#ifndef GA_BLOCKS_TABLE_H\n#define GA_BLOCKS_TABLE_H\n#include \"blockFunctions.h\"\n\nBlockType BlockTypesTable[] = \n{\n");

	numModules = getNumModules();

	for (l=0; l < numModules; ++l)
	{
		module = getNthModuleName(l);
		name = module;

		numInterfaces = getNumSymbolsInInterfaceOf(module);
		interfaceNames = getSymbolNamesInInterfaceOf(module);

		N = getStoichiometryMatrix (module);
		m = getStoichiometryMatrixNumRows (module);
		n = getStoichiometryMatrixNumColumns (module);
		rates = getReactionRates (module);

		numSpecies = getNumSymbolsOfType (module, 1);
		speciesNames = getSymbolNamesOfType (module, 1);
		speciesValues = getSymbolEquationsOfType (module, 1);

		numParams = getNumSymbolsOfType (module, 16);
		params = getSymbolNamesOfType (module, 16);
		paramValues = getSymbolEquationsOfType (module, 16);

		if (numParams < 1 || m < 1 || n < 1 || numInterfaces < 1) continue;

		printf("module: %s\n",name);

		fprintf(fp1,"void %s_init(Block * block) {\n",name);

		for (i=0; i < numParams; ++i)
			fprintf(fp1,"\tblock->params[%i] = %s;\n",i,paramValues[i]);

		for (i=0,j=0,k=0; i < numSpecies; ++i)
			if (stringInList(speciesNames[i],interfaceNames,numInterfaces))
			{
				fprintf(fp1,"\tblock->initValsExternals[%i] = %s;\n",j,speciesValues[i]);
				++j;
			}
			else
			{
				fprintf(fp1,"\tblock->initValsInternals[%i] = %s;\n",k,speciesValues[i]);
				++k;
			}

		fprintf(fp1,"}\n\nvoid %s_rates(double t, double * y, double * r, Block * block) {\n",name);

		for (i=0; i < numParams; ++i)
			fprintf(fp1,"\tdouble %s = block->params[%i];\n",params[i],i);

		for (i=0,j=0,k=0; i < numSpecies; ++i)
			if (stringInList(speciesNames[i],interfaceNames,numInterfaces))
			{
				fprintf(fp1,"\tdouble %s = y[ block->externals[%i] ];\n",speciesNames[i],j);
				++j;
			}
			else
			{
				fprintf(fp1,"\tdouble %s = y[ block->internals[%i] ];\n",speciesNames[i],k);
				++k;
			}

		for (i=0; i < n; ++i)
			fprintf(fp1,"\tr[%i] = %s;\n",i,rates[i]);

		fprintf(fp1,"}\n\nvoid %s_stoic(Matrix * m, Block * block, int i) {\n",name);

		for (i=0,p=0,k=0; i < m; ++i)
		{
			if (stringInList(speciesNames[i],interfaceNames,numInterfaces))
			{
				for (j=0; j < n; ++j)
					if (N[i][j] != 0)
						fprintf(fp1,"\tvalueAt(*m, block->externals[%i], i+%i) += %lf;\n",p,j,N[i][j]);
				++p;
			}
			else
			{
				for (j=0; j < n; ++j)
					if (N[i][j] != 0)
						fprintf(fp1,"\tvalueAt(*m, block->internals[%i], i+%i) += %lf;\n",k,j,N[i][j]);
				++k;
			}

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

		fprintf(fp2,"    {\"%s\",&%s_stoic,&%s_rates,&%s_init,%i,%i,%i,%i,1,%s_paramLowerBound,%s_paramUpperBound},\n",name,name,name,name,n,numInterfaces,numSpecies-numInterfaces,numParams,name,name);
	}

	fprintf(fp2,"\n    {0,0,0,0,0,0,0,0,0,0}\n};\n\n#endif\n");

	fclose(fp1);
	fclose(fp2);

	freeAll();
}

int main(int args, char ** argv)
{
	convertToBlock(argv[1],argv[2],argv[3]);
}
