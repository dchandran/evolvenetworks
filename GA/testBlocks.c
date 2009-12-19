#include "blocks.h"

int main(int args, char** argv)
{
	int i,j,sz;
	double * y = 0;
	Matrix M;
	System S;
	GApopulation P;

	initMTrand();
/*
	S.numBlocks = 2;
	S.blocks = (Block**)malloc(2*sizeof(Block*));
	S.blocks[0] = (Block*)malloc(sizeof(Block));
	S.blocks[1] = (Block*)malloc(sizeof(Block));

	S.blocks[0]->type = getBlockTypeIndex("enzyme_catalysis");
	S.blocks[0]->params = (double*)malloc(numParams(S.blocks[0])*sizeof(double));
	S.blocks[0]->internals = (int*)malloc(numInternals(S.blocks[0])*sizeof(int));
	S.blocks[0]->externals = (int*)malloc(numExternals(S.blocks[0])*sizeof(int));
	S.blocks[0]->initVals = (double*)malloc(numSpecies(S.blocks[0])*sizeof(double));

	S.blocks[1]->type = getBlockTypeIndex("uniuni");
	S.blocks[1]->params = (double*)malloc(numParams(S.blocks[1])*sizeof(double));
	S.blocks[1]->internals = (int*)malloc(numInternals(S.blocks[1])*sizeof(int));
	S.blocks[1]->externals = (int*)malloc(numExternals(S.blocks[1])*sizeof(int));
	S.blocks[1]->initVals = (double*)malloc(numSpecies(S.blocks[1])*sizeof(double));

	S.numSpecies = numExternals(S.blocks[0]);
	initializeSystem(&S);

	for (i=0; i < S.numSpecies; ++i)
		S.blocks[0]->externals[i] = i;

	S.blocks[1]->externals[0] = 0;
	S.blocks[1]->externals[1] = S.numSpecies-1;

	//y = simulateODE(&S,y0,10.0,0.1);
	sz = 100;
	y = simulateStochastic(&S,100.0,&sz);

	if (y)
	{
		for (i=0; i < sz; ++i)
		{
			for (j=0; j < 6; ++j)
				printf("%lf\t",getValue(y,6,i,j));
			printf("\n");
		}
		free(y);
	}
*/
	
	return 0;
}

