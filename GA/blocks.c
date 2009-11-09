#include "blocks.h"

static int MIN_SIZE = 3;
static int MAX_SIZE = 20;
static int PROB_PARAM_CHANGE = 1;
static int PROB_REWIRE = 1;

static void freeBlock(Block * block)
{
	if (!block) return;
	
	if (block->internals)
		free(block->internals);
		
	if (block->inputs)
		free(block->inputs);
	
	if (block->outputs)
		free(block-outputs);
	
	if (block->params)
		free(block->params);
	
	free(block);
}

static void copyBlock(Block * block, Block * block2)
{
	int i;
	int n;

	block2->type = block->type;

	n = BlockTypesTable[ block->type ].numInternals;
	block2->internals = (int*)malloc(n * sizeof(int));
	for (i=0; i < n; ++i)
		block2->internals[i] = block->internals[i];

	n = BlockTypesTable[ block->type ].numInputs;
	block2->inputs = (int*)malloc(n * sizeof(int));
	for (i=0; i < n; ++i)
		block2->inputs[i] = block->inputs[i];

	n = BlockTypesTable[ block->type ].numOutputs;
	block2->outputs = (int*)malloc(n * sizeof(int));
	for (i=0; i < n; ++i)
		block2->outputs[i] = block->outputs[i];

	n = BlockTypesTable[ block->type ].numParameters;
	block2->params = (double*)malloc(n * sizeof(double));
	for (i=0; i < n; ++i)
		block2->params[i] = block->params[i];

	return block2;
}

static System * cloneSystem(System * s)
{
	int i;
	int numBlocks = s->numBlocks;
	System * s2 = (System*)malloc(sizeof(System));
	
	s2->blocks = (Block*)malloc(numBlocks * sizeof(Block*));
	
	for (i=0; i < numBlocks; ++i)
		copyBlock(s->blocks[i],s2->blocks[i]);

	s2->numBlocks = s->numBlocks;
	s2->numSpecies = s->numSpecies;
	
	return s2;
}

static int totalInternalSpecies(System * s)
{
	int i = 0, sz = 0;
	
	for (i=0; i < s->size; ++i)
		sz += BlockTypesTable[ s->blocks[i]->type ].numInternals;
	
	return sz;
}

static void pruneSystem(System * s)
{
	int i,j,k,l,m,n,b0,b1,b2,total;
	
	total = s->numSpecies;
	
	for (i=0; i < s->numSpecies; ++i)
	{
		b0 = b1 = b2 = 0;
		for (j=0; j < s->numBlocks; ++j)
		{
			n = BlockTypesTable[ s->blocks[j]->type ].numInternals; 
			for (k=0; k < n; ++k)
				if (s->blocks[j]->internals[k] == i)
				{
					b0 = 1;
					break;
				}
			if (b0) break;
			
			n = BlockTypesTable[ s->blocks[j]->type ].numInputs;
			for (k=0; k < n; ++k)
				if (s->blocks[j]->inputs[k] == i)
				{
					b1 = 1;
					break;
				}
			n = BlockTypesTable[ s->blocks[j]->type ].numOutputs;
			for (k=0; k < n; ++k)
				if (s->blocks[j]->outputs[k] == i)
				{
					b2 = 1;
					break;
				}
		}
		
		if (!b0 && (!b1 || !b2))
		{
			--total;
			n = BlockTypesTable[ s->blocks[j]->type ].numInternals; 
			for (k=0; k < n; ++k)
			{
				if (s->blocks[j]->internals[k] == i)
					s->blocks[j]->internals[k] = -1;
				
				if (s->blocks[j]->internals[k] > i)
					s->blocks[j]->internals[k] -= 1;
			}
			if (b0) break;
			
			n = BlockTypesTable[ s->blocks[j]->type ].numInputs;
			for (k=0; k < n; ++k)
			{
				if (s->blocks[j]->inputs[k] == i)
					s->blocks[j]->inputs[k] = -1;
					
				if (s->blocks[j]->inputs[k] > i)
					s->blocks[j]->inputs[k] -= 1;
			}
			
			n = BlockTypesTable[ s->blocks[j]->type ].numOutputs;
			for (k=0; k < n; ++k)
			{
				if (s->blocks[j]->outputs[k] == i)
					s->blocks[j]->outputs[k] = -1;
				
				if (s->blocks[j]->outputs[k] > i)
					s->blocks[j]->outputs[k] -= 1;
			}
		}
	}
	s->numSpecies = total;
}

static System * randomSubsystem(System * s, double prob)
{
	int i, int i2;
	int numBlocks = s->numBlocks, numBlocks2;
	System * s2 = (System*)malloc(sizeof(System));
	
	numBlocks2 = (int)(prob * (double)numBlocks);
	
	s2->blocks = (Block*)malloc(numBlocks2 * sizeof(Block*));
	
	for (i=0, i2=0; i < numBlocks && i2 < numBlocks2; ++i)
	{
		if (mtrand() < prob)
		{
			copyBlock(s->blocks[i],s2->blocks[i2]);
			++i2;
		}
		
		if (i >= numBlocks)
			i = 0;
	}

	s2->numBlocks = numBlocks2;
	s2->numSpecies = s->numSpecies;
	
	pruneSystem(s2);
	return s2;
}

static GAindividual * GAcrossoverBlocks(GAindividual * X, GAindividual * Y)
{
	System * s1 = randomSubsystem((System*)X, 0.5);
	System * s2 = randomSubsystem((System*)Y, 0.5);
	System * s3;
	int sz1 = s1->numSpecies;
	int i,j,n;
	
	for (i=0; i < s2->numBlocks; ++i)
	{
		n = BlockTypesTable[ s2->blocks[i]->type ].numInternals;
		for (j=0; j < n; ++j)
			s2->internals[j] += sz1;
		
		n = BlockTypesTable[ s2->blocks[i]->type ].numInputs;
		for (j=0; j < n; ++j)
			s2->inputs[j] += sz1;
		
		n = BlockTypesTable[ s2->blocks[i]->type ].numOutputs;
		for (j=0; j < n; ++j)
			s2->outputs[j] += sz1;
	}
	
	s3 = (System*)malloc(sizeof(System));
	
	s3->numBlocks = s1->numBlocks + s2->numBlocks;
	s3->numSpecies = s1->numSpecies + s2->numSpecies;
	
	s3->blocks = (Block*)malloc(s3->numBlocks * sizeof(Block));
	
}

static GAindividual * GAmutateBlocks(GAindividual * X, GAindividual * Y)
{
	
}

void setMutateParameterProb(double d)
{
	PROB_PARAM_CHANGE = d;
}

void setMutateStructureProb(double d)
{
	PROB_REWIRE = d;
}

int isNullBlockType(BlockType type)
{
	return (type.name == 0 || type.numReactions < 1);
}

int isNullBlock(Block * block) 
{ 
	int i = block->type;
	return isNullBlockType(BlockTypesTable[i]); 
}

static int numBlockTypes()
{
	int i = 0;
	while (!isNullBlockType(BlockTypesTable[i]))
		++i;
	return i;
}

void addBlockType(int i)
{
	int sz = numBlockTypes();
	if (i < 0 || i >= sz) return;
	
	
}

GApopulation evolveNetworks(int initialPopulationSize, int finalPopulationSize, GAfitnessFunc fitness, GAcallback callback)
{
}
