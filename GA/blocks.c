#include "blocks.h"
#include "functions.h"

static int MIN_SIZE = 3;  //minimum allowed size of a system (size = num. of blocks)
static int MAX_SIZE = 20;  //maximum allowed size of a system (size = num. of blocks)
static int PROB_PARAM_CHANGE = 1; //prob. of changing parameter (during mutation)
static int PROB_REWIRE = 1;  //prob. of rewiring vs. changing parameter (during mutation)
static int** FIXED_PARAMS = 0; //indicates whether or not to mutate specific parameters
static int* ALLOWED_BLOCKS = 0; //indicates which block types to use

int numBlockTypes()
{
	int i = 0;
	while (BlockTypesTable && !isNullBlockType(BlockTypesTable[i]))
		++i;
	return i;
}

int numInputs(Block * block)
{
	return BlockTypesTable[ block->type ].numInputs;
}

int numOutputs(Block * block)
{
	return BlockTypesTable[ block->type ].numOutputs;
}

int numParams(Block * block)
{
	return BlockTypesTable[ block->type ].numParams;
}

int numInternals(Block * block)
{
	return BlockTypesTable[ block->type ].numInternals;
}

int numReactions(Block * block)
{
	return BlockTypesTable[ block->type ].numReactions;
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

static void initialzeArrays()
{
	int total, i, j;
	
	total = numBlockTypes();
	FIXED_PARAMS = (int**)malloc(total * sizeof(int*));
	ALLOWED_BLOCKS = (int*)malloc(total * sizeof(int));
	
	for (i=0; i < total; ++i)
	{
		ALLOWED_BLOCKS[i] = 1; //all blocks are allowed by default
		
		FIXED_PARAMS[i] = (int*)malloc( BlockTypesTable[i].numParams * sizeof(int) );
		for (j=0; j < BlockTypesTable[i].numParams; ++j)
		{
			FIXED_PARAMS[i][j] = 0; //no parameters are fixed by default
		}
	}
}

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

static void freeSystemOfBlocks(GAindividual X)
{
	SystemOfBlocks * s = (SystemOfBlocks*)X;
	int i;
	int numBlocks;

	numBlocks = s->numBlocks;
	for (i=0; i < numBlocks; ++i)
		freeBlock(s->blocks[i]);
	free(s>blocks);
	free(s);
}

static SystemOfBlocks * cloneSystemOfBlocks(GAindividual X)
{
	SystemOfBlocks * s = (SystemOfBlocks*)X;
	int i;
	int numBlocks;
	SystemOfBlocks * s2;

	numBlocks = s->numBlocks;
	s2 = (SystemOfBlocks*)malloc(sizeof(SystemOfBlocks));
	
	s2->blocks = (Block*)malloc(numBlocks * sizeof(Block*));
	
	for (i=0; i < numBlocks; ++i)
		copyBlock(s->blocks[i],s2->blocks[i]);

	s2->numBlocks = s->numBlocks;
	s2->numSpecies = s->numSpecies;
	
	return s2;
}

static int totalInternalSpecies(SystemOfBlocks * s)
{
	int i = 0, sz = 0;
	
	for (i=0; i < s->size; ++i)
		sz += BlockTypesTable[ s->blocks[i]->type ].numInternals;
	
	return sz;
}

static void pruneSystemOfBlocks(SystemOfBlocks * s)
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

static SystemOfBlocks * randomSubsystem(SystemOfBlocks * s, double prob)
{
	int i, int i2;
	int numBlocks = s->numBlocks, numBlocks2;
	SystemOfBlocks * s2 = (SystemOfBlocks*)malloc(sizeof(SystemOfBlocks));
	
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
	
	pruneSystemOfBlocks(s2);
	return s2;
}

static GAindividual * GAcrossoverBlocks(GAindividual X, GAindividual Y)
{
	SystemOfBlocks * s1 = randomSubsystem((SystemOfBlocks*)X, 0.5);
	SystemOfBlocks * s2 = randomSubsystem((SystemOfBlocks*)Y, 0.5);
	SystemOfBlocks * s3;
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
	
	s3 = (SystemOfBlocks*)malloc(sizeof(SystemOfBlocks));
	
	s3->numBlocks = s1->numBlocks + s2->numBlocks;
	s3->numSpecies = s1->numSpecies + s2->numSpecies;
	
	s3->blocks = (Block*)malloc(s3->numBlocks * sizeof(Block));
	
}

static GAindividual GAmutateBlocks(GAindividual X)
{
	SystemOfBlocks * s = (SystemOfBlocks*)X;
	int i,j,n;
	
	i = (int)(mtrand() * s->numBlocks);
	
	n = BlockTypesTable[ s->blocks[i]->type ].numParams;
	j = (int)(mtrand() * n)
	
	s->blocks[i]->params[j] *= (2.0 * mtrand());
	
	return (GAindividual)s;
	
	//mutate a parameter
	BlockTypesTable[ block->type ].numInputs;
	int * inputs;
	int * outputs;
	int * internals;
	double * params;
}

static int stringsAreSame(const char * a, const char * b)
{
	int i=0;	
	if (!a || !b) return 0;
	
	while (a[i] && b[i] && a[i]==b[i]) ++i;
	
	return (!a[i] && !b[i] && i > 0);
}

int addBlockTypeByIndex(int i)
{
	if (i < 0 || i >= numBlockTypes()) return 0;	
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	ALLOWED_BLOCKS[i] = 1;
}

int addBlockType(const char * name)
{
	int i;
	int sz = numBlockTypes();
	
	for (i=0; i < sz; ++i)
		if (stringsAreSame(name,BlockTypesTable[i].name))
			break;

	addBlockTypeByIndex(i);
}

int removeBlockTypeByIndex(int i)
{
	if (i < 0 || i >= numBlockTypes()) return 0;	
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	ALLOWED_BLOCKS[i] = 0;
}

int removeBlockType(const char * name)
{
	int i;
	int sz = numBlockTypes();
	
	for (i=0; i < sz; ++i)
		if (stringsAreSame(name,BlockTypesTable[i].name))
			break;
	
	removeBlockTypeByIndex(i);
}

GApopulation evolveNetworks(int initialPopulationSize, int finalPopulationSize, GAfitnessFunc fitness, GAcallback callback)
{
	return 0;
}
