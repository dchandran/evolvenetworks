#include "blocks.h"
#include "functions.h"

static int MIN_SIZE = 1;  //minimum allowed size of a system (size = num. of blocks)
static int MAX_SIZE = 20;  //maximum allowed size of a system (size = num. of blocks)

static double PROB_PARAM_CHANGE = 0.6; //prob. of changing parameter (during mutation)
static double PROB_REWIRE = 0.3;  //prob. of rewiring vs. changing parameter (during mutation)
static double PROB_NEW_BLOCK = 0.05;  //prob. of adding new block (during mutation)
static double PROB_DEL_BLOCK = 0.05;  //prob. of removing a block (during mutation)
static int MUTATION_RATE = 1; //number of changes made during a single mutation event. Cannot be < 1

static double** FIXED_PARAM_VALUES = 0; //values for specific fixed parameters
static int** FIXED_PARAMS = 0; //indicates whether or not to mutate specific parameters
static int** FIXED_INPUTS = 0; //indicates whether or not to rewire specific inputs
static int** FIXED_OUTPUTS = 0; //indicates whether or not to rewire specific outputs
static int* ALLOWED_BLOCKS = 0; //indicates which block types to use
static int NO_SAME_INPUT_OUTPUT = 1; //whether the same species can be an input and output of same block

/*******************
* helper functions
********************/

static int stringsAreSame(const char * a, const char * b) // a == b
{
	int i=0;	
	if (!a || !b) return 0;
	
	while (a[i] && b[i] && a[i]==b[i]) ++i;
	
	return (!a[i] && !b[i] && i > 0);
}

/*************
* initialize
**************/

static void initialzeArrays() //initialize the above arrays
{
	int total, i, j;
	
	total = numBlockTypes();
	FIXED_PARAM_VALUES = (double**)malloc(total * sizeof(double*));
	FIXED_PARAMS = (int**)malloc(total * sizeof(int*));
	FIXED_INPUTS = (int*)malloc(total * sizeof(int));
	FIXED_OUTPUTS = (int*)malloc(total * sizeof(int));
	ALLOWED_BLOCKS = (int*)malloc(total * sizeof(int));
	
	for (i=0; i < total; ++i)
	{
		ALLOWED_BLOCKS[i] = 1; //all blocks are allowed by default
		
		FIXED_PARAM_VALUES[i] = (double*)malloc( BlockTypesTable[i].numParams * sizeof(double) );
		FIXED_PARAMS[i] = (int*)malloc( BlockTypesTable[i].numParams * sizeof(int) );
		FIXED_INPUTS[i] = (int*)malloc( BlockTypesTable[i].numParams * sizeof(int) );
		FIXED_OUTPUTS[i] = (int*)malloc( BlockTypesTable[i].numParams * sizeof(int) );
		
		//nothing fixed by default
		for (j=0; j < BlockTypesTable[i].numParams; ++j)
		{
			FIXED_PARAMS[i][j] = 0;
			FIXED_PARAM_VALUES[i][j] = 1.0;
		}
		
		for (j=0; j < BlockTypesTable[i].numInputs; ++j)
			FIXED_INPUTS[i][j] = 0;
	
		for (j=0; j < BlockTypesTable[i].numOutputs; ++j)
			FIXED_OUTPUTS[i][j] = 0;
	}
}

/**********************************************************************
* Basic functions for getting information about blocks and block types
***********************************************************************/

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

int isNullBlockType(BlockType type)
{
	return (type.name == 0 || type.numReactions < 1);
}

int isNullBlock(Block * block) 
{ 
	int i = block->type;
	return isNullBlockType(BlockTypesTable[i]); 
}

int getBlockTypeIndex(const char * name)
{
	int i;
	int sz = numBlockTypes();
	
	for (i=0; i < sz; ++i)
		if (stringsAreSame(name,BlockTypesTable[i].name))
			break;
	
	if (i >= sz) i = -1;
	return i;
}

const char* getBlockTypeName(int i)
{
	int sz = numBlockTypes();
	
	if (i < 0 || i >= sz) i = 0;
	return BlockTypesTable[i].name;
}

/************************
* set evolution settings
*************************/

void setMutationRate(int r)
{
	if (r > 0)
		MUTATION_RATE = r;
}

void allowSameInputAndOutput()
{
	NO_SAME_INPUT_OUTPUT = 0;
}

void disallowSameInputAndOutput()
{
	NO_SAME_INPUT_OUTPUT = 1;
}

void setMutateParameterProb(double d)
{
	if (d >= 0.0 && d <= 1.0) 
	{
		PROB_REWIRE = 1.0 - d;
		PROB_PARAM_CHANGE = d;
	}
}

void setRewiringProb(double d)
{
	setMutateParameterProb(1.0 - d);
}

void allowParameterChange()
{
	int total, i, j;	
	total = numBlockTypes();
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	for (i=0; i < total; ++i)	
		for (j=0; j < BlockTypesTable[i].numParams; ++j)
				FIXED_PARAMS[i][j] = 0;
				
	if (PROB_REWIRE == 1.0)
		PROB_REWIRE = 0.3;

	PROB_PARAM_CHANGE = 1 - PROB_REWIRE;
}

void disallowParameterChange()
{
	int total, i, j;	
	total = numBlockTypes();
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	for (i=0; i < total; ++i)	
		for (j=0; j < BlockTypesTable[i].numParams; ++j)
				FIXED_PARAMS[i][j] = 1;
	
	setMutateParameterProb(0.0);
}

void allowParameterChangeFor(const char* name)
{
	int j, i = getBlockTypeIndex(name);
	
	if (i < 0) return;
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	for (j=0; j < BlockTypesTable[i].numParams; ++j)
		FIXED_PARAMS[i][j] = 0;
}

void disallowParameterChangeFor(const char* name)
{
	int j, i = getBlockTypeIndex(name);
	
	if (i < 0) return;
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	for (j=0; j < BlockTypesTable[i].numParams; ++j)
		FIXED_PARAMS[i][j] = 1;
}

void fixParameter(const char* name, int j, double value)
{
	int i = getBlockTypeIndex(name);
	
	if (i < 0 || j < 0 || (j >= BlockTypesTable[i].numParams)) return;
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	FIXED_PARAMS[i][j] = 1;
	FIXED_PARAM_VALUES[i][j] = value;
}

void allowRewiring()
{
	int total, i, j;	
	total = numBlockTypes();
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	for (i=0; i < total; ++i)	
	{
		for (j=0; j < BlockTypesTable[i].numInputs; ++j)
				FIXED_INPUTS[i][j] = 0;
		
		for (j=0; j < BlockTypesTable[i].numOutputs; ++j)
				FIXED_OUTPUTS[i][j] = 0;
	}
	
	if (PROB_PARAM_CHANGE == 1.0)
		PROB_PARAM_CHANGE = 0.7;

	PROB_REWIRE = 1 - PROB_PARAM_CHANGE;
	
}

void disallowRewiring()
{
	int total, i, j;	
	total = numBlockTypes();
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	for (i=0; i < total; ++i)	
	{
		for (j=0; j < BlockTypesTable[i].numInputs; ++j)
				FIXED_INPUTS[i][j] = 1;
		
		for (j=0; j < BlockTypesTable[i].numOutputs; ++j)
				FIXED_OUTPUTS[i][j] = 1;
	}
	
	PROB_PARAM_CHANGE == 1.0;
	PROB_REWIRE = 0.0;
}

void allowRewiringFor(const char * name)
{
	int j, i = getBlockTypeIndex(name);
	
	if (i < 0) return;
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	for (j=0; j < BlockTypesTable[i].numInputs; ++j)
		FIXED_INPUTS[i][j] = 0;
		
	for (j=0; j < BlockTypesTable[i].numOutputs; ++j)
		FIXED_OUTPUTS[i][j] = 0;
}

void disallowRewiringFor(const char * name)
{
	int j, i = getBlockTypeIndex(name);
	
	if (i < 0) return;
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	for (j=0; j < BlockTypesTable[i].numInputs; ++j)
		FIXED_INPUTS[i][j] = 1;
		
	for (j=0; j < BlockTypesTable[i].numOutputs; ++j)
		FIXED_OUTPUTS[i][j] = 1;
}

void setInputFixed(const char * name, int j, int fixed)
{
	int i = getBlockTypeIndex(name);
	
	if (i < 0 || j < 0 || j > BlockTypesTable[i].numInputs) return;
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	FIXED_INPUTS[i][j] = fixed;
}

void setOutputFixed(const char * name, int input, int fixed)
{
	int i = getBlockTypeIndex(name);
	
	if (i < 0 || j < 0 || j > BlockTypesTable[i].numOutputs) return;
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	FIXED_OUTPUTS[i][j] = fixed;
}

int addBlockTypeByIndex(int i)
{
	if (i < 0 || i >= numBlockTypes()) return 0;	
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	ALLOWED_BLOCKS[i] = 1;
}

int addBlockType(const char * name)
{
	addBlockTypeByIndex( getBlockTypeIndex(name) );
}

int removeBlockTypeByIndex(int i)
{
	if (i < 0 || i >= numBlockTypes()) return 0;	
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	ALLOWED_BLOCKS[i] = 0;
}

int removeBlockType(const char * name)
{
	removeBlockTypeByIndex(getBlockTypeIndex(name));
}

void setAddBlockProb(double d)
{
	if (d >= 0.0 && d <= 1.0)
		PROB_NEW_BLOCK = d;
}

void setRemoveBlockProb(double d)
{
	if (d >= 0.0 && d <= 1.0)
		PROB_DEL_BLOCK = d;
}

void setSizeRange(int min, int max)
{
	if (min > 0 && max > 0)
	{
		MIN_SIZE = min;
		MAX_SIZE = max;
	}
}

/*********************************************
* for use by GA (see ga.h)
**********************************************/

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

static Block * copyBlock(Block * block)
{
	int i;
	int n;
	Block * block2 = (Block*)malloc(sizeof(Block));

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
	free(s->blocks);
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
	
	s2->blocks = (Block**)malloc(numBlocks * sizeof(Block*));
	
	for (i=0; i < numBlocks; ++i)
		s2->blocks[i] = copyBlock(s->blocks[i]);

	s2->numBlocks = s->numBlocks;
	s2->numSpecies = s->numSpecies;
	
	return s2;
}

/*****************************************************************************************
* crossover two systems by taking random blocks from each and rewiring the severed links
******************************************************************************************/

static int totalInternalSpecies(SystemOfBlocks * s)
{
	int i = 0, sz = 0;
	
	for (i=0; i < s->size; ++i)
		sz += BlockTypesTable[ s->blocks[i]->type ].numInternals;
	
	return sz;
}

static void pruneSystemOfBlocks(SystemOfBlocks * s) //find unconnected inputs/outputs
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

static SystemOfBlocks * randomSubsystem(SystemOfBlocks * s, double prob) //get random subset of blocks
{
	int i, int i2;
	int numBlocks = s->numBlocks, numBlocks2;
	SystemOfBlocks * s2 = (SystemOfBlocks*)malloc(sizeof(SystemOfBlocks));
	
	numBlocks2 = (int)(prob * (double)numBlocks);
	
	s2->blocks = (Block**)malloc(numBlocks2 * sizeof(Block*));
	
	for (i=0, i2=0; i < numBlocks && i2 < numBlocks2; ++i)
	{
		if (mtrand() < prob)
		{
			s2->blocks[i2] = copyBlock(s->blocks[i]);
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

static GAindividual * GAcrossoverBlocks(GAindividual X, GAindividual Y) //place two subsets together
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
	
	s3->blocks = (Block**)malloc(s3->numBlocks * sizeof(Block));
	
}

/*******************
* mutation
********************/

static GAindividual GAmutateBlocksH(GAindividual X)  //rewire or change parameter
{
	SystemOfBlocks * S = (SystemOfBlocks*)X;
	int i,j,k,k0,k1,m0,m1,n,a0,a1;
	
	int species = s->numSpecies,
		blocks = s->numBlocks,
		type = s->type;
	
	Block * oldBlocks;
	
	if (mtrand() < PROB_PARAM_CHANGE) //mutate parameter
	{	
		k0 = (int)(mtrand() * numBlocks);
		k1 = k0;
		
		while (1)
		{	
			n = numParams(S->blocks[k1]);
			m0 = (int)(mtrand() * n);
			
			m1 = m0;
			while (FIXED_PARAMS[type][k3])
			{
				++m1;
				if (m1 >= n) m1 = 0; //cycle
				if (m1 == m0) break;
			}
			
			if (!FIXED_PARAMS[type][m1] || m1!=m0)
			{
				S->blocks[k1]->params[m1] *= 2.0 * mtrand();  //mutate (k0,k3)
				return (GAindividual)S;
			}
			
			++k1;
			if (k1 >= blocks) k1 = 0;
			if (k0==k1) break;
		}
	}
	
	//add or remove block?
	if (mtrand() < (PROB_NEW_BLOCK + PROB_DEL_BLOCK))
	{
		if (mtrand() < (PROB_NEW_BLOCK/(PROB_NEW_BLOCK + PROB_DEL_BLOCK)))
		{
			oldBlocks = S->blocks;
			S->blocks = (Block**)malloc(
		}
		else
		{
		}
	}
	
	//rewire
	
	//first block and input
		
		a0 = a1 = -1;
		k0 = (int)(mtrand() * numBlocks);
		k1 = k0;
		
		while (1)
		{	
			n = numInputs(S->blocks[k1]);
			m0 = (int)(mtrand() * n);
			
			m1 = m0;
			while (FIXED_INPUTS[type][k3])
			{
				++m1;
				if (m1 >= n) m1 = 0; //cycle
				if (m1 == m0) break;
			}
			
			if (!FIXED_PARAMS[type][m1] || m1!=m0)
			{
				S->blocks[k1]->params[m1] *= 2.0 * mtrand();  //mutate (k0,k3)
				return (GAindividual)S;
			}
			
			++k1;
			if (k1 >= blocks) k1 = 0;
			if (k0==k1) break;
		}
	
	//second block and output
	if (a0 >= 0 && a1 >= 0)
	{
		
	}
	
	return (GAindividual)S;
}

static GAindividual GAmutateBlocks(GAindividual X)  //rewire or change parameter n times
{
	int i;
	
	if ((PROB_PARAM_CHANGE + PROB_REWIRE) > 0.0)
		for (i=0; i < MUTATION_RATE; ++i)
			GAmutateBlocksH(X);

	return X;
}

/**************************
* main evolution function
**************************

GApopulation evolveNetworks(int initialPopulationSize, int finalPopulationSize, GAfitnessFunc fitness, GAcallback callback)
{
	return 0;
}
