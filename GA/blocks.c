#include <stdlib.h>
#include <stdio.h>
#include "blocks.h"
#include "functions.h"

static int MIN_SIZE = 1;  //minimum allowed size of a system (size = num. of blocks)
static int MAX_SIZE = 20;  //maximum allowed size of a system (size = num. of blocks)

static double PROB_PARAM_CHANGE = 0.0; //prob. of changing parameter (during mutation)
static double PROB_REWIRE = 0.3;  //prob. of rewiring vs. changing parameter (during mutation)
static double PROB_NEW_BLOCK = 0.0;  //prob. of adding new block (during mutation)
static double PROB_DEL_BLOCK = 0.0;  //prob. of removing a block (during mutation)
static int MUTATION_RATE = 1; //number of changes made during a single mutation event. Cannot be < 1

static double** FIXED_PARAM_VALUES = 0; //values for specific fixed parameters
static int** FIXED_PARAMS = 0; //indicates whether or not to mutate specific parameters
static int** FIXED_INPUTS = 0; //indicates whether or not to rewire specific inputs
static int** FIXED_OUTPUTS = 0; //indicates whether or not to rewire specific outputs
static int* ALLOWED_BLOCKS = 0; //indicates which block types to use
static int NO_SAME_INPUT_OUTPUT = 0; //whether the same species can be an input and output of same block

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
	FIXED_INPUTS = (int**)malloc(total * sizeof(int*));
	FIXED_OUTPUTS = (int**)malloc(total * sizeof(int*));
	ALLOWED_BLOCKS = (int*)malloc(total * sizeof(int));
	
	for (i=0; i < total; ++i)
	{
		ALLOWED_BLOCKS[i] = 1; //all blocks are allowed by default
		
		FIXED_PARAM_VALUES[i] = (double*)malloc( BlockTypesTable[i].numParams * sizeof(double) );
		FIXED_PARAMS[i] = (int*)malloc( BlockTypesTable[i].numParams * sizeof(int) );
		FIXED_INPUTS[i] = (int*)malloc( BlockTypesTable[i].numInputs * sizeof(int) );
		FIXED_OUTPUTS[i] = (int*)malloc( BlockTypesTable[i].numOutputs * sizeof(int) );
		
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
	
	PROB_PARAM_CHANGE = 1.0;
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

void setOutputFixed(const char * name, int j, int fixed)
{
	int i = getBlockTypeIndex(name);
	
	if (i < 0 || j < 0 || j > BlockTypesTable[i].numOutputs) return;
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	
	FIXED_OUTPUTS[i][j] = fixed;
}

void addBlockTypeByIndex(int i)
{
	if (i < 0 || i >= numBlockTypes()) return;	
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	ALLOWED_BLOCKS[i] = 1;
}

void addBlockType(const char * name)
{
	addBlockTypeByIndex( getBlockTypeIndex(name) );
}

void removeBlockTypeByIndex(int i)
{
	if (i < 0 || i >= numBlockTypes()) return;	
	
	if (!ALLOWED_BLOCKS) initialzeArrays();
	ALLOWED_BLOCKS[i] = 0;
}

void removeBlockType(const char * name)
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

int systemSize(System* s)
{
	int i,n,total=0;
	
	n = s->numBlocks;
	for (i=0; i < n; ++i)
		total += numReactions(s->blocks[i]);

	return total;
}

double paramUpperBound(Block * block, int k)
{
	if (!block || block->type > numBlockTypes()) return 1.0;

	if (BlockTypesTable[block->type].paramsUpperBound && 
		k >= 0 &&
		k < BlockTypesTable[block->type].numParams)
		return BlockTypesTable[block->type].paramsUpperBound[k];
	return 1.0;
}

double paramLowerBound(Block * block, int k)
{
	if (!block || block->type > numBlockTypes()) return 0.0;

	if (BlockTypesTable[block->type].paramsLowerBound && 
		k >= 0 &&
		k < BlockTypesTable[block->type].numParams)
		return BlockTypesTable[block->type].paramsLowerBound[k];
	return 0.0;
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
		free(block->outputs);
	
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

	n = numInternals(block);
	block2->internals = (int*)malloc(n * sizeof(int));
	for (i=0; i < n; ++i)
		block2->internals[i] = block->internals[i];

	n = numInputs(block);
	block2->inputs = (int*)malloc(n * sizeof(int));
	for (i=0; i < n; ++i)
		block2->inputs[i] = block->inputs[i];

	n = numOutputs(block);
	block2->outputs = (int*)malloc(n * sizeof(int));
	for (i=0; i < n; ++i)
		block2->outputs[i] = block->outputs[i];

	n = numParams(block);
	block2->params = (double*)malloc(n * sizeof(double));
	for (i=0; i < n; ++i)
		block2->params[i] = block->params[i];

	return block2;
}

static void freeSystem(GAindividual X)
{
	System * s = (System*)X;
	int i;
	int numBlocks;

	numBlocks = s->numBlocks;
	for (i=0; i < numBlocks; ++i)
		freeBlock(s->blocks[i]);
	free(s->blocks);
	free(s);
}

static System * cloneSystem(GAindividual X)
{
	System * s = (System*)X;
	int i;
	int numBlocks;
	System * s2;

	numBlocks = s->numBlocks;
	s2 = (System*)malloc(sizeof(System));
	
	s2->blocks = (Block**)malloc(numBlocks * sizeof(Block*));
	
	for (i=0; i < numBlocks; ++i)
		s2->blocks[i] = copyBlock(s->blocks[i]);

	s2->numBlocks = s->numBlocks;
	s2->numSpecies = s->numSpecies;
	
	return s2;
}

/***********************************************************************
* crossover two systems by taking random blocks from each and rewiring the severed links
************************************************************************/

static int totalInternalSpecies(System * s)
{
	int i = 0, sz = 0;
	
	for (i=0; i < s->numBlocks; ++i)
		sz += numInternals(s->blocks[i]);
	
	return sz;
}

static void pruneSystem(System * s) //find unconnected inputs/outputs
{
	int i,j,k,n,b1,b2,total,maxi;
	int * isUsed;
	total = s->numSpecies;
	
	isUsed = (int*)malloc(total * sizeof(int));
	
	for (i=0; i < s->numSpecies; ++i)
	{
		b1 = b2 = 0;
		for (j=0; j < s->numBlocks; ++j)
		{
			n = numInputs(s->blocks[j]);
			for (k=0; k < n; ++k)
				if (s->blocks[j]->inputs[k] == i)
				{
					b1 = 1;
					break;
				}
			n = numOutputs(s->blocks[j]);
			for (k=0; k < n; ++k)
				if (s->blocks[j]->outputs[k] == i)
				{
					b2 = 1;
					break;
				}
		}
		
		k = (int)(b1 && b2);
		isUsed[i] = k;
		total -= (1 - k);
	}
	
	maxi = 0;
	for (i=0; i < s->numSpecies; ++i)
	{
		for (j=0; j < s->numBlocks; ++j)
		{
			n = numInputs(s->blocks[j]);
			for (k=0; k < n; ++k)
				if (s->blocks[j]->inputs[k] == i)
					s->blocks[j]->inputs[k] = isUsed[i] ? maxi : -1;
					
			n = numOutputs(s->blocks[j]);
			for (k=0; k < n; ++k)
				if (s->blocks[j]->outputs[k] == i)
					s->blocks[j]->outputs[k] = isUsed[i] ? maxi : -1;
		}
		
		maxi += isUsed[i];
	}
	s->numSpecies = total;
}

static System * randomSubsystem(System * s, double prob) //get random subset of blocks
{
	int i, i2;
	int numBlocks = s->numBlocks, numBlocks2;
	System * s2 = (System*)malloc(sizeof(System));
	
	numBlocks2 = (int)(prob * (double)numBlocks);
	
	s2->blocks = (Block**)malloc(numBlocks2 * sizeof(Block*));
	
	for (i=0, i2=0; i2 < numBlocks2; ++i)
	{
		if (i >= numBlocks)
			i = 0;
		if (mtrand() < prob)
		{
			s2->blocks[i2] = copyBlock(s->blocks[i]);
			++i2;
		}
	}

	s2->numBlocks = numBlocks2;
	s2->numSpecies = s->numSpecies;
	
	pruneSystem(s2);
	return s2;
}

static void printB(Block * block, int n)
{
	int i;
	
	int in = numInputs(block),
		out = numOutputs(block);
	
	for (i=0; i < in; ++i)
		printf("%i ",block->inputs[i]);
	
	printf(" -- [%i] -- ",n);
	
	for (i=0; i < out; ++i)
		printf("%i ",block->outputs[i]);

	printf("\n");
}

static void reassignInputsOutputs(Block * block, int numSpecies)
{
	int i,j,n1,n2,b;

	n1 = numInputs(block);
	n2 = numOutputs(block);

	if (numSpecies <= n2) return;

	for (i=0; i < n1; ++i)
	{
		b = 1;
		while (b)
		{
			b = 0;

			for (j=0; j < n2; ++j)
				if ((block->inputs[i] == block->outputs[j]))
				{
					b = 1;
					break;
				}

			if (b == 1)
				block->inputs[i] = (int)(mtrand() * numSpecies);
		}
	}
}

static void printB(Block * block, int n);

static GAindividual * crossoverBlocks(GAindividual X, GAindividual Y) //place two subsets together
{
	System * s1 = randomSubsystem((System*)X, 0.6);
	System * s2 = randomSubsystem((System*)Y, 0.6);
	System * s3;
	
	int sz1 = s1->numSpecies,
		sz2 = s2->numSpecies;
	int i,j,n;
	
	for (i=0; i < s1->numBlocks; ++i)
	{
		n = numInputs(s1->blocks[i]);
		for (j=0; j < n; ++j)
			if (s1->blocks[i]->inputs[j] < 0)
				s1->blocks[i]->inputs[j] = sz1 + (int)(mtrand() * sz2);
		
		n = numOutputs(s1->blocks[i]);
		for (j=0; j < n; ++j)
			if (s1->blocks[i]->outputs[j] < 0)
				s1->blocks[i]->outputs[j] = sz1 + (int)(mtrand() * sz2);
	}

	for (i=0; i < s2->numBlocks; ++i)
	{
		n = numInputs(s2->blocks[i]);
		for (j=0; j < n; ++j)
			if (s2->blocks[i]->inputs[j] < 0)
				s2->blocks[i]->inputs[j] = (int)(mtrand() * sz1);
			else
				s2->blocks[i]->inputs[j] += sz1;

		n = numOutputs(s2->blocks[i]);
		for (j=0; j < n; ++j)
			if (s2->blocks[i]->outputs[j] < 0)
				s2->blocks[i]->outputs[j] = (int)(mtrand() * sz1);
			else
				s2->blocks[i]->outputs[j] += sz1;
	}
	
	s3 = (System*)malloc(sizeof(System));
	s3->numBlocks = s1->numBlocks + s2->numBlocks;
	s3->numSpecies = s1->numSpecies + s2->numSpecies;
	s3->blocks = (Block**)malloc(s3->numBlocks * sizeof(Block));

	for (i=0; i < s1->numBlocks; ++i)
	{
		s3->blocks[i] = s1->blocks[i];
		printB(s3->blocks[i],i);
	}

	for (i = 0; i < s2->numBlocks; ++i)
		s3->blocks[i + s1->numBlocks] = s2->blocks[i];

	s1->blocks = 0;
	s2->blocks = 0;
	free(s1);
	free(s2);

	return (void*)s3;
}

/*******************
* mutation
********************/

static Block * randomBlock()
{
	int i,k1,k2;
	int n = numBlockTypes();
	double d1,d2;
	Block * block = (Block*)malloc(sizeof(Block));
	k1 = (int)(mtrand()*n);

	if (!ALLOWED_BLOCKS) initialzeArrays();

	k2 = k1;
	while (!ALLOWED_BLOCKS[k1])
	{
		++k1;
		if (k1 >= n) k1 = 0;
		if (k1 == k2) break;
	}

	block->type = k1;

	block->inputs = (int*)malloc( numInputs(block) * sizeof(int) );
	block->outputs = (int*)malloc( numOutputs(block) * sizeof(int) );
	block->internals = (int*)malloc( numInternals(block) * sizeof(int) );

	n = numParams(block);
	block->params = (double*)malloc( n * sizeof(double) );

	for (i=0; i < n; ++i)
	{
		d1 = paramLowerBound(block,i);
		d2 = paramUpperBound(block,i);
		block->params[i] = d1 + mtrand() * (d2 - d1);
	}

	return block;
}

static GAindividual mutateBlocksH(GAindividual X)  //rewire or change parameter
{
	System * S = (System*)X;
	int i,j,k,k0,k1,m0,m1,n;
	
	int numBlocks = S->numBlocks,
		type = 0;
	
	Block ** oldBlocks, *block;
	
	if (mtrand() < PROB_PARAM_CHANGE) //mutate parameter
	{	
		k0 = (int)(mtrand() * numBlocks);
		k1 = k0;
		
		while (1)
		{	
			n = numParams(S->blocks[k1]);
			m0 = (int)(mtrand() * n);
			type = S->blocks[k1]->type;
			m1 = m0;
			while (FIXED_PARAMS[type][k1])
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
			if (k1 >= numBlocks) k1 = 0;
			if (k0==k1) break;
		}
	}
	
	//add or remove block?
	if ((mtrand() < (PROB_NEW_BLOCK + PROB_DEL_BLOCK)) && (S->numBlocks > MIN_SIZE || S->numBlocks < MAX_SIZE))
	{
		if (mtrand() < (PROB_NEW_BLOCK/(PROB_NEW_BLOCK + PROB_DEL_BLOCK)) && S->numBlocks < MAX_SIZE)
		{
			oldBlocks = S->blocks;
			S->blocks = (Block**)malloc((S->numBlocks+1)*sizeof(Block*));
			for (i=0; i < S->numBlocks; ++i)
				S->blocks[i] = oldBlocks[i];
			free(oldBlocks);
			
			block = randomBlock();

			n = numInputs(block);
			for (i=0; i < n; ++i)
				block->inputs[i] = (int)(mtrand() * S->numSpecies);

			n = numOutputs(block);
			for (i=0; i < n; ++i)
				block->outputs[i] = (int)(mtrand() * S->numSpecies);

			if (NO_SAME_INPUT_OUTPUT)
				reassignInputsOutputs(block,S->numSpecies);

			S->blocks[S->numBlocks+1] = block;
			S->numBlocks++;
		}
		else
			if (S->numBlocks > MIN_SIZE)
			{
				k = (int)(mtrand() * S->numBlocks);
				oldBlocks = S->blocks;
				S->blocks = (Block**)malloc((S->numBlocks-1)*sizeof(Block*));
				for (i=0,j=0; i < S->numBlocks; ++i)
					if (i != k)
					{
						S->blocks[j] = oldBlocks[i];
						++j;
					}
					else
						freeBlock(oldBlocks[i]);
				free(oldBlocks);
				
				S->numBlocks--;
				pruneSystem(S);

				for (i=0; i < S->numBlocks; ++i)
				{
					block = S->blocks[i];
					n = numInputs(block);
					for (i=0; i < n; ++i)
						if (block->inputs[i] < 0)
							block->inputs[i] = (int)(mtrand() * S->numSpecies);
	
					n = numOutputs(block);
					for (i=0; i < n; ++i)
						if (block->outputs[i] < 0)
							block->outputs[i] = (int)(mtrand() * S->numSpecies);
					if (NO_SAME_INPUT_OUTPUT)
						reassignInputsOutputs(block,S->numSpecies);
				}
			}
	}
	
	//rewire
	
	if (mtrand() < 0.5) //input or output
	{	
		k0 = (int)(mtrand() * numBlocks);
		k1 = k0;
		
		while (1)
		{	
			n = numInputs(S->blocks[k1]);
			m0 = (int)(mtrand() * n);
			type = S->blocks[k1]->type;
			m1 = m0;
			while (FIXED_INPUTS[type][m1])
			{
				++m1;
				if (m1 >= n) m1 = 0; //cycle
				if (m1 == m0) break;
			}

			++k1;
			if (k1 >= numBlocks) k1 = 0;
			if (k0==k1) break;
			if (!FIXED_INPUTS[type][m1])
			{
				S->blocks[k1]->inputs[m1] = (int)(mtrand() * S->numSpecies);
				if (NO_SAME_INPUT_OUTPUT)
					reassignInputsOutputs(S->blocks[k1],S->numSpecies);
				return (GAindividual)S;
			}
		}
	}
	else
	{
		k0 = (int)(mtrand() * numBlocks);
		k1 = k0;
		
		while (1)
		{	
			n = numOutputs(S->blocks[k1]);
			m0 = (int)(mtrand() * n);
			type = S->blocks[k1]->type;
			m1 = m0;
			while (FIXED_OUTPUTS[type][m1])
			{
				++m1;
				if (m1 >= n) m1 = 0; //cycle
				if (m1 == m0) break;
			}

			++k1;
			if (k1 >= numBlocks) k1 = 0;
			if (k0==k1) break;
			if (!FIXED_OUTPUTS[type][m1])
			{
				S->blocks[k1]->outputs[m1] = (int)(mtrand() * S->numSpecies);
				if (NO_SAME_INPUT_OUTPUT)
					reassignInputsOutputs(S->blocks[k1],S->numSpecies);
				return (GAindividual)S;
			}
		}
	}
	
	return (GAindividual)S;
}

static GAindividual mutateBlocks(GAindividual X)  //rewire or change parameter n times
{
	int i;
	
	if ((PROB_PARAM_CHANGE + PROB_REWIRE) > 0.0)
		for (i=0; i < MUTATION_RATE; ++i)
			mutateBlocksH(X);

	return X;
}

/**************************
* main evolution function
***************************/

static System * randomSystem(int numBlocks, int numSpecies)
{
	int i,j,n;
	System * S = (System*)malloc(sizeof(System));
	S->blocks = (Block**)malloc(numBlocks * sizeof(Block*));

	S->numBlocks = numBlocks;
	S->numSpecies = numSpecies;

	for (i=0; i < numBlocks; ++i)
	{
		S->blocks[i] = randomBlock();

		n = numInputs(S->blocks[i]);
		for (j=0; j < n; ++j)
			S->blocks[i]->inputs[j] = (int)(mtrand() * numSpecies);

		n = numOutputs(S->blocks[i]);
		for (j=0; j < n; ++j)
			S->blocks[i]->outputs[j] = (int)(mtrand() * numSpecies);

		if (NO_SAME_INPUT_OUTPUT)
			reassignInputsOutputs(S->blocks[i],numSpecies);
	}

	return S;
}

GApopulation evolveNetworks(int initialPopulationSize, int finalPopulationSize, GAFitnessFunc fitness, GACallbackFunc callback)
{
	return 0;
}

int main(int args, char** argv)
{
	int i;
	initMTrand();
	
	System * S1 = randomSystem(5,5);
	System * S2 = randomSystem(6,6);
	System * S3 = 0;
	
	for (i=0; i < S1->numBlocks; ++i)
		printB(S1->blocks[i],i);
	
	printf("\n\n");
	
	for (i=0; i < S2->numBlocks; ++i)
		printB(S2->blocks[i],i);
	
	printf("\n\n");
	
	for (i=0; i<2; ++i)
	{
		mutateBlocks(S1);
		mutateBlocks(S2);
	}
	
	printf("\nMutants\n");
	
	for (i=0; i < S1->numBlocks; ++i)
		printB(S1->blocks[i],i);
	
	printf("\n\n");
	
	for (i=0; i < S2->numBlocks; ++i)
		printB(S2->blocks[i],i);
	
	S3 = (System*)crossoverBlocks(S1,S2);
	
	printf("\nnum species=%i\n",S3->numSpecies);
	
	for (i=0; i < S3->numBlocks; ++i)
		printB(S3->blocks[i],i);
	
	freeSystem(S3);
	freeSystem(S1);
	freeSystem(S2);

	return 0;
}
