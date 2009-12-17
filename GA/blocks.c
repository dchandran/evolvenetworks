#include <stdlib.h>
#include <stdio.h>
#include "blocks.h"
#include "blocksTable.h"

static int MIN_SIZE = 1;  //minimum allowed size of a system (size = num. of blocks)
static int MAX_SIZE = 20;  //maximum allowed size of a system (size = num. of blocks)
static int PERCENT_OVERLAP = 0.2; //expected percent of two blocks that will overlap

static double PROB_PARAM_CHANGE = 0.5; //prob. of changing parameter (during mutation)
static double PROB_REWIRE = 0.3;  //prob. of rewiring vs. changing parameter (during mutation)
static double PROB_NEW_BLOCK = 0.1;  //prob. of adding new block (during mutation)
static double PROB_DEL_BLOCK = 0.1;  //prob. of removing a block (during mutation)
static int MUTATION_RATE = 1; //number of changes made during a single mutation event. Cannot be < 1

static int** FIXED_PARAMS = 0; //indicates whether or not to mutate specific parameters
static int** FIXED_IO = 0; //indicates whether or not to rewire specific inputs/outputs
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
	FIXED_PARAMS = (int**)malloc(total * sizeof(int*));
	FIXED_IO = (int**)malloc(total * sizeof(int*));
	ALLOWED_BLOCKS = (int*)malloc(total * sizeof(int));

	for (i=0; i < total; ++i)
	{
		ALLOWED_BLOCKS[i] = 1; //all blocks are allowed by default

		FIXED_PARAMS[i] = (int*)malloc( BlockTypesTable[i].numParams * sizeof(int) );
		FIXED_IO[i] = (int*)malloc( BlockTypesTable[i].numExternals * sizeof(int) );

		//nothing fixed by default
		for (j=0; j < BlockTypesTable[i].numParams; ++j)
		{
			FIXED_PARAMS[i][j] = 0;
		}

		for (j=0; j < BlockTypesTable[i].numExternals; ++j)
			FIXED_IO[i][j] = 0;
	}
}

void clearArrays()
{
    int total, i, j;

	total = numBlockTypes();

	for (i=0; i < total; ++i)
	{
		ALLOWED_BLOCKS[i] = 1; //all blocks are allowed by default

        if (FIXED_PARAMS[i])
            free (FIXED_PARAMS[i]);

        if (FIXED_IO[i])
            free (FIXED_IO[i]);
	}

	if (FIXED_PARAMS)
        free(FIXED_PARAMS);

	if (FIXED_IO)
        free(FIXED_IO);

	if (ALLOWED_BLOCKS)
        free(ALLOWED_BLOCKS);
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

int numSpecies(Block * block)
{
    return (numExternals(block) + numInternals(block));
}

int numExternals(Block * block)
{
	return BlockTypesTable[ block->type ].numExternals;
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

void fixParameter(const char* name, int j)
{
	int i = getBlockTypeIndex(name);

	if (i < 0 || j < 0 || (j >= BlockTypesTable[i].numParams)) return;

	if (!ALLOWED_BLOCKS) initialzeArrays();

	FIXED_PARAMS[i][j] = 1;
}

void allowRewiring()
{
	int total, i, j;
	total = numBlockTypes();

	if (!ALLOWED_BLOCKS) initialzeArrays();

	for (i=0; i < total; ++i)
	{
		for (j=0; j < BlockTypesTable[i].numExternals; ++j)
				FIXED_IO[i][j] = 0;
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
		for (j=0; j < BlockTypesTable[i].numExternals; ++j)
				FIXED_IO[i][j] = 1;
	}

	PROB_PARAM_CHANGE = 1.0;
	PROB_REWIRE = 0.0;
}

void allowRewiringFor(const char * name)
{
	int j, i = getBlockTypeIndex(name);

	if (i < 0) return;

	if (!ALLOWED_BLOCKS) initialzeArrays();

	for (j=0; j < BlockTypesTable[i].numExternals; ++j)
		FIXED_IO[i][j] = 0;
}

void disallowRewiringFor(const char * name)
{
	int j, i = getBlockTypeIndex(name);

	if (i < 0) return;

	if (!ALLOWED_BLOCKS) initialzeArrays();

	for (j=0; j < BlockTypesTable[i].numExternals; ++j)
		FIXED_IO[i][j] = 1;
}

void setExternalFixed(const char * name, int j, int fixed)
{
	int i = getBlockTypeIndex(name);

	if (i < 0 || j < 0 || j > BlockTypesTable[i].numExternals) return;

	if (!ALLOWED_BLOCKS) initialzeArrays();

	FIXED_IO[i][j] = fixed;
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

void setPercentOverlap(double p)
{
    if (p > 0.0 && p < 1.0)
        PERCENT_OVERLAP = p;
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

void setParamUpperBound(Block * block, int k, double d)
{
	if (!block || block->type > numBlockTypes()) return;

	if (BlockTypesTable[block->type].paramsUpperBound &&
		k >= 0 &&
		k < BlockTypesTable[block->type].numParams)
		BlockTypesTable[block->type].paramsUpperBound[k] = d;
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

void setParamLowerBound(Block * block, int k, double d)
{
	if (!block || block->type > numBlockTypes()) return;

	if (BlockTypesTable[block->type].paramsLowerBound &&
		k >= 0 &&
		k < BlockTypesTable[block->type].numParams)
		BlockTypesTable[block->type].paramsLowerBound[k] = d;
}

/*********************************************
* for use by GA (see ga.h)
**********************************************/

static void freeBlock(Block * block)
{
	if (!block) return;

	if (block->internals)
		free(block->internals);

	if (block->externals)
		free(block->externals);

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

	n = numExternals(block);
	block2->externals = (int*)malloc(n * sizeof(int));
	for (i=0; i < n; ++i)
		block2->externals[i] = block->externals[i];

	n = numParams(block);
	block2->params = (double*)malloc(n * sizeof(double));
	for (i=0; i < n; ++i)
		block2->params[i] = block->params[i];

    n = numSpecies(block);
	block2->initVals = (double*)malloc(n * sizeof(double));
	for (i=0; i < n; ++i)
		block2->initVals[i] = block->initVals[i];

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

static GAindividual cloneSystem(GAindividual X)
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

	return (GAindividual)s2;
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
			n = numExternals(s->blocks[j]);
			for (k=0; k < n; ++k)
				if (s->blocks[j]->externals[k] == i)
				{
					b1 = 1;
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
			n = numExternals(s->blocks[j]);
			for (k=0; k < n; ++k)
				if (s->blocks[j]->externals[k] == i)
					s->blocks[j]->externals[k] = isUsed[i] ? maxi : -1;
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

static void reassignInputsOutputs(Block * block, int numSpecies)
{
	int i,j,n,b;

	n = numExternals(block);

	if (numSpecies <= n) return;

	for (i=0; i < n; ++i)
	{
		b = 1;
		while (b)
		{
			b = 0;

			for (j=(i+1); j < n; ++j)
				if ((block->externals[i] == block->externals[j]))
				{
					b = 1;
					break;
				}

			if (b == 1)
				block->externals[j] = (int)(mtrand() * numSpecies);
		}
	}
}

static GAindividual crossoverBlocks(GAindividual X, GAindividual Y) //place two subsets together
{
	System * s1 = randomSubsystem((System*)X, 0.6);
	System * s2 = randomSubsystem((System*)Y, 0.6);
	System * s3;

	int sz1 = s1->numSpecies,
		sz2 = s2->numSpecies;
	int i,j,n;

	for (i=0; i < s1->numBlocks; ++i)
	{
		n = numExternals(s1->blocks[i]);
		for (j=0; j < n; ++j)
			if (s1->blocks[i]->externals[j] < 0)
				s1->blocks[i]->externals[j] = sz1 + (int)(mtrand() * sz2);
	}

	for (i=0; i < s2->numBlocks; ++i)
	{
		n = numExternals(s2->blocks[i]);
		for (j=0; j < n; ++j)
			if (s2->blocks[i]->externals[j] < 0)
				s2->blocks[i]->externals[j] = (int)(mtrand() * sz1);
			else
				s2->blocks[i]->externals[j] += sz1;
	}

	s3 = (System*)malloc(sizeof(System));
	s3->numBlocks = s1->numBlocks + s2->numBlocks;
	s3->numSpecies = s1->numSpecies + s2->numSpecies;
	s3->blocks = (Block**)malloc(s3->numBlocks * sizeof(Block*));

	for (i=0; i < s1->numBlocks; ++i)
	{
		s3->blocks[i] = s1->blocks[i];
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

	k1 = numInternals(block);
	k2 = numExternals(block);

	block->externals = (int*)malloc( k1 * sizeof(int) );
	block->internals = (int*)malloc( k2 * sizeof(int) );
	block->initVals = (double*)malloc( (k1 + k2) * sizeof(double));

	n = numParams(block);
	block->params = (double*)malloc( n * sizeof(double) );

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

			n = numExternals(block);
			for (i=0; i < n; ++i)
				block->externals[i] = (int)(mtrand() * S->numSpecies);

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
					n = numExternals(block);
					for (i=0; i < n; ++i)
						if (block->externals[i] < 0)
							block->externals[i] = (int)(mtrand() * S->numSpecies);
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
			n = numExternals(S->blocks[k1]);
			m0 = (int)(mtrand() * n);
			type = S->blocks[k1]->type;
			m1 = m0;
			while (FIXED_IO[type][m1])
			{
				++m1;
				if (m1 >= n) m1 = 0; //cycle
				if (m1 == m0) break;
			}

			++k1;
			if (k1 >= numBlocks) k1 = 0;
			if (k0==k1) break;
			if (!FIXED_IO[type][m1])
			{
				S->blocks[k1]->externals[m1] = (int)(mtrand() * S->numSpecies);
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
			n = numExternals(S->blocks[k1]);
			m0 = (int)(mtrand() * n);
			type = S->blocks[k1]->type;
			m1 = m0;
			while (FIXED_IO[type][m1])
			{
				++m1;
				if (m1 >= n) m1 = 0; //cycle
				if (m1 == m0) break;
			}

			++k1;
			if (k1 >= numBlocks) k1 = 0;
			if (k0==k1) break;
			if (!FIXED_IO[type][m1])
			{
				S->blocks[k1]->externals[m1] = (int)(mtrand() * S->numSpecies);
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

/****************************
* simulation functions
****************************/

Matrix getStoichiometryMatrix(System * S)
{
	Matrix N;
	int numSpecies = 0, numReacs = 0, numInt = 0, numExt = 0;
	double * values = 0;
	int i,j,k,n;
	StiochiometryFunction f;

	N.rownames = N.colnames = 0;
	N.values = 0;

	for (i = 0; i < S->numBlocks; ++i)
	{
		numInt += numInternals(S->blocks[i]);
		numReacs += numReactions(S->blocks[i]);
	}

	k = S->numSpecies;
	for (i = 0; i < S->numBlocks; ++i)
	{
		n = numInternals(S->blocks[i]);
		for (j = 0; j < n; ++j, ++k)
			S->blocks[i]->internals[j] = k;
	}

	numSpecies = S->numSpecies + numInt;

	N.cols = numReacs;
	N.rows = numSpecies;
	N.values = (double*)malloc(numReacs * numSpecies * sizeof(double));
	for (i = 0; i < (numReacs*numSpecies); ++i)
		N.values[i] = 0.0;
	k = 0;

	for (i = 0; i < S->numBlocks; ++i)
	{
		f = BlockTypesTable[ S->blocks[i]->type ].stoic;
		f(&N,S->blocks[i],k);
		k += BlockTypesTable[ S->blocks[i]->type ].numReactions;
	}

	return N;
}

void getRates(double time, double* conc, double* rates, void * s)
{
	int i,k,n;
	RatesFunction f;
	System * S = (System*)s;

	k = 0;
	for (i = 0; i < S->numBlocks; ++i)
	{
		f = BlockTypesTable[ S->blocks[i]->type ].rates;
		f(time, conc, (rates + k), S->blocks[i]);
		k += BlockTypesTable[ S->blocks[i]->type ].numReactions;
	}
}

/*graphviz format*/
void printSystem(GAindividual s, FILE * fp)
{
    int i,j,n;
    System * S = (System*)s;

    for (i=0; i < S->numBlocks; ++i)
    {
        n = numExternals(S->blocks[i]);
        for (j=0; j < n; ++j)
            fprintf(fp,"M%i -- x%i",i,S->blocks[i]->externals[j]);
    }
}

void printSystemStats(GAindividual s,FILE * fp)
{
    System * S = (System*)s;
    fprintf(fp,"%i\t%i",S->numBlocks,S->numSpecies);
}

double * getInitialValues(System * S)
{
    int i,j,k1,k2,n,t;
    double * y = 0;

    n = S->numSpecies;

    for (i=0; i < S->numBlocks; ++i)
        n += numInternals(S->blocks[i]);

    y = (double*)malloc(n * sizeof(double));
    for (i=0; i < n; ++i)
        y[i] = 0.0;

    for (i=0, k2=0; i < S->numBlocks; ++i)
    {
        initializeBlock(S->blocks[i]);
        n = numExternals(S->blocks[i]);
        k1=0;

        for (j=0; j < n; ++j, ++k1)
            y[ S->blocks[i]->externals[j] ] = S->blocks[i]->initVals[k1];

        n = numInternals(S->blocks[i]);

        for (j=0; j < n; ++j, ++k1)
            y[ S->blocks[i]->internals[j] ] = S->blocks[i]->initVals[k1];
    }

    return y;
}

double * simulateStochastic(System * S, double time, int * sz)
{
	Matrix N = getStoichiometryMatrix(S);
	double * y, * y0;
	int n;

	y0 = getInitialValues(S);

	y = SSA(N.rows, N.cols, N.values , &getRates, y0, 0.0, time, 100000, sz, S);

    free(y0);
	free(N.values);

	return y;
}

double * simulateODE(System * S, double time, double dt)
{
	Matrix N = getStoichiometryMatrix(S);
	double * y, * y0;

	y0 = getInitialValues(S);

	y = ODEsim2(N.rows, N.cols, N.values , &getRates, y0, 0.0, time, dt, S);

    free(y0);
	free(N.values);

	return y;
}

/**************************
* main evolution function
***************************/

static System * randomSystem(int numBlocks)
{
	int i,j,n,numSpecies;
	System * S = (System*)malloc(sizeof(System));
	S->blocks = (Block**)malloc(numBlocks * sizeof(Block*));

	S->numBlocks = numBlocks;
	numSpecies = 0;

	for (i=0; i < numBlocks; ++i)
	{
		S->blocks[i] = randomBlock();
		numSpecies += numExternals(S->blocks[i]);
    }

    numSpecies *= (1.0 - PERCENT_OVERLAP);
    S->numSpecies = numSpecies;

    for (i=0; i < numBlocks; ++i)
	{
        n = numExternals(S->blocks[i]);

	    for (j=0; j < n; ++j)
			S->blocks[i]->externals[j] = (int)(mtrand() * numSpecies);

		if (NO_SAME_INPUT_OUTPUT)
			reassignInputsOutputs(S->blocks[i],numSpecies);
	}

	return S;
}

void initializeBlock(Block * block)
{
	BlockTypesTable[ block->type ].init(block);
}

void initializeSystem(System * S)
{
	int i;
	for (i=0; i < S->numBlocks; ++i)
		initializeBlock(S->blocks[i]);
}

GApopulation evolveNetworks(GAFitnessFunc fitness, int initialPopulationSize, int finalPopulationSize, int iter, GACallbackFunc callback)
{
    int i;
	GApopulation P;

	if (!ALLOWED_BLOCKS) initialzeArrays();

	P = (GApopulation)malloc(initialPopulationSize * sizeof(GAindividual));

	for (i=0; i < initialPopulationSize; ++i)
        P[i] = (GApopulation)randomSystem(MIN_SIZE + mtrand() * (MAX_SIZE - MIN_SIZE));

    printf("here1\n");
	GAinit( &freeSystem, &cloneSystem , fitness, &crossoverBlocks, &mutateBlocks, &GArouletteWheelSelection, callback);
	return GArun(P,initialPopulationSize,finalPopulationSize,iter);
}
