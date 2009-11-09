#ifndef GA_BLOCK_BASED_EVOLUTION
#define GA_BLOCK_BASED_EVOLUTION

#include "ga.h"
//#include "cvodesim.h"
//#include "ssa.h"

/*! \brief a 2D matrix with rownames and colnames. 
The 2D matrix is stored as a single array for efficient memory management.
Use valueAt(M,i,j) macro to get i,j-th element of matrix M. */
typedef struct 
{
	int rows, cols;
	double * values;
	char** rownames;
	char** colnames;
} 
Matrix;

#ifndef valueAt
#define valueAt(M, i, j) ( M.values[ (i)*(M.cols) + (j) ] )
#endif

/*! \brief a block is one or more biological reactions representing some well defined process. */
typedef struct
{
	int type;
	int * inputs;
	int * outputs;
	int * internals;
	double * params;
}
Block;

typedef void (*StiochiometryFunction)(Matrix*,Block*);
typedef void (*RatesFunction)(double*,double*,Block*);

/*! \brief a block type. */
typedef struct
{
	char * name;
	StiochiometryFunction stoic;
	RatesFunction rates;
	int numReactions;
	int numExternals;
	int numInternals;
	int numParameters;
	double * paramsLowerBound;
	double * paramsUpperBound;
}
BlockType;

/*! \brief a system is defined as a set of blocks. */
typedef struct
{
	Block * blocks;
	int numBlocks;
	int numSpecies;
}
System;

void setMutateParameterProb(double);
void setMutateStructureProb(double);

void addBlockType(int);

GApopulation evolveNetworks(int initialPopulationSize, int finalPopulationSize, GAfitnessFunc fitness, GAcallback callback);

int isNullBlockType(BlockType);

int isNullBlock(Block *);

static BlockType BlockTypesTable[] =
{
	{"single reactant and product", &uniuni_stoic, &uniuni_rates, 1, 1, 1, 0, 1, 0, 0},
	{"dimerization", &biuni_stoic, &biuni_rates, 1, 2, 1, 0, 1, 0, 0},
	{"dissociation", &unibi_stoic, &unibi_rates, 1, 1, 2, 0, 1, 0, 0},
	{"two reactants and products", &bibi_stoic, &bibi_rates, 1, 2, 2, 0, 1, 0, 0},
	{"enzymatic reaction", &enzyme_stoic, &enzyme_rates, 4, 2, 1, 1, 4, 0, 0},
	{0,0,0,0,0,0,0,0,0,0} //NULL type
};

#endif
