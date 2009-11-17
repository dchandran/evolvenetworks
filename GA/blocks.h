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

/*! \brief function type for getting the stoichiometry matrix of a block 
*	\param Matrix the output stoichiometry matrix
*	\param Block * the block
*/
typedef void (*StiochiometryFunction)(Matrix*,Block*);

/*! \brief function type for calculating the rate vector for a block at the given concentration values
*	\param double* vector of concentration values
*	\param double* output vector of rate values
*	\param Block * the block
*/
typedef void (*RatesFunction)(double*,double*,Block*);

/*! \brief a block type */
typedef struct
{
	char * name;
	StiochiometryFunction stoic;
	RatesFunction rates;
	int numReactions;
	int numInputs;
	int numOutputs;
	int numInternals;
	int numParameters;
	double * paramsLowerBound;
	double * paramsUpperBound;
}
BlockType;

/*! \brief a system is defined as a set of blocks and a set number of molecular species */
typedef struct
{
	Block * blocks;
	int numBlocks;
	int numSpecies;
}
SystemOfBlocks;

/*! \brief (used during mutation events) set the probability for mutating a parameter value */
void setMutateParameterProb(double);

/*! \brief (used during mutation events) set the probability for rewiring two blocks in a system*/
void setMutateStructureProb(double);

/*! \brief evolve a population of Systems for optimizing the given fitness function
*	\param int initial size of population (number of systems)
*	\param int final size of population (number of systems)
*	\param GAFitnessFunc the fitness function (see ga.h)
*	\param GACallbackFunc optional callback funtion (see ga.h)
*	\return GApopulation a set of Systems sorted by fitness, null terminated. use GAfree to remove.
*/
GApopulation evolveNetworks(int initialPopulationSize, int finalPopulationSize, GAFitnessFunc fitness, GACallbackFunc callback);

/*! \brief check is a block type is a null type
*	\return int 0 or 1*/
int isNullBlockType(BlockType);

/*! \brief check is a block's type is a null type
*	\return int 0 or 1*/
int isNullBlock(Block *);

/*! \brief get the number of total block types that exist
*	\return int number of block types*/
int numBlockTypes();

/*! \brief (used during initialization) set the type of blocked to be used in a system
*	\param char* name of the block type, e.g. "one to one"
*	\return int 0 or 1 indicating whether this block name exists*/
int addBlockType(const char* name);

/*! \brief (used during initialization) set the type of blocked to be used in a system
*	\param int index of the block type
*	\return int 0 or 1 indicating whether this block type exists*/
int addBlockTypeByIndex(int);

/*! \brief (used during initialization) disallow a type of blocked from being used in a system
*	\param char* name of the block type, e.g. "one to one"
*	\return int 0 or 1 indicating whether this block name exists*/
int removeBlockType(const char* name);

/*! \brief (used during initialization) disallow a type of blocked from being used in a system
*	\param int index of a block type
*	\return int 0 or 1 indicating whether this block type exists*/
int removeBlockTypeByIndex(int);

#endif
