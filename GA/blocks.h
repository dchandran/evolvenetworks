/*!
  \file    blocks.h
  \author: Deepak Chandran (dchandran1@gmail.com)
  \brief   A framework for constructing and evolving (see ga.h) modular biological systems

   This file defines the structs called Block and System. A Block represents a single biological process
   that is defined by a stoichiometry matrix and reaction rate equations. Each block contains unique
   parameter values, inputs, outputs, and internal molecules. The input and output molecules are used to connect
   two or more blocks to each other. Each block has a BlockType, which is a struct that stores the stoichiometry and 
   rate function for blocks of that type. It also stores the number of inputs and outputs and other information
   about all blocks of that type. 
   
   A System is composed of a set number of molecules and a set number of interacting blocks. By generating
   the rates and stoichiometry matrices for each block, the entire system can be reduced to a single
   dynamical system.
	
**/

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

/*! \brief function type for printing a block
*	\param FILE* output
*	\param Block * the block
*/
typedef void (*PrintFunction)(FILE*,Block*);


/*! \brief a block type */
typedef struct
{
	char * name;
	StiochiometryFunction stoic;
	RatesFunction rates;
	PrintFunction print;
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
	Block ** blocks;
	int numBlocks;
	int numSpecies;
}
SystemOfBlocks;

/*!
    \name main function
	\{
*/

/*! \brief evolve a population of Systems for optimizing the given fitness function
*	\param int initial size of population (number of systems)
*	\param int final size of population (number of systems)
*	\param GAFitnessFunc the fitness function (see ga.h)
*	\param GACallbackFunc optional callback funtion (see ga.h)
*	\return GApopulation a set of Systems sorted by fitness, null terminated. use GAfree to remove.
*/
GApopulation evolveNetworks(int initialPopulationSize, int finalPopulationSize, GAFitnessFunc fitness, GACallbackFunc callback);

/*! \}
    \name getting information about blocks
	\{
*/

/*! \brief number of inputs in the given block
*	\return int */
int numInputs(Block * block);

/*! \brief number of outputs in the given block
*	\return int */
int numOutputs(Block * block);

/*! \brief number of parameters in the given block
*	\return int */
int numParams(Block * block);

/*! \brief number of internal molecules in the given block
*	\return int */
int numInternals(Block * block);

/*! \brief number of reactions in the given block
*	\return int */
int numReactions(Block * block);

/*! \brief check is a block type is a null type
*	\return int 0 or 1*/
int isNullBlockType(BlockType);

/*! \brief check is a block's type is a null type
*	\return int 0 or 1*/
int isNullBlock(Block *);

/*! \brief get the number of total block types that exist
*	\return int number of block types*/
int numBlockTypes();

/*! \brief get the index of the block type with given name
*	\param const char* name
*	\return int */
int getBlockTypeIndex(const char* name);

/*! \brief get the name of the block type with given index
*	\param int index
*	\return char* name */
const  char* getBlockTypeName(int);

/*! \}
    \name evolution settings
	\{
*/

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

/*! \brief (used during mutation events) set the probability for mutating a parameter value
*	\param double probability
*/
void setMutateParameterProb(double);

/*! \brief (used during mutation events) set the probability for adding a new block. 
*	Important: this probability is multiplied by (1 - parameter mutation prob.)
*	\param double probability
*/
void setAddBlockProb(double);

/*! \brief (used during mutation events) set the probability for removing a new block 
*	Important: this probability is multiplied by (1 - parameter mutation prob.)
*	\param double probability
*/
void setRemoveBlockProb(double);

/*! \brief size limits of a system (size = number of blocks)
*	\param int min
*	\param int max
*/
void setSizeRange(int min, int max);

/*! \brief (used during mutation events) set the probability for rewiring two blocks in a system
*	\param double probability
*/
void setRewiringProb(double);

/*! \brief (used during mutation events) allow for mutating of parameter values */
void allowParameterChange();

/*! \brief (used during mutation events) disallow for mutating of parameter values */
void disallowParameterChange();

/*! \brief (used during mutation events) allow for mutating of all parameter values for a particular block type 
*	\param const char* name of block type
*/
void allowParameterChangeFor(const char* name);

/*! \brief (used during mutation events) disallow for mutating of all parameter values for a particular block type
*	\param const char* name of block type
*/
void disallowParameterChangeFor(const char* name);

/*! \brief (used during mutation events) fix a particular parameter of a particular block type (i.e. it will not change). 
	Note: calling allowParameterChange() after this function will remove the fixed-ness
*	\param const char* name of block type
* 	\param int parameter index
*	\param double value of parameter
*/
void fixParameter(const char* name, int param, double value);

/*! \brief (used during mutation events) allow rewiring */
void allowRewiring();

/*! \brief (used during mutation events) disallow rewiring */
void disallowRewiring();

/*! \brief (used during mutation events) allow rewiring for a particular block type
*	\param const char* name of block type
*/
void allowRewiringFor(const char *);

/*! \brief (used during mutation events) disallow rewiring for a particular block type
 *	\param const char* name of block type
*/
void disallowRewiringFor(const char *);

/*! \brief (used during mutation events) disallow rewiring for a particular input in a particular block type
 *	\param const char* name of block type
 *	\param int index of input
 *	\param int 1=fix 0=free to rewire
*/
void setInputFixed(const char * name, int input, int fixed);

/*! \brief (used during mutation events) disallow rewiring for a particular input in a particular block type
 *	\param const char* name of block type
 *	\param int index of input
 *	\param int 1=fix 0=free to rewire
*/
void setOutputFixed(const char * name, int input, int fixed);

/*! \brief (used during mutation events) set the number of changes to the system that occur during a single mutation event
*	\param int must be positive
*/
void setMutationRate(int);

/*! \brief (used during mutation events) allow the same molecule to act as input and output for the same block*/
void allowSameInputAndOutput();

/*! \brief (used during mutation events) disallow the same molecule to act as input and output for the same block*/
void disallowSameInputAndOutput();
#endif
