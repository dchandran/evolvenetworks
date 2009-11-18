/*!
  \file    functions.h
  \author: Deepak Chandran (dchandran1@gmail.com)
  \brief   A framework for constructing and evolving (see ga.h) modular biological systems

   This file supplements the blocks implementation described in blocks.h
   This file defines the rates and stoichiometry matrix functions for the different 
   block types. The different block types are also constructed in this file.
   
   New block types can be added by addeding a new row to the BlockTypesTable
   at the end of this file. New functions need to be constructed calculating for the new 
   block type's rates and stoichiometry matrix. See some of the other block types
   listed in the BlockTypesTable as an example. See the definition of BlockType in blocks.h
   to see what constitutes a block type.
	
**/

#ifndef BLOCKTYPES_FUNCTIONS_H
#define BLOCKTYPES_FUNCTIONS_H

#include "blocks.h"

/******************************************************************
* simple mass action with combinations of 1 and 2 reactants and products
*******************************************************************/

static void print_ma(FILE * file, Block* block)
{
	
}

static void uniuni_stoic( Matrix * matrix, Block * block )
{
	valueAt(*matrix, block->inputs[0], 0) = -1.0;
	valueAt(*matrix, block->outputs[0], 0) = 1.0;
}

static void uniuni_rates( double t, double * rates, double * conc, Block * block )
{
	rates[0] = block->params[0] * conc[ block->inputs[0] ];
}

static void biuni_stoic( Matrix * matrix, Block * block )
{
	valueAt(*matrix, block->inputs[0], 0) = -1.0;
	valueAt(*matrix, block->inputs[1], 0) = -1.0;
	valueAt(*matrix, block->outputs[0], 0) = 1.0;
}

static void biuni_rates( double t, double * rates, double * conc, Block * block )
{
	rates[0] = block->params[0] * conc[ block->inputs[0] ] * conc[ block->inputs[1] ];
}

static void unibi_stoic( Matrix * matrix, Block * block )
{
	valueAt(*matrix, block->inputs[0], 0) = -1.0;
	valueAt(*matrix, block->outputs[0], 0) = 1.0;
	valueAt(*matrix, block->outputs[1], 0) = 1.0;
}

static void unibi_rates( double t, double * rates, double * conc, Block * block )
{
	rates[0] = block->params[0] * conc[ block->inputs[0] ];
}

static void bibi_stoic( Matrix * matrix, Block * block )
{
	valueAt(*matrix, block->inputs[0], 0) = -1.0;
	valueAt(*matrix, block->inputs[1], 0) = -1.0;
	valueAt(*matrix, block->outputs[0], 0) = 1.0;
	valueAt(*matrix, block->outputs[1], 0) = 1.0;
}

static void bibi_rates( double t, double * rates, double * conc, Block * block )
{
	rates[0] = block->params[0] * conc[ block->inputs[0] ] * conc[ block->inputs[1] ];
}

static void print_enzyme(FILE * file, Block * block)
{
}

static void enzyme_stoic( Matrix * matrix, Block * block )
{
	valueAt(*matrix, block->inputs[0], 0) = -1.0;   //s + e -> es
	valueAt(*matrix, block->inputs[1], 0) = -1.0;
	valueAt(*matrix, block->internals[0], 0) = 1.0;
	
	valueAt(*matrix, block->inputs[0], 0) = 1.0;   //es -> s + e
	valueAt(*matrix, block->inputs[1], 0) = 1.0;
	valueAt(*matrix, block->internals[0], 0) = -1.0;
	
	valueAt(*matrix, block->outputs[0], 0) = 1.0;  //es -> p + e
	valueAt(*matrix, block->inputs[1], 0) = 1.0;
	valueAt(*matrix, block->internals[0], 0) = -1.0;
	
	valueAt(*matrix, block->outputs[0], 0) = -1.0;  //p + e -> es
	valueAt(*matrix, block->inputs[1], 0) = -1.0;
	valueAt(*matrix, block->internals[0], 0) = 1.0;
}

static void enzyme_rates( double t, double * rates, double * conc, Block * block )
{
	rates[0] = block->params[0] * conc[ block->inputs[0] ] * conc[ block->inputs[1] ];
	rates[1] = block->params[1] * conc[ block->internals[0] ];
	rates[2] = block->params[2] * conc[ block->internals[0] ];
	rates[3] = block->params[3] * conc[ block->inputs[0] ] * conc[ block->outputs[0] ];
}

/*! \brief  This lookup table is used to get the details of a Block type from its index, or type. 
	Each Block contains a "type" field, which is an index in this array. 
	Each row defines a Block type using the following information:
	 name,
	 stoichiometry function,
	 rates function,
	 num. reactions,
	 num. inputs,
	 num. outputs,
	 num. internal variables,
	 num. parameters,
	 (optional) parameter lower limits,
	 (optional) parameter upper limits
*/
BlockType BlockTypesTable[] =
{
	{"single reactant and product\0", &uniuni_stoic, 	&uniuni_rates,	&print_ma,	1, 	1, 	1, 	0, 	1, 	0, 	0},
	{"dimerization\0", 				&biuni_stoic, 	&biuni_rates, 	&print_ma,	 	1, 	2, 	1, 	0, 	1, 	0, 	0},
	{"dissociation\0", 				&unibi_stoic, 	&unibi_rates, 	&print_ma,	 	1, 	1, 	2, 	0, 	1, 	0, 	0},
	{"two reactants and products\0", 	&bibi_stoic, 	&bibi_rates, 	&print_ma,	 	1, 	2, 	2, 	0, 	1, 	0, 	0},
	{"enzymatic reaction\0", 			&enzyme_stoic, 	&enzyme_rates, 	&print_enzyme,	 	4, 	2, 	1, 	1, 	4, 	0, 	0},
	{0,0,0,0,0,0,0,0,0,0} //NULL type to mark the end of array
};

#endif
