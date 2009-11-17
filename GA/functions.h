#include "blocks.h"

static void uniuni_stoic( Matrix * matrix, Block * block )
{
	valueAt(matrix, block->inputs[0], 0) = -1.0;
	valueAt(matrix, block->outputs[0], 0) = 1.0;
}

static void uniuni_rates( double * rates, double * conc, Block * block )
{
	rates[0] = block->params[0] * conc[ block->inputs[0] ];
}

static void biuni_stoic( Matrix * matrix, Block * block )
{
	valueAt(matrix, block->inputs[0], 0) = -1.0;
	valueAt(matrix, block->inputs[1], 0) = -1.0;
	valueAt(matrix, block->outputs[0], 0) = 1.0;
}

static void biuni_rates( double * rates, double * conc, Block * block )
{
	rates[0] = block->params[0] * conc[ block->inputs[0] ] * conc[ block->inputs[1] ];
}

static void unibi_stoic( Matrix * matrix, Block * block )
{
	valueAt(matrix, block->inputs[0], 0) = -1.0;
	valueAt(matrix, block->outputs[0], 0) = 1.0;
	valueAt(matrix, block->outputs[1], 0) = 1.0;
}

static void unibi_rates( double * rates, double * conc, Block * block )
{
	rates[0] = block->params[0] * conc[ block->inputs[0] ];
}

static void bibi_stoic( Matrix * matrix, Block * block )
{
	valueAt(matrix, block->inputs[0], 0) = -1.0;
	valueAt(matrix, block->inputs[1], 0) = -1.0;
	valueAt(matrix, block->outputs[0], 0) = 1.0;
	valueAt(matrix, block->outputs[1], 0) = 1.0;
}

static void bibi_rates( double * rates, double * conc, Block * block )
{
	rates[0] = block->params[0] * conc[ block->inputs[0] ] * conc[ block->inputs[1] ];
}

static void enzyme_stoic( Matrix * matrix, Block * block )
{
	valueAt(matrix, block->inputs[0], 0) = -1.0;   //s + e -> es
	valueAt(matrix, block->inputs[1], 0) = -1.0;
	valueAt(matrix, block->internals[0], 0) = 1.0;
	
	valueAt(matrix, block->inputs[0], 0) = 1.0;   //es -> s + e
	valueAt(matrix, block->inputs[1], 0) = 1.0;
	valueAt(matrix, block->internals[0], 0) = -1.0;
	
	valueAt(matrix, block->outputs[0], 0) = 1.0;  //es -> p + e
	valueAt(matrix, block->inputs[1], 0) = 1.0;
	valueAt(matrix, block->internals[0], 0) = -1.0;
	
	valueAt(matrix, block->outputs[0], 0) = -1.0;  //p + e -> es
	valueAt(matrix, block->inputs[1], 0) = -1.0;
	valueAt(matrix, block->internals[0], 0) = 1.0;
}

static void enzyme_rates( double * rates, double * conc, Block * block )
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
	{"single reactant and product\0", &uniuni_stoic, 	&uniuni_rates, 	1, 	1, 	1, 	0, 	1, 	0, 	0},
	{"dimerization\0", 				&biuni_stoic, 	&biuni_rates, 	1, 	2, 	1, 	0, 	1, 	0, 	0},
	{"dissociation\0", 				&unibi_stoic, 	&unibi_rates, 	1, 	1, 	2, 	0, 	1, 	0, 	0},
	{"two reactants and products\0", 	&bibi_stoic, 	&bibi_rates, 	1, 	2, 	2, 	0, 	1, 	0, 	0},
	{"enzymatic reaction\0", 			&enzyme_stoic, 	&enzyme_rates, 	4, 	2, 	1, 	1, 	4, 	0, 	0},
	{0,0,0,0,0,0,0,0,0,0} //NULL type to mark the end of array
};