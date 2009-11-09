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


