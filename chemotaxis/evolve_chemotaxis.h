/********************************************************************************************************

	This is a program that uses my genetic algorithm library to evolve networks that are able to mimic
	chemotaxis is cells.
	
	gcc evolve_chemotaxis.c protein_network.c ga.c mtrand.c ssa.c cvodesim.c loops.c -lm -lcvode
		
*********************************************************************************************************/

#ifndef EVOLVE_CHEMOTAXIS_GA
#define EVOLVE_CHEMOTAXIS_GA

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ga.h"

/*
 * uncomment one of the includes and one of the typedefs.
 * This allows other reaction network models 
 * to be used for evolving chemotaxis.
*/

#include "reactionNetwork.h"

/*!
 *\brief chemotaxis network is comprised of reaction network and co-ordinate information
*/
typedef struct
{
	ReactionNetwork * network;
	double x,y,angle,xmax,xmin,ymax,ymin;
}
chemotaxis_network;

/*!
 * \brief ODE function
 * \param double time
 * \param double vector of concentration values
 * \param double vector of derivatives values
 * \param void* network
*/

void chemotaxisODE(double time,double* u,double* du,void * p);
/*!
 * \brief SSA propensity function
 * \param double time
 * \param double vector of concentration values
 * \param double vector of rates
 * \param void* network
*/
void chemotaxisPropensity(double time,double* u,double* v,void * p);

#endif

