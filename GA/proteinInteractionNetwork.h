/********************************************************************************************************

Copyright (C) 2009 Deepak Chandran

	This file defines a chemical reaction network where each protein can switch from active to inactive form
	via Michaelis Menten type kinetics. The number of enzymes affecting each transition can be many.

	Example reaction:
		A --> A';  rate = vmax1 * A * P1 / (Km1 + A)
		A' --> A;  rate = vmax2 * (total - A)  * P2 /(Km2 + (total - A) )
 
	The functions in this file are designed to be used with the GA library that I have written. 
	
	This file defines the following functions that are required in the GA:
		1) a struct that represents an "individual"
		2) a mutation function to randomly alter an individual
		3) a crossover function to make a new individual from two individuals
	
	In addition, the file defines:
		1) ODE function to simulate a network
		2) Propensity function to simulate a network stochastically
	
	The program using this file must define:
		1) a function that returns the fitness of a network
		
*********************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ga.h"
#include "mtrand.h"
#include "cvodesim.h"   //for ODE simulation
#include "ssa.h"        //for stochastic simulation

#ifndef PROTEIN_INTERACTION_NETWORK_FOR_GA
#define PROTEIN_INTERACTION_NETWORK_FOR_GA

/*! \brief Struct used to store a set of regulators for one protein
* \ingroup proteinnetwork
*/
typedef struct
{
	int size;
	int * proteins; //index (equivalent to name) of each regulator
	double * Km;  //Michaelis-Menten constant for each regulator
	double * Vmax;  //Vmax for each regulator. Negative Vmax means negative regulator
}
Regulators;

/*! \brief Protein interaction network
* \ingroup proteinnetwork
*/
typedef struct
{
	int species;
	double * totals; //the total concentration of each protein
	Regulators * regulators; //the set of regulators for each protein
} 
ProteinInteractionNetwork;

/*****************************************************
   @name  Functions needed by GA
******************************************************/

/*! \brief Free an individual from memory
 * \param GAindividual a single individual
 * \ingroup proteinnetwork
*/
void deleteProteinInteractionNetwork(GAindividual individual);

/*! \brief Make a copy of an individual and return the memory pointer
 * \param GAindividual target individual
 * \ingroup proteinnetwork
*/
GAindividual cloneProteinInteractionNetwork(GAindividual individual);

/*! \brief Combine two individuals to generate a new individual
 * \param GAindividual parent individual 1
 * \param GAindividual parent individual 2
 * \return GAindividual pointer to an individual (can be the same as one of the parents)
 * \ingroup proteinnetwork
*/
GAindividual crossoverProteinInteractionNetwork(GAindividual individualA, GAindividual individualB);

/*! \brief Change an individual randomly to generate a new individual
 * \param GAindividual parent individual
 * \return GAindividual pointer to an individual (can be the same as one of the parents)
 * \ingroup proteinnetwork
*/
GAindividual mutateProteinInteractionNetwork(GAindividual individual);

/*****************************************************
   @name Functions for simulating and printing the network defined above
******************************************************/

/*! \brief Propensity function to be used by the SSA function (see ssa.h)
 * \param double time
 * \param double* values for variables in the system
 * \param double* rates (output)
 * \param GAindividual network
 * \ingroup proteinnetwork
*/
void ratesForProteinInteractionNetwork(double,double*,double*,GAindividual);
/*! \brief Get the stoichiometry matrix for a network
 * \param GAindividual network
 * \return double* linearized stoichiometry matrix
 * \ingroup proteinnetwork
*/
double * stoichiometryForProteinInteractionNetwork(GAindividual);
/*! \brief Print a network
 * \param GAindividual network
 * \ingroup proteinnetwork
*/
void printProteinInteractionNetwork(GAindividual);

/*****************************************************
   @name Functions for initializing a GA
******************************************************/
/*! \brief Set the initial (average) parameters for generating random networks.
 * \param double the average Ka value for random networks
 * \param double the average Vmax value for random networks
 * \param double the average Total (conservation law) for random networks
 * \ingroup proteinnetwork
*/
void setParametersForProteinInteractionNetwork(double, double , double );
/*! \brief Set the initial (average) parameters for generating random networks.
 * \param int average number of species in the networks
 * \param int regulations per species
 * \ingroup proteinnetwork
*/
void setSizeForProteinInteractionNetwork(int, int);
/*! \brief Set parameters for the mutation and crossover functions. First four arguments must add to 1.
 * \param double probability of rewiring the network during mutation
 * \param double probability of changing network parameter during mutation
 * \param double probability of changing a conservation law during mutation
 * \param double probability of adding/removing a node during mutation
 * \param double probability of crossover
 * \ingroup proteinnetwork
*/
void setMutationAndCrossoverRatesForProteinInteractionNetwork(double, double, double, double, double);
/*! \brief Creates an array of randomized mass action networks.
 * \return GAindividual* GApopulation of random networks
 * \ingroup proteinnetwork
*/
GApopulation randomProteinInteractionNetworks(int);

/*! \brief Make a new empty network
 * \param int number of species
 * \param int number of reactions
 * \return ProteinInteractionNetwork* network
 * \ingroup proteinnetwork
*/
ProteinInteractionNetwork * newForProteinInteractionNetwork(int,int);

#endif

