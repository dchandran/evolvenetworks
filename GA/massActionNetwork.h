/********************************************************************************************************

Copyright (C) 2009 Deepak Chandran

	This file defines a chemical reaction network using mass-action kinetics. The functions in this file
	are designed to be used with the GA library that I have written. 
	
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

#ifndef MASS_ACTION_FOR_GA
#define MASS_ACTION_FOR_GA

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ga.h"
#include "mtrand.h"
#include "cvodesim.h"   //for ODE simulation
#include "ssa.h"        //for stochastic simulation

/*! \brief Mass action network
* \ingroup massaction
*/
typedef struct
{
	int * r1;    //first reactant (can be -1 if none)
	int * r2;    //second reactant (can be -1 if none)
	int * p1;    //first product (can be -1 if none)
	int * p2;    //second product (can be -1 if none)
	double * k;  //rate constant for each reaction
	int reactions;    //number of reactions;
	int species;    //number of species;
	int * fixed;  //array of size=species. 1= ith species is fixed
} 
MassActionNetwork;

/*****************************************************
   @name  Functions needed by GA
******************************************************/

/*! \brief Free an individual from memory
 * \param GAindividual a single individual
 * \ingroup massaction
*/
void deleteMassActionNetwork(GAindividual individual);

/*! \brief Make a copy of an individual and return the memory pointer
 * \param GAindividual target individual
 * \ingroup massaction
*/
GAindividual cloneMassActionNetwork(GAindividual individual);

/*! \brief Combine two individuals to generate a new individual
 * \param GAindividual parent individual 1
 * \param GAindividual parent individual 2
 * \return GAindividual pointer to an individual (can be the same as one of the parents)
 * \ingroup massaction
*/
GAindividual crossoverMassActionNetwork(GAindividual individualA, GAindividual individualB);

/*! \brief Change an individual randomly to generate a new individual
 * \param GAindividual parent individual
 * \return GAindividual pointer to an individual (can be the same as one of the parents)
 * \ingroup massaction
*/
GAindividual mutateMassActionNetwork(GAindividual individual);

/*****************************************************
   @name  Functions for simulating and printing the network defined above
******************************************************/
/*! \brief Propensity function to be used by the SSA function (see ssa.h)
 * \param double time
 * \param double* values for variables in the system
 * \param double* rates (output)
 * \param GAindividual the network
 * \ingroup massaction
*/
void ratesForMassActionNetwork(double,double*,double*,GAindividual);

/*! \brief Get the stoichiometry matrix for a network
 * \param GAindividual network
 * \return double* linearized stoichiometry matrix
 * \ingroup massaction
*/
double * stoichiometryForMassActionNetwork(GAindividual);

/*! \brief Print a network
 * \param GAindividual network
 * \ingroup massaction
*/
void printMassActionNetwork(GAindividual);
/*! \brief get the number of variables in the network.
    This is equal to the number of rows in the stoichiometry matrix
 \param MassActionNetwork network
 \ingroup massaction
*/
int getNumSpeciesForMassActionNetwork(GAindividual);
/*! \brief get the number of reactions in the network. 
    This is equal to the number of columns in the stoichiometry matrix
 \param MassActionNetwork network
 \ingroup massaction
*/
int getNumReactionsForMassActionNetwork(GAindividual);
/*! \brief set a species as a fixed (constant, boundary) species
 \param MassActionNetwork network
 \param int index of species that should be set at fixed
 \param int value = 0 or 1, where 1 = fixed
 \ingroup massaction
*/
void setFixedSpeciesForMassActionNetwork(GAindividual, int,int);


/*****************************************************
    @name Functions for initializing a GA
******************************************************/

/*! \brief Set parameters for randomly generating mass-action networks
 * \param double percent of reactions having two reactants instead of one
 * \param double percent of reactions having two products instead of one
 * \param double percent of reactions representing creation or degradation (i.e. no products or reactants)
 * \param double average rate constant
 * \ingroup massaction
*/
void setParametersForMassActionNetwork(double , double , double, double );
/*! \brief Set parameters for randomly generating mass-action networks
 * \param int average number of species
 * \param int average number of reactions
 * \ingroup massaction
*/
void setSizeForMassActionNetwork(int,int);
/*! \brief Set parameters for the mutation and crossover functions
 * \param double probability of mutating a coefficient (default = 0.5) vs. adding/removing a reaction
 * \param double probability of crossover (default = 1.0, i.e. always)
 * \ingroup massaction
*/
void setMutationAndCrossoverRatesForMassActionNetwork(double, double);
/*! \brief Creates an array of randomized mass action networks.
 * \param int number of networks in population
 * \return GAindividual* GApopulation of random networks
 * \ingroup massaction
*/
GApopulation randomMassActionNetworks(int);
/*! \brief Make a new empty network
 * \param int number of species
 * \param int number of reactions
 * \return MassActionNetwork* network
 * \ingroup massaction
*/
MassActionNetwork * newMassActionNetwork(int,int);

#endif

