/*!
  \file    massActionNetwork.h
  \author: Deepak Chandran (dchandran1@gmail.com)
  \brief   evolve mass action networks

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
		
*/

#ifndef MASS_ACTION_FOR_GA
#define MASS_ACTION_FOR_GA

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ga.h"
#include "mtrand.h"
#include "cvodesim.h"   //for ODE simulation
#include "ssa.h"        //for stochastic simulation

/*! \brief 
Mass action network is defined using a set of reactions. 
Each reaction has a maximum of two reactants (reactant1 and reactant2) 
and a maximum of two products (product1 and product2). The reactant1 and reactant2
arrays store index values of the molecular species. The 
index value can range from 0 to (species-1), where (species) is
the number of molecules in this system. A value
of -1 is used to indicate no reactant or no product.
Each reaction also has a rate constant, k. The default reaction
rate is determined by the product of the rate constant and
the reactant concentrations. 
* \ingroup massaction
*/
typedef struct
{
	int * reactant1;    //first reactant (can be -1 if none)
	int * reactant2;    //second reactant (can be -1 if none)
	int * product1;    //first product (can be -1 if none)
	int * product2;    //second product (can be -1 if none)
	double * k;  //rate constant for each reaction
	int reactions;    //number of reactions;
	int species;    //number of species;
} 
MassActionNetwork;

/*!
   \name  Functions needed by GA
   \{
*/

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

/*!
   \}
   \name  Functions for simulating and printing the network defined above
   \{
*/
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
 * \param FILE* output stream
 * \param GAindividual network
 * \ingroup massaction
*/
void printMassActionNetwork(FILE *, GAindividual);
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

/*! \}
    @name Functions for initializing a GA
\{*/

/*! \brief Set probabilities for different reaction types mass-action networks. These probabilities will also
           be used by the mutation function
 * \param double percent of initial reactions having one reactant and one product
 * \param double percent of initial reactions having one reactant and two products
 * \param double percent of initial reactions having two reactants and one product
 * \param double percent of initial reactions having two reactant and two product
 * \param double percent of initial reactions representing creation (i.e. no reactants)
 * \param double percent of initial reactions representing degradation (i.e. no products)
 * \ingroup massaction
*/
void setDistributionOfMassActionNetwork(double uni_uni, double uni_bi, double bi_uni, double bi_bi, double no_reactant, double no_product);
/*! \brief range for the reaction rate constant when randomly generating mass-action networks and mutating
 * \param double min rate constant
 * \param double max rate constant
 * \ingroup massaction
*/
void setRateConstantForMassActionNetwork(double min_rate_constant, double max_rate_constant);
/*! \brief Set parameters for randomly generating mass-action networks
 * \param int min number of species
 * \param int max number of species
 * \param int min number of reactions
 * \param int max number of reactions
 * \ingroup massaction
*/
void setSizeForMassActionNetwork(int min_species,int max_species, int min_reactions, int max_reactions);
/*! \brief Set parameters for the mutationr functions
 * \param double probability of mutating a coefficient (default = 0.5) 
 * \param double probability of remove a reaction during mutation. This may also remove some species.
 * \param double probability of adding a reaction during mutation. This may also add species.
 * \ingroup massaction
*/
void setMutationRatesForMassActionNetwork(double prob_mutate_constants, double prob_mutate_remove_reaction, double prob_mutate_add_reaction);
/*! \brief Set crossover probability
 * \param double probability of crossover (default = 1.0, i.e. always)
 * \ingroup massaction
*/
void setCrossoverRateForMassActionNetwork(double crossover_prob);
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
MassActionNetwork * newMassActionNetwork(int species,int reactions);

/*!\}*/

#endif

