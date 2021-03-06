/*!
  \file    proteinNetwork.h
  \author: Deepak Chandran (dchandran1@gmail.com)
  \brief    evolve gene regulatory networks

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
	int * fixed;  //array of size=species. 1= ith species is fixed
} 
ProteinInteractionNetwork;

/*!
   \name  Functions needed by GA
\{*/

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

/*!\}
   \name Functions for simulating and printing the network defined above
\{*/

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
 * \param FILE* output stream
 * \param GAindividual network
 * \ingroup proteinnetwork
*/
void printProteinInteractionNetwork(FILE *, GAindividual);
/*! \brief get the number of variables in the network.
    This is equal to the number of rows in the stoichiometry matrix
 \param ProteinInteractionNetwork network
 \ingroup proteinnetwork
*/
int getNumSpeciesForProteinInteractionNetwork(GAindividual);
/*! \brief get the number of reactions in the network. 
    This is equal to the number of columns in the stoichiometry matrix
 \param ProteinInteractionNetwork network
 \ingroup proteinnetwork
*/
int getNumReactionsForProteinInteractionNetwork(GAindividual);
/*! \brief set a species as a fixed (constant, boundary) species
 \param ProteinInteractionNetwork network
 \param int index of species that should be set at fixed
 \param int value = 0 or 1, where 1 = fixed
 \ingroup proteinnetwork
*/
void setFixedSpeciesForProteinInteractionNetwork(GAindividual, int, int);

/*!\}
   @name Functions for initializing a GA
\{*/
/*! \brief Set the min and max parameters for generating and mutating networks.
 * \param double the min Ka value for random networks
 * \param double the max Ka value for random networks
 * \param double the min Vmax value for random networks
 * \param double the max Vmax value for random networks
 * \param double the min Total (conservation law) for random networks
 * \param double the max Total (conservation law) for random networks
 * \ingroup proteinnetwork
*/
void setRateConstantsForProteinInteractionNetwork(double min_ka, double max_ka, double min_vmax, double max_vmax, double min_total, double max_total);
/*! \brief Set the initial (average) parameters for generating random networks.
 * \param int min number of species in the networks
 * \param int max number of species in the networks
 * \param int min regulations per species
 * \param int max regulations per species
 * \ingroup proteinnetwork
*/
void setSizeForProteinInteractionNetwork(int, int, int, int);
/*! \brief Set parameters for the mutation functions. Arguments must add to 1.
 * \param double probability of rewiring the network during mutation
 * \param double probability of changing network parameter during mutation
 * \param double probability of changing a conservation law during mutation
 * \param double probability of adding/removing a node during mutation
 * \ingroup proteinnetwork
*/
void setMutationRatesForProteinInteractionNetwork(double rewire, double parameter, double total, double addremove);
/*! \brief Set crossover probability.
 * \param double probability of crossover
 * \ingroup proteinnetwork
*/
void setCrossoverRateForProteinInteractionNetwork(double);
/*! \brief Creates an array of randomized protein networks.
 * \return GAindividual* GApopulation of random networks
 * \ingroup proteinnetwork
*/
GApopulation randomProteinInteractionNetworks(int);

/*!\}*/
#endif

