/********************************************************************************************************

Copyright (C) 2009 Deepak Chandran

	This file defines a chemical reaction network that is very similar to mass-action network ,except
	that it contains an additional enzyme for each reaction with a single reactant and product. 
	For such reactions, the rate expression used is Michaeilis-Menten rather than mass-action.

	The struct defined in this file build on the MassActionNetwork struct and functions.

	Example reaction:
		A --> B;  rate = vmax1 * A * E / (Km1 + A)
		
*********************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "massActionNetwork.h" //builds on mass action network

#ifndef ENZYME_CATALYZED_AND_MASS_ACTION_NETWORK_FOR_GA
#define ENZYME_CATALYZED_AND_MASS_ACTION_NETWORK_FOR_GA

/*! \brief 
Enzyme catalyzed network is defined using a MassActionNetwork struct
and an array of enzymes. The enzymes array stores a set of indices
corresponding to the enzyme for that reaction. The default rate expression
and stoichiometry is the same as that for MassActionNetwork. However, for
reaction with a single reactant, the rate expression becomes Michaelis-Menten, 
i.e. hyperbolic in the reactant. All the functions build on the existing
functions of MassActionNetwork. 
	\ingroup modifiedmassaction
*/
typedef struct
{
	MassActionNetwork * massActionNetwork;
	int * enzymes;
	double * Km;
}
EnzymeNetwork;

/*****************************************************
   @name  Functions needed by GA
******************************************************/

/*! \brief Free an individual from memory
 * \param GAindividual a single individual
 * \ingroup modifiedmassaction
*/
void deleteEnzymeNetwork(GAindividual individual);

/*! \brief Make a copy of an individual and return the memory pointer
 * \param GAindividual target individual
 * \ingroup modifiedmassaction
*/
GAindividual cloneEnzymeNetwork(GAindividual individual);

/*! \brief Combine two individuals to generate a new individual
 * \param GAindividual parent individual 1
 * \param GAindividual parent individual 2
 * \return GAindividual pointer to an individual (can be the same as one of the parents)
 * \ingroup modifiedmassaction
*/
GAindividual crossoverEnzymeNetwork(GAindividual individualA, GAindividual individualB);

/*! \brief Change an individual randomly to generate a new individual
 * \param GAindividual parent individual
 * \return GAindividual pointer to an individual (can be the same as one of the parents)
 * \ingroup modifiedmassaction
*/
GAindividual mutateEnzymeNetwork(GAindividual individual);

/*****************************************************
   @name Functions for simulating and printing the network defined above
******************************************************/

/*! \brief Propensity function to be used by the SSA function (see ssa.h)
 * \param double time
 * \param double* values for variables in the system
 * \param double* rates (output)
 * \param GAindividual network
 * \ingroup modifiedmassaction
*/
void ratesForEnzymeNetwork(double,double*,double*,GAindividual);
/*! \brief Get the stoichiometry matrix for a network
 * \param GAindividual network
 * \return double* linearized stoichiometry matrix
 * \ingroup modifiedmassaction
*/
double * stoichiometryForEnzymeNetwork(GAindividual);
/*! \brief Print a network
 * \param FILE* output stream
 * \param GAindividual network
 * \ingroup modifiedmassaction
*/
void printEnzymeNetwork(FILE *, GAindividual);
/*! \brief get the number of variables in the network.
    This is equal to the number of rows in the stoichiometry matrix
 \param EnzymeNetwork network
 \ingroup modifiedmassaction
*/
int getNumSpeciesForEnzymeNetwork(GAindividual);
/*! \brief get the number of reactions in the network. 
    This is equal to the number of columns in the stoichiometry matrix
 \param EnzymeNetwork network
 \ingroup modifiedmassaction
*/
int getNumReactionsForEnzymeNetwork(GAindividual);
/*! \brief set a species as a fixed (constant, boundary) species
 \param EnzymeNetwork network
 \param int index of species that should be set at fixed
 \param int value = 0 or 1, where 1 = fixed
 \ingroup modifiedmassaction
*/
void setFixedSpeciesForEnzymeNetwork(GAindividual, int, int);

/*****************************************************
   @name Functions for initializing a GA
******************************************************/
/*! \brief Set the initial (average) parameters for generating random networks.
 * \param double percent of initial reactions having one reactant and one product
 * \param double percent of initial reactions having one reactant and two products
 * \param double percent of initial reactions having two reactants and one product
 * \param double percent of initial reactions having two reactant and two product
 * \param double percent of initial reactions representing creation (i.e. no reactants)
 * \param double percent of initial reactions representing degradation (i.e. no products)
 * \ingroup modifiedmassaction
*/
void setDistributionOfEnzymeNetwork(double uni_uni, double uni_bi, double bi_uni, double bi_bi, double no_reactant, double no_product);
/*! \brief Set the initial (average) parameters for generating random networks.
 * \param double average rate constant (Vmax for enzymatic reactions)
 * \param double average half-saturation point for enzyme reactions (Km)
 * \ingroup modifiedmassaction
*/
void setRateConstantsForEnzymeNetwork(double avg_rate_constant, double avg_km);

/*! \brief Set the initial (average) parameters for generating random networks.
 * \param int average number of species
 * \param int average number of reactions
 * \ingroup modifiedmassaction
*/
void setSizeForEnzymeNetwork(int, int);
/*! \brief Set parameters for the mutation and crossover functions. Arguments must add to 1.
 * \param double probability of mutating an enzyme (default = 0.2) 
 * \param double probability of mutating a rate constant (default = 0.5) 
 * \param double probability of remove a reaction during mutation. This may also remove some species.
 * \param double probability of adding a reaction during mutation. This may also add species.
 * \ingroup modifiedmassaction
*/
void setMutationRatesForEnzymeNetwork(double, double, double, double);
/*! \brief Set crossover probabilty
 * \param double probability of crossover (default = 1.0, i.e. always)
 * \ingroup modifiedmassaction
*/
void setCrossoverRateForEnzymeNetwork(double);
/*! \brief Creates an array of randomized enzyme networks.
 * \return GAindividual* GApopulation of random networks
 * \ingroup modifiedmassaction
*/
GApopulation randomEnzymeNetworks(int);

/*! \brief Make a new empty network
 * \param int number of species
 * \param int number of reactions
 * \return EnzymeNetwork* network
 * \ingroup modifiedmassaction
*/
EnzymeNetwork * newEnzymeNetwork(int,int);

#endif

