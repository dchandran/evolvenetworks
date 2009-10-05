/*!
  \file    enzymeNetwork.h
  \author: Deepak Chandran (dchandran1@gmail.com)
  \brief   evolve networks with enzymatic reactions

Copyright (C) 2009 Deepak Chandran

	This file defines a chemical reaction network that is very similar to mass-action network ,except
	that it contains an additional enzyme for each reaction with a single reactant and product. 
	For such reactions, the rate expression used is a reversible Hill equation 	(from 
	"The reversible Hill equation: how to incorporate cooperative enzymes into metabolic models" by
	Jan-Hendrik S. Hofmeyr and Athel Cornish-Bowden)

	The struct defined in this file build on the MassActionNetwork struct and functions.

	Example reaction:
		A --> B;  rate = kcat * E * (A - B/A/Keq) * (A/A_half + B/B_half) / (1 + (A/A_half + B/B_half))
		
*********************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "massActionNetwork.h" //builds on mass action network

#ifndef ENZYME_CATALYZED_AND_MASS_ACTION_NETWORK_FOR_GA
#define ENZYME_CATALYZED_AND_MASS_ACTION_NETWORK_FOR_GA

/*! \brief 
Enzyme catalyzed network is defined using a MassActionNetwork struct
and additional arrays for specifying enzyme kinetics. The default rate expression
and stoichiometry is the same as that for MassActionNetwork, except for
reactions with a single reactant and product. In those cases, the rate expression 
used is a reversible Hill equation (from "The reversible Hill equation: how to
incorporate cooperative enzymes into metabolic models" by
Jan-Hendrik S. Hofmeyr and Athel Cornish-Bowden)
\ingroup modifiedmassaction
*/
typedef struct
{
	MassActionNetwork * massActionNetwork;
	/*! \brief enzymes (index values) for each uni-uni reaction*/
	int * enzymes;
	/*! \brief enzyme[i] is an activator if alpha[i] > 1 or an inhibitor if alpha[i] < 1*/
	double * alpha;
	/*! \brief equilibrium ratio of product to substrate*/
	double * Keq;
	/*! \brief half saturation point for substrate*/
	double * S_half;
	/*! \brief half saturation point for product*/
	double * P_half;
	/*! \brief number of modification sites*/
	double * h;
}
EnzymeNetwork;

/*!
   \name  Functions needed by GA
   \{
*/

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

/*! \}
    \name Functions for simulating and printing the network defined above
	\{
*/

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

/*! \}
   \name Functions for initializing a GA
   \{
*/
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
/*! \brief Set the min and max parameters for generating and mutating networks.
 * \param double min rate constant (Vmax for enzymatic reactions)
 * \param double max rate constant (Vmax for enzymatic reactions)
 * \param double min equilibrium point (Keq) in log scale
 * \param double max equilibrium point (Keq) in log scale
 * \param double min number of modification sites for an enzyme (h) in log scale
 * \param double max number of modification sites for an enzyme (h) in log scale
 * \param double min half-saturation point for substrate
 * \param double max half-saturation point for substrate
 * \param double min half-saturation point for product
 * \param double max half-saturation point for product
 * \ingroup modifiedmassaction
*/
void setRateConstantsForEnzymeNetwork(double min_kcat, double max_kcat, 
									  double min_log_keq, double max_log_keq, 
									  double min_log_alpha, double max_log_alpha, 
									  double min_h, double max_h, 
									  double min_s_half, double max_s_half, 
									  double min_p_half, double max_p_half);

/*! \brief Set the initial (average) parameters for generating random networks.
 * \param int min number of species
 * \param int min number of species
 * \param int average number of reactions
 * \ingroup modifiedmassaction
*/
void setSizeForEnzymeNetwork(int, int,int, int);
/*! \brief Set parameters for the mutation and crossover functions. Arguments must add to 1.
 * \param double probability of mutating an enzyme
 * \param double probability of mutating a rate constant (vmax for enzyme reaction)
 * \param double probability of mutating an equilibrium ratio
 * \param double probability of mutating the activation and inhibition profile
 * \param double probability of mutating the number of modification sites on an enzyme
 * \param double probability of mutating the substrate half saturation point
 * \param double probability of mutating the product half saturation point
 * \param double probability of remove a reaction during mutation. This may also remove some species.
 * \param double probability of adding a reaction during mutation. This may also add species.
 * \ingroup modifiedmassaction
*/
void setMutationRatesForEnzymeNetwork(double enzyme, 
									  double k_cat, 
									  double k_eq, 
									  double alpha,
									  double h, 
									  double s_half, 
									  double p_half, 
									  double remove, 
									  double add);
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

/*!\}*/
