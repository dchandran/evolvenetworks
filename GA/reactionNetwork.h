/*******************************************************

	Copyright (C) 2009 Deepak Chandran
	
	This file provides a convenient wrapper for the funtions in:
	
	massActionNetwork.h
	proteinInteractionNetwork.h
	geneRegulationNetwork.h
	
	Each file above contains the function needed for running 
	evolving networks using genetic algorithm (ga.h). This
	file defines a single network structure called "ReactionNetwork"
	that can contain any of the three above network types. 
	
	The mutation, crossover, simulation, and all the other functions
	from the three networks above are used to perform the same
	functions for this combined network.

********************************************************/

#include "ga.h"
#include "massActionNetwork.h"
#include "proteinInteractionNetwork.h"
#include "geneRegulationNetwork.h"

#define MASS_ACTION_NETWORK 1
#define PROTEIN_INTERACTION_NETWORK 2
#define GENE_REGULATION_NETWORK 4

#define MASS_ACTION_NETWORK_INDEX 0
#define PROTEIN_INTERACTION_NETWORK_INDEX 1
#define GENE_REGULATION_NETWORK_INDEX 2

/*! \brief A generic network that can contain either a 
	mass-action network, protein-interaction network, or 
	gene-regulation network
 \ingroup genericNetwork
*/

typedef struct 
{
	//MASS_ACTION_NETWORK, PROTEIN_INTERACTION_NETWORK, or GENE_REGULATION_NETWORK
	int type; 
	
	//the network itself
	GAindividual network;
}
ReactionNetwork;

/***********************************************************
  @name Get network information
************************************************************/

/*! \brief print the reaction network. 
    This function simply calls the print functions of one 
	of the individual networks.
 \param ReactionNetwork network
 \ingroup genericNetwork
*/
void printNetwork(ReactionNetwork *);
/*! \brief get the number of variables (e.g. genes, molecular species, etc.) in the network.
    This is equal to the number of rows in the stoichiometry matrix
 \param ReactionNetwork network
 \ingroup genericNetwork
*/
int getNumSpecies(ReactionNetwork *);
/*! \brief get the number of reactions in the network. 
    This is equal to the number of columns in the stoichiometry matrix
 \param ReactionNetwork network
 \ingroup genericNetwork
*/
int getNumReactions(ReactionNetwork *);
/*! \brief get the stoichiometry matrix for the network
 \param ReactionNetwork network
 \return double* linearized 2D matrix -- use getValue(i,j)
 \ingroup genericNetwork
*/
double* getStoichiometryMatrix(ReactionNetwork *);
/*! \brief get the reaction rates at a given point
 \param ReactionNetwork network
 \param double* concentration values at which the rate should be calculated
 \return double* values for the rates at this point
 \ingroup genericNetwork
*/
double* getReactionRates(ReactionNetwork *, double*);

/***********************************************************
  @name Simulation functions
************************************************************/

/*! \brief set the rates calculating function for the given type of network
 \param int type of the network
 \param PropensityFunction rates function (see cvodesim.h)
 \ingroup genericNetwork
*/
void setRatesFunction( int, PropensityFunction );

/*! \brief set the stoichiometry calculating function for the given type of network
 \param int type of the network
 \param double*(*f) stoichiometry function (see cvodesim.h)
 \ingroup genericNetwork
*/
void setStoichiometryFunction( int, double* (*f)(GAindividual) );

/*! \brief simulate using a system of ordinary differential equations (ODE) (see cvodesim.h)
 \param ReactionNetwork network to simulate
 \param double* initial concentrations
 \param double total simulation time
 \param double step-size for simulation, 0-1
 \return double* linearized 2D array, use getValue (see cvodesim.h)
 \ingroup genericNetwork
*/
double * simulateNetworkODE( ReactionNetwork *, double*, double, double );

/*! \brief get steady state of a system using a system (see cvodesim.h)
 \param ReactionNetwork network to simulate
 \param double* initial concentrations
 \return double* steady state values
 \ingroup genericNetwork
*/
double * networkSteadyState( ReactionNetwork *, double* );

/*! \brief simulate using stochastic simulation algorithm (see ssa.h)
 \param ReactionNetwork network to simulate
 \param double* initial concentrations
 \param double total simulation time
 \param int* returns the final size (rows) of the matrix here
 \param double step-size for simulation, 0-1
 \return double* linearized 2D array, use getValue (see cvodesim.h)
 \ingroup genericNetwork
*/
double * simulateNetworkStochastically( ReactionNetwork *, double*, double, int* );

/**************************************
  @name Functions for GA
***************************************/

/*! \brief set the fitness function for the GA. Use before calling evolveNetworks.
 * \param GAFitnessFnc function pointer (cannot be 0)
 \ingroup genericNetwork
*/
void setFitnessFunction(GAFitnessFnc);

/*! \brief set the crossover function for one network type. Use before calling evolveNetworks.
 \param int type of the network, e.g. e.g. MASS_ACTION_NETWORK | PROTEIN_INTERACTION_NETWORK
 \param GACrossoverFnc function pointer
 \ingroup genericNetwork
*/
void setCrossoverFunction( int , GACrossoverFnc);

/*! \brief set the mutation function for one network type. Use before calling evolveNetworks.
 \param int type of the network, e.g. e.g. MASS_ACTION_NETWORK | PROTEIN_INTERACTION_NETWORK
 \param GACrossoverFnc function pointer
 \ingroup genericNetwork
*/
void setMutationFunction( int , GAMutateFnc );

/*! \brief set the proportional of the population that one type of networks occupies. Use before calling evolveNetworks.
 \param int type of the network, e.g. e.g. MASS_ACTION_NETWORK | PROTEIN_INTERACTION_NETWORK
 \param double probability, 0-1
 \ingroup genericNetwork
*/
void setNetworkTypeProbability(int, double);

/*! \brief use this function to make the GA use only one type of network. This is the same as setting 
     one network type's probability to 1 and the others to 0 using setNetworkTypeProbability. Use before calling evolveNetworks.
 \param int type of the network, e.g. e.g. MASS_ACTION_NETWORK | PROTEIN_INTERACTION_NETWORK
 \ingroup genericNetwork
*/
void setNetworkType(int);

/*! \brief set the average network size of the initial population. Use before calling evolveNetworks.
     The "size" can mean different things for different types of networks. In general, it corresponds to the number
	 of rows and columns in the stoichiometry matrix.
 \param int number of variables (e.g, species or genes)
 \param int number of reactions (e.g, enzymatic or transcriptional)
 \ingroup genericNetwork
*/
void setInitialNetworkSize(int,int);

/*! \brief print the reaction network. 
    This function simply calls the print functions of one 
	of the individual networks.
 \param int number of individuals in the initial population (use large number here)
 \param int number of individuals in each successive population (use relatively small number for speed)
 \param int number of generations for evolution
 \param GACallbackFnc a callback function (optional, use 0 for none)
 \return GApopulation the final evolved population of networks
 \ingroup genericNetwork
*/
GApopulation evolveNetworks(int,int,int,GACallbackFnc);

/***********************************************************
  @name functions for allocating and deallocating networks
************************************************************/

/*! \brief get a random set of networks (used inside evolveNetworks)
	\param int number of networks
	\return GApopulation
	\ingroup genericNetwork
*/
GApopulation randomNetworks(int);

/*! \brief deallocate the space occupied by a network
	\param GAindividual must be a ReactionNetwork* casted as GAindividual
	\ingroup genericNetwork
*/
void deleteNetwork(GAindividual);

/*! \brief get a copy of a network
	\param GAindividual must be a ReactionNetwork* casted as GAindividual
	\return GAindividual clone
	\ingroup genericNetwork
*/
GAindividual cloneNetwork(GAindividual);

/*! \brief mutate the network
	\param GAindividual must be a ReactionNetwork* casted as GAindividual
	\return GAindividual mutant of the input network
	\ingroup genericNetwork
*/

GAindividual mutateNetwork(GAindividual);

/*! \brief crossover two networks if they are of the same type. Otherwise, same as mutateNetwork
	\param GAindividual must be a ReactionNetwork* casted as GAindividual
	\param GAindividual must be a ReactionNetwork* casted as GAindividual
	\return GAindividual combination of the two input networks if same network type, 
	        otherwise just a mutant of the input network
	\ingroup genericNetwork
*/
GAindividual crossoverNetwork(GAindividual, GAindividual);
