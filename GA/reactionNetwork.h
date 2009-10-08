/*!
  \file    reactionNetwork.h
  \author: Deepak Chandran (dchandran1@gmail.com)
  \brief    evolve any of the other biological networks

	Copyright (C) 2009 Deepak Chandran
	
	This file provides a convenient wrapper for the funtions in:
	
	massActionNetwork.h
	proteinInteractionNetwork.h
	enzymeNetwork.h
	geneRegulationNetwork.h
	
	Each file above contains the function needed for running 
	evolving networks using genetic algorithm (ga.h). This
	file defines a single network structure called "ReactionNetwork"
	that can contain any of the three above network types. 
	
	The mutation, crossover, simulation, and all the other functions
	from the three networks above are used to perform the same
	functions for this combined network.

********************************************************/

#ifndef GENERIC_REACTION_NETWORK_GA_WRAPPER_H
#define GENERIC_REACTION_NETWORK_GA_WRAPPER_H

#include "ga.h"
#include "massActionNetwork.h"
#include "enzymeNetwork.h"
#include "proteinInteractionNetwork.h"
#include "geneRegulationNetwork.h"

#define MASS_ACTION_NETWORK 0
#define ENZYME_NETWORK 1
#define PROTEIN_INTERACTION_NETWORK 3
#define GENE_REGULATION_NETWORK 2


/*! \brief A generic network that can contain either a 
	mass-action network, enzyme network, protein-interaction network, or 
	gene-regulation network. The type information is used
	to call the corresponding functions that are already 
	defined in each of the other network types.
 \ingroup genericNetwork
*/

typedef struct 
{
	/*! \brief MASS_ACTION_NETWORK, PROTEIN_INTERACTION_NETWORK, or GENE_REGULATION_NETWORK*/
	int type; 
	
	/*! \brief the network itself*/
	GAindividual network;

	/*! \brief network ID*/
	int id;

	/*! \brief parent network IDs for following evolutionary lineage (null terminated)*/
	int * parents;

	/*! \brief initial values*/
	double * initialValues;
}
ReactionNetwork;

/*!
  \name Get network information
  \{
*/

/*! \brief print the reaction network to stdout.
    This function calls the print functions of one 
	of the individual networks.
 \param FILE* where to print (e.g. stdout)
 \param ReactionNetwork network
 \ingroup genericNetwork
*/
void printNetwork(FILE * stream, GAindividual);
/*! \brief print the reaction network to a file. 
    This function calls the print functions of one 
	of the individual networks.
 \param ReactionNetwork network
 \param char* file name
 \ingroup genericNetwork
*/
void printNetworkToFile( char * filename, GAindividual);
/*! \brief get the number of variables (e.g. genes, molecular species, etc.) in the network.
    This is equal to the number of rows in the stoichiometry matrix
 \param ReactionNetwork network
 \ingroup genericNetwork
*/
int getNumSpecies(GAindividual);
/*! \brief get the number of reactions in the network. 
    This is equal to the number of columns in the stoichiometry matrix
 \param ReactionNetwork network
 \ingroup genericNetwork
*/
int getNumReactions(GAindividual);
/*! \brief get the stoichiometry matrix for the network
 \param ReactionNetwork network
 \return double* linearized 2D matrix -- use getValue(i,j)
 \ingroup genericNetwork
*/
double* getStoichiometryMatrix(GAindividual);
/*! \brief get the reaction rates at a given point
 \param ReactionNetwork network
 \param double* concentration values at which the rate should be calculated
 \return double* values for the rates at this point
 \ingroup genericNetwork
*/
double* getReactionRates(GAindividual, double*);

/*!
  \}
  @name Related to lineage tracking
  \{
*/
/*! \brief turn on lineage tracking. This will track the parents of each individual when crossover occurs
 \ingroup genericNetwork
*/
void lineageTrackingON();

/*! \brief turn on lineage tracking. 
	This will prevent tracking of the parents of each individual when crossover occurs
	ReactioNetwork's parent field will be 0.
 \ingroup genericNetwork
*/
void lineageTrackingOFF();

/*! \brief set the ID of this individual
 \param ReactionNetwork network
 \param int ID for this individual
 \ingroup genericNetwork
*/
void setID(GAindividual,int);

/*! \brief get the ID of this individual
 \param ReactionNetwork network
 \return int ID for this individual
 \ingroup genericNetwork
*/
int getID(GAindividual);

/*! \brief get the ID for all ancestors of this individual
 \param ReactionNetwork network
 \return int* NULL TERMINATED array with IDs of parents. DO NOT FREE
 \ingroup genericNetwork
*/
int* getParentIDs(GAindividual);

/*!
  \}
  @name Simulation functions
  \{
*/

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

/*! \brief set the initial values of the variables in the network
 \param ReactionNetwork network
 \param double* initial concentrations
 \ingroup genericNetwork
*/
void setInitialValues( GAindividual, double *);

/*! \brief get the initial values of the variables in the network
 \param ReactionNetwork network
 \return double* initial concentrations
 \ingroup genericNetwork
*/
double * getInitialValues( GAindividual );

/*! \brief simulate using a system of ordinary differential equations (ODE) (see cvodesim.h)
 \param ReactionNetwork network to simulate
 \param double total simulation time
 \param double step-size for simulation, 0-1
 \return double* linearized 2D array, use getValue (see cvodesim.h)
 \ingroup genericNetwork
*/
double * simulateNetworkODE( GAindividual, double, double );

/*! \brief get steady state of a system using a system (see cvodesim.h)
 \param ReactionNetwork network to simulate
 \return double* steady state values
 \ingroup genericNetwork
*/
double * networkSteadyState( GAindividual );

/*! \brief simulate using stochastic simulation algorithm (see ssa.h)
 \param ReactionNetwork network to simulate
 \param double total simulation time
 \param int* returns the final size (rows) of the matrix here
 \param double step-size for simulation, 0-1
 \return double* linearized 2D array, use getValue (see cvodesim.h)
 \ingroup genericNetwork
*/
double * simulateNetworkStochastically( GAindividual, double, int* );

/*!
  \}
  @name Functions for GA
  \{
*/

/*! \brief set the fitness function for the GA. Use before calling evolveNetworks.
	IMPORTANT: higher fitness = better. If you want to minimize a function, simple invert
	the value, e.g. 1/(1+x) or something similar. Otherwise, you will need to override the selection function using
	GAsetSelectionFunction, because the default selection function assumes higher fitness is better. 
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

/*! \brief set the min and max network sizes. Use before calling evolveNetworks.
     The "size" can mean different things for different types of networks. In general, it corresponds to the number
	 of rows and columns in the stoichiometry matrix. In somecases, the columns might be twice the given size, due
	 to degradation reactions. 
 \param int min number of variables (e.g, species or genes)
 \param int max number of variables (e.g, species or genes)
 \param int min number of reactions (e.g, enzymatic or transcriptional)
 \param int max number of reactions (e.g, enzymatic or transcriptional)
 \ingroup genericNetwork
*/
void setNetworkSize(int min_vars,int max_vars,int min_reactions,int max_reactions);

/*! \brief Run the genetic algorithm. 
     Calls GArun from ga.h using a random population created using randomNetworks()
 \param int number of individuals in the initial population (use large number here)
 \param int number of individuals in each successive population (use relatively small number for speed)
 \param int number of generations for evolution
 \param GACallbackFnc a callback function (optional, use 0 for none)
 \return GApopulation the final evolved population of networks. Population is sorted by fitness of individuals.
		The first individual in the population (index 0) will be the best individual.
 \ingroup genericNetwork
*/
GApopulation evolveNetworks(int init_popSz,int final_popSz,int iterations,GACallbackFnc callback);

/*!
  \}
  @name functions for creating and removing networks
  \{
*/

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

/*!
  \}
  @name functions related to mutation and crossover
  \{
*/

/*! \brief set the probability of crossover (default = 1.0, i.e. always)
	\param double value between 0 and 1
	\ingroup genericNetwork
*/
void setCrossoverRate(double);

/*! \brief set the average initial value for the variables in a network. 
		The initial values are used to initialize a simulation.
	\param double positive real
	\ingroup genericNetwork
*/
void setAverageInitialValue(double);

/*! \brief set the probability of mutating the initial values
	\param double value between 0 - 1 ( default = 0.2 )
	\ingroup genericNetwork
*/
void setMutationRateOfInitialValues(double);

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

/*!
  \}
  \name Premade fitness functions
  \{
*/

/*! \brief compares the steady state values of the network to an input/output table provided. 
           Sets the first n input species as fixed species, and compares the output
		   again each of the other species.
		   returns the fitness value 1/(1 + sum_of_square_differences), where
		   sum_of_square_differences is the sum of square differences between the desired
		   outputs and the best matching species in the network.
	\param GAindividual must be a ReactionNetwork* casted as GAindividual
	\param double** input/output table. 
			Number of inputs + number of outputs <= number of species in network.
	\param int number of rows in the input table
	\param int number of input columns (first set of columns)
	\param int number of output columns (last set of columns)
	\param int 1=use correlation, not absolute differences. 0 = use absolute differences.
	\param double ** can be 0. if non-zero, this matrix MUST BE THE SAME SIZE as the input table. 
			The output from the network will be placed in this table. This is for the purpose of
			comparing the results against the original table.
	\return double fitness score = 1/(1 + sum_of_square_differences)
	\ingroup genericNetwork
*/
double compareSteadyStates(GAindividual, double **, int , int, int, int, double ** );

/*!
  \}
  \name Keeping a log file
  \{
*/


/*! \brief enable automatic log file. The log file will store information about every generation
	as well as the final result of the GA. Setting the file will enable automatic log
	\param char* log file name
*/
void enableLogFile(char * filename);

/*! \brief disable automatic log file. */
void disableLogFile();

/*! \brief All arguments must be 0 or 1 (Boolean). This function allows the evolution experiment
	to report information during each generation. The information to report can be configured
	using the arguments. The GA uses a custom callback routine to write the information to 
	the log file.
	\param int report the best fitness score during each generation
	\param int print the best network's script during each generation
	\param int report the best network's size during each generation
	\param int report the best network's parents during each generation
	\param int report all fitness values during each generation
	\param int report all networks' parents during each generation
*/
void configureContinuousLog(int bestNetworkFitness, 
							int bestNetworkScript,
							int bestNetworkSize, 
							int bestNetworkLineage,
							int allFitness,
							int allNetworkLineage );

/*! \brief All arguments must be 0 or 1 (Boolean). This function allows the evolution experiment
	to report information at the end of the GA.
	\param int report the best fitness score at the end
	\param int print the best network's script at the end
	\param int report the best network's size at the end
	\param int report the best network's parents at the end
	\param int report all fitness value at the end
	\param int report all networks' parents at the end
	\param int report the random number generator's seeds (useful for duplicating the experiment)
*/
void configureFinalLog(int bestNetworkFitness, 
							int bestNetworkScript,
							int bestNetworkSize, 
							int bestNetworkLineage,
							int allFitness, 
							int allNetworkLineage,
							int seeds);

/*! \brief set parameters for networkSteadyState()
	\param double the allowed error. When the sum of squares of all derivatives is
			between two time points that is delta apart is within this error range, 
			then the system is considered to be at steady state.
			default = 1.0E-3 (quite good for the default delta)
	\param double the time separation for checking error tolerance.
			default = 0.1
	\param double maximum time allows. after this time limit, a 0 is returned, i.e. no steady state
			default = 10000.0 (you may want to lower this if you want speed).
*/
void configureSteadyStateFunction(double tolerance, 
									double delta,
									double maxTime);


/*!\}*/

#endif
