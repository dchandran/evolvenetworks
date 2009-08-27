/********************************************************************************************************

Copyright (C) 2009 Deepak Chandran

	This file defines a gene regulatory network using fractional saturation models for the kinetics. 
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

#ifndef GENE_REGULATORY_FOR_GA
#define GENE_REGULATORY_FOR_GA

/*! \brief
  A complex composed of one or more transcription factors
*/
typedef struct
{
	int * TFs;
	int size;
}
TFComplex;

/*! \brief Gene regulatory network
	\ingroup geneticnetwork
*/
typedef struct
{
	int species;  //number of genes
	int numComplexes; //size of each array defined in this struct
	
	TFComplex * complexes; //a list of complexes
	int * targetGene; //the gene that is regulated by the complexes
	double * Ka; //the association constant (if negative, then it is a repressor)
	
	double * Vmax; //the Vmax for each protein
	double * degradation; //degradation rate for each protein
	
	int * fixed;  //array of size=species. 1= ith species is fixed
} 
GeneRegulationNetwork;

/*****************************************************
   @name  Functions needed by GA
******************************************************/

/*! \brief Free an individual from memory
 * \param GAindividual a single individual
 * \ingroup geneticnetwork
*/
void deleteGeneRegulationNetwork(GAindividual individual);

/*! \brief
 * Make a copy of an individual and return the memory pointer
 * \param GAindividual target individual
 * \ingroup geneticnetwork
*/
GAindividual cloneGeneRegulationNetwork(GAindividual individual);

/*! \brief
 * combine two individuals to generate a new individual
 * \param GAindividual parent individual 1
 * \param GAindividual parent individual 2
 * \return pointer to an individual (can be the same as one of the parents)
 * \ingroup geneticnetwork
*/
GAindividual crossoverGeneRegulationNetwork(GAindividual individualA, GAindividual individualB);

/*! \brief
 * Change an individual randomly to generate a new individual
 * \param GAindividual parent individual
 * \return pointer to an individual (can be the same as one of the parents)
 * \ingroup geneticnetwork
*/
GAindividual mutateGeneRegulationNetwork(GAindividual individual);

/*****************************************************
   @name Functions for simulating and printing the network defined above
******************************************************/
/*! \brief
 * Propensity function to be used by the SSA function (see ssa.h)
 * \param double time
 * \param double values for variables (species)
 * \param double rates (output)
 * \param GAindividual network
 * \ingroup geneticnetwork
*/
void ratesForGeneRegulationNetwork(double,double*,double*,GAindividual);
/*! \brief
 * Get the stoichiometry matrix for a network
 * \param GAindividual network
 * \return linearized stoichiometry matrix
 * \ingroup geneticnetwork
*/
double * stoichiometryForGeneRegulationNetwork(GAindividual);
/*! \brief
 * Print a network
 * \param GAindividual network
 * \ingroup geneticnetwork
*/
void printGeneRegulationNetwork(GAindividual);
/*! \brief get the number of variables in the network.
    This is equal to the number of rows in the stoichiometry matrix
 \param GeneRegulationNetwork network
 \ingroup geneticnetwork
*/
int getNumSpeciesForGeneRegulationNetwork(GAindividual);
/*! \brief get the number of reactions in the network. 
    This is equal to the number of columns in the stoichiometry matrix
 \param GeneRegulationNetwork network
 \ingroup geneticnetwork
*/
int getNumReactionsForGeneRegulationNetwork(GAindividual);
/*! \brief set a gene as a fixed (constant, boundary) species
 \param GeneRegulationNetwork network
 \param int index of gene that should be set at fixed
 \param int value = 0 or 1, where 1 = fixed
 \ingroup geneticnetwork
*/
void setFixedSpeciesForGeneRegulationNetwork(GAindividual, int, int);

/*****************************************************
   @name Functions for initializing a GA
******************************************************/
/*! \brief Set the initial parameters for generating random networks.
 * \param int maximum number of transcription factors in a complex (i.e. 2 means that only monomers or dimers allowed)
 * \param double the average Ka value for random networks
 * \param double the average maximum production rate for a gene
 * \param double the average degradation rate constant for a gene product
 * \ingroup geneticnetwork
*/
void setParametersForGeneRegulationNetwork(int, double, double, double);

/*! \brief Set the initial parameters for generating random networks.
 * \param int average number of genes in the networks
 * \param int average number of regulations per gene in the networks
 * \ingroup geneticnetwork
*/
void setSizeForGeneRegulationNetwork(int,int);

/*! \brief Set parameters for the mutation and crossover functions. First four arguments must add to 1.
 * \param double probability of mutating a random Ka value
 * \param double probability of mutating a degradation value or Vmax value
 * \param double probability of mutating a component of a random complex
 * \param double probability of adding a gene
 * \param double probability of removing a gene
 * \param double probability of crossover
 * \ingroup geneticnetwork
*/
void setMutationAndCrossoverRatesForGeneRegulationNetwork(double, double, double, double, double, double);

/*! \brief
 * Creates an array of randomized genetic networks.
 * \param int number of networks in population
 * \return GAindividual* GApopulation of random networks
 * \ingroup geneticnetwork
*/
GApopulation randomGeneRegulationNetworks(int);

/*! \brief
 * Make a new empty network
 * \param number of species
 * \param number of regulations
 * \return network
 * \ingroup geneticnetwork
*/
GeneRegulationNetwork * newGeneRegulationNetwork(int,int);

#endif

