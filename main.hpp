/**
 * @mainpage Genetic Algorithm and Simulation Pack
 * @brief 

stochastic simulation -- user needs to specify propensity and stoichiometry matrix

deterministic simulation -- user needs to specify function with differential equations

multi-cell stochastic simulation -- uses logistic equations to model cell growth, with a network inside each cell. 
user needs to specify propensity and stoichiometry matrix

genetic algorithm -- users needs to specify fitness, mutation, crossover functions

mass action genetic algorithm -- evolves a mass-action network; user provides fitness function

protein network genetic algorithm -- evolves a protein-network (enzyme kinetics); user provides fitness function

gene regulatory network genetic algorithm -- evolves a genetic network; user provides fitness function

reaction network genetic algorithm -- evolves any of the above network types; user provides fitness function

find bistability -- finds parameters to make a system bistable; user provided differential equation function

negative and positive loops -- find negative and positive loops in a square adjacency matrix such as a Jacobian
 
*/

/*! \defgroup gillespie Stochastic simulation of single cell and population of cells
   \brief
    Defines standard Gillespie stochastic simulation algorithm, which accepts 
	propensity function and stoichiometry matrix from user
	
	Defines a simulator that performs a multi-cell stochastic simulation using 
	a slight variant of the Gillespie algorithm. 
	
	The algorithm has two levels of events (or reactions):
	
	Level 1: Cell birth and death
	Level 2: reactions inside each cell
	
	An overview of the algorithm:
	
		current time: t0
		Calculate next time for cell birth or cell death: t1
	    Calculate probability of cell death and event
		
		select an event: cell birth or cell death
		
		update cell count
		
		from t0 until t1:
		    for each cell:
			    initialize concentrations for that cell
				do Gillespie simulation
				store final concentrations for that cell
				
	
	Output: all values are aligned to a grid. 

 */

/*! \defgroup cvodewrapper A wrapper over the Sundials CVODE numerical integrator
   \brief 

A wrapper for the Sundials CVODE numerical integration library
The wrapper is only for initial value problems. 
Each wrapper function has two versions: 
	one that accepts a function with differential equation definitions
	one that accepts a stoichiometry matrix and a propensity function
The wrapper provides support for event functions and user defined structs as parameter.

 */

/*! \defgroup ga Generic genetic algorithm
   \brief 
   	This library provides the basic functions for running a Genetic Algorithm (GA), but
	the user is required to setup the fitness function and related functions.
	
	The user MUST define:
		1) a struct that represents an "individual"
		2) a function that returns the fitness of an individual
		3) a mutation function that randomly alters an individual
		4) a function to free an individual
		5) a function to clone an individual
	
	The following function definition is optional but highly recommended:
		1) a crossover function to make a new individual from two individuals
	
	The following functions definitions are entirely optional:
		1) A function that selects individuals using fitness as probabilities (library provides one by default)
		2) A callback function can be used to examine or terminate the GA at any iteration
		
	The main functions are: GAinit and GArun
 */

/*! \defgroup geneticnetwork Evolve bi-molecular network
   \brief 
   	This file defines a chemical reaction network using fractional saturation rate laws. 
	The 
	
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
 
 /*! \defgroup proteinnetwork Evolve protein interaction network
   \brief
   	This file defines a chemical reaction network where each protein can switch from active to inactive form
	via Michaelis Menten type kinetics. The number of enzymes affecting each transition can be many. 
	The functions in this file are designed to be used with the GA library that I have written. 
	
	This file defines the following functions that are required in the GA:
		1) a struct that represents an "individual"
		2) a mutation function to randomly alter an individual
		3) a crossover function to make a new individual from two individuals
	
	In addition, the file defines the following functions for simulating a network:
		1) stoichiometry function 
		2) rates function
	
	The program using this file must define:
		1) a function that returns the fitness of a network
 */

/*! \defgroup massaction Evolve bi-molecular network
   \brief 
   	This file defines a chemical reaction network using mass-action kinetics. The functions in this file
	are designed to be used with the GA library that I have written. 
	
	This file defines the following functions that are required in the GA:
		1) a struct that represents an "individual"
		2) a mutation function to randomly alter an individual
		3) a crossover function to make a new individual from two individuals
	
	In addition, the file defines the following functions for simulating a network:
		1) stoichiometry function 
		2) rates function
	
	The program using this file must define:
		1) a function that returns the fitness of a network
 */
 
 /*! \defgroup genericNetwork Evolve any of the other networks tyes
   \brief 
   	This file defines a generic reaction network using either of the other network architectures. 
	The functions in this file are designed to be used with the GA library that I have written. 
	
	This file defines the following functions that are required in the GA:
		1) a struct that represents an "individual"
		2) numerous functions that utilize the existing networks

	The program using this file must define:
		1) a function that returns the fitness of a network

	Any of the functions in the library can be replaced using the setXXX functions
 */

/*! \defgroup bistable Find parameters for bistability
	\brief
	This file uses genetic algorithm optimization method to find parameters for a system of ODE such 
	that the ODE system will have two or more stable states. 
	
	The objective function tries to find a saddle or unstable point in the ODE system by doing
	reverse time or partial-reverse time simulation. Once a saddle or unstable point is found, the algorithm
	solves the IVP problem starting around the unstable or saddle point, hoping to reach two stable points. 
	
	The unstable or saddle points are found by optimizing: 
	
		p (parameters of the system) and a (vector of reals)
		
	where the ODE system is:
	
		dx/dt = a * f(x,p)
		
	where f(x,p) is the original system defined by the user.
	The alpha vector alters the stability of the critical points, allowing the algorithm to search
	for saddle points of the original system by searching for stable points of this new transformed system.
*/

/*! \defgroup loops Find positive and negative loops in adjacency matrix (Jacobian)
	\brief
	This file contains an algorithm for finding loops in an adjacency matrix. 
	The edges can be negative or positive, such as in a Jacobian matrix. 
	
	The final output will list negative and positive loops. A negative loop is where the total
	product of all the edges in the loop is negative. Positive loop is a loop that is not negative.
	
	Originally, this code was used to find loops in the Jacobian matrix of an ODE. Loops
	in the Jacobian can indicate locations of positive feedback or negative feedback. However, 
	the code is genetic for any square matrix. 
*/
