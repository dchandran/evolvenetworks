/**
 * @mainpage Genetic Algorithm and Simulation Pack
 * @brief 

	The entire library can be divided into two main sections and a third set containing several different algorithms. 
	
	set 1: simulation algorithms (mostly wrappers of existing libraries such as lapack and sundials)
	
	set 2: genetic algorithm that can be used to evolve modular biological networks			
		example use cases include: evolving oscillatory networks, evolving logic gates, and evolving noise-supressing networks
	
	set 3: several small algorithms such as:
			1. optimizing a system's parameters for bistability
			2. evolution of chemotaxis network
			3. information theoretic approach for finding important parameters of a system
 
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

 /*! \defgroup gablocks Evolving modular biological systems
   \brief 
   This section of the library defines the structs called Block and System. A Block represents a single biological process
   that is defined by a stoichiometry matrix and reaction rate equations. Each block contains unique
   parameter values, inputs, outputs, and internal molecules. The input and output molecules are used to connect
   two or more blocks to each other. Each block has a BlockType, which is a struct that stores the stoichiometry and 
   rate function for blocks of that type. It also stores the number of inputs and outputs and other information
   about all blocks of that type. 
   
   A System is composed of a set number of molecules and a set number of interacting blocks. By generating
   the rates and stoichiometry matrices for each block, the entire system can be reduced to a single
   dynamical system.
   
   The types of Blocks, or modules, can be extended by modifying functions.h   
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
