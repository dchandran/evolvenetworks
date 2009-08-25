Network Evolution Library (2009)
author: Deepak Chandran (dchandran1@gmail.com)
LICENSE: open source BSD

--------------------------------------------------------
--------------------------------------------------------
How to Build the project
--------------------------------------------------------
--------------------------------------------------------

1. install CMake. 
2. Run the CMake GUI
3. select the network evolution directory as the Source and Build directory
4. "Configure"
5. Select the C compiler that you would like to use *
6. May need to "Configure" again
6. "Generate"

* 
Note: the CMake project has been tested using gcc in Unix and MinGW in Windows. Visual Studio "should" work if it is configured for C code

-----------------------
Without using CMake:
-----------------------

If you want to build one of the example C programs in the examples/ folder, you will need to compile all the files in the following folders:
external
simulation
GA

You may use any standard C compiler and compile all the files in the above folders along with one of the C files in the examples folder. You can then replace the example code with your own code.



--------------------------------------------------------
--------------------------------------------------------
What libraries get generated
--------------------------------------------------------
--------------------------------------------------------

ode -- used for deterministic simulation, steady state analysis, and Jacobian (see cvodesim.h)
ssa -- stochastic simulation algorithm (see ssa.h)
ga -- the basic genetic algorithm functions
netga -- genetic algorithm specifically for evolving biological networks
gabistabilize -- find parameters for making a system bistable
loops -- find loops in an interaction matrix

--------------------------------------------------------
--------------------------------------------------------
What executables get generated
--------------------------------------------------------
--------------------------------------------------------

--------------------------------------
simple examples of using netga
--------------------------------------

evolveOscillator, evolveNoiseDamper

--------------------------------------
more complex use of netga and ga
--------------------------------------

evolveChemotaxis
bistabilize

--------------------------------------
finding loops in a Jacobian matrix
--------------------------------------

findLoops
