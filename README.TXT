Network Evolution Library (2009)
Author: Deepak Chandran (dchandran1@gmail.com)
LICENSE: open source BSD (see COPYRIGHT.TXT)

--------------------------------------------------------
--------------------------------------------------------
How to Build the project
--------------------------------------------------------
--------------------------------------------------------

1. install CMake (www.cmake.org) 
2. Run the CMake GUI, located in the CMake/bin folder. The GUI will ask you for the Source and Build directories.
3. Select the network evolution library directory as the Source directory
4. Select the (network evolution library)/BUILD directory as the BUILD directory (create this folder)
5. Click the "Configure" button
6. Select the C compiler that you would like to use *
7. May need to "Configure" again during the first build
8. After configure is done, click "Generate"
9. CMake has not created the platform-specific makefile, which can be used to automatically compile everything. 
10. Go to the (network evolution library)/BUILD folder and type "make" (or the equivalent command for other compilers, such as nmake or mingw32-make)
11. The (network evolution library)/BUILD/bin will contain all the binaries and (network evolution library)/BUILD/lib will contain all the libraries. 

* Note: the CMake project has been tested using gcc in Unix and MinGW in Windows. Visual Studio will also work for most of the code, but there may be issues when compiling the /GUI classes which use Qt. 

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
