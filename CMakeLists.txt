cmake_minimum_required(VERSION 2.0)

PROJECT( NETWORKEVOLUTION )

#---------------------------------------------------------------------------
# Doxygen documentation
#---------------------------------------------------------------------------
OPTION( DOXYGENDOC "Generate Doxygen Documentation" OFF )

IF( DOXYGENDOC )

  FIND_PACKAGE( Doxygen )
    IF( DOXYGEN_FOUND )
      ADD_SUBDIRECTORY( documentation )
    ELSE( DOXYGEN_FOUND )
      MESSAGE( SEND_ERROR
        "Doxygen MUST be installed 
          (visit http://www.stack.nl/~dimitri/doxygen/)" )
    ENDIF( DOXYGEN_FOUND )
	
ENDIF( DOXYGENDOC )

#---------------------------------------------------------------------------

INCLUDE_DIRECTORIES( BEFORE
    ${NETWORKEVOLUTION_SOURCE_DIR}/external/
	${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode
	${NETWORKEVOLUTION_SOURCE_DIR}/external/optim/
	${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/
	${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/cvode
	${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/nvec_ser
	${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/sundials
	${NETWORKEVOLUTION_SOURCE_DIR}/simulation
	${NETWORKEVOLUTION_SOURCE_DIR}/GA
	${NETWORKEVOLUTION_SOURCE_DIR}/bistability
	${NETWORKEVOLUTION_SOURCE_DIR}/loops
	${NETWORKEVOLUTION_SOURCE_DIR}/chemotaxis
	${NETWORKEVOLUTION_SOURCE_DIR}/lib
	${NETWORKEVOLUTION_SOURCE_DIR}/bin
)
LINK_DIRECTORIES(
	${NETWORKEVOLUTION_SOURCE_DIR}/lib
	${NETWORKEVOLUTION_SOURCE_DIR}/bin
)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# LIBRARIES
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

#ADD_SUBDIRECTORY( external )

#---------------------------------------------------------------------------
#Nelder-Mead algorithm
#---------------------------------------------------------------------------

FILE( GLOB OPTIM_SRC ${NETWORKEVOLUTION_SOURCE_DIR}/external/optim/*.c )

ADD_LIBRARY( optim 
  STATIC
  ${OPTIM_SRC}
)

SET_TARGET_PROPERTIES(optim PROPERTIES LINKER_LANGUAGE C)

#------------------------------
# Output paths
#------------------------------

SET(EXECUTABLE_OUTPUT_PATH
 ${NETWORKEVOLUTION_BINARY_DIR}/bin
 CACHE PATH "Single output directory for building all executables."
)

SET(LIBRARY_OUTPUT_PATH
 ${NETWORKEVOLUTION_BINARY_DIR}/lib
 CACHE PATH "Single output directory for building all libraries."
)

#---------------------------------------------------------------------------
#Wrapper functions for Sundials CVODE
#---------------------------------------------------------------------------

FILE( GLOB CVODE_SRC ${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/cvode/*.c )
FILE( GLOB NVECSER_SRC ${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/nvec_ser/*.c )
FILE( GLOB SUNDIALS_SRC ${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/sundials/*.c )

ADD_LIBRARY( ode
  STATIC
  ${CVODE_SRC}
  ${SUNDIALS_SRC}
  ${NVECSER_SRC}
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/cvodesim.c
)
SET_TARGET_PROPERTIES(ode PROPERTIES LINKER_LANGUAGE C)

#---------------------------------------------------------------------------
#stochastic simulation
#---------------------------------------------------------------------------

#Gillespie

ADD_LIBRARY( ssa
  STATIC
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/ssa.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/mtrand.c
)

#Multi-cell stochastic simulation

ADD_LIBRARY( cells_ssa
  STATIC
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/cells_ssa.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/mtrand.c
)

#---------------------------------------------------------------------------
# Genetic algorithm (GA)
#---------------------------------------------------------------------------

#generic GA

ADD_LIBRARY( ga
  STATIC
  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/ga.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/mtrand.c
)

#Evolve network library

ADD_LIBRARY( netga
  STATIC
  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/massActionNetwork.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/proteinInteractionNetwork.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/enzymeNetwork.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/geneRegulationNetwork.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/reactionNetwork.c
)

TARGET_LINK_LIBRARIES( netga
  netga
  ga
  ode
  ssa
)

#---------------------------------------------------------------------------
# Bistability finder
#---------------------------------------------------------------------------

ADD_LIBRARY( gaBistabilize
  STATIC
  ${NETWORKEVOLUTION_SOURCE_DIR}/bistability/ga_bistable.c
)

TARGET_LINK_LIBRARIES( gaBistabilize
  gaBistabilize
  netga
)

#---------------------------------------------------------------------------
# Loops in interaction matrix finder
#---------------------------------------------------------------------------

ADD_LIBRARY( loops
  STATIC
  ${NETWORKEVOLUTION_SOURCE_DIR}/loops/loops.c
)

TARGET_LINK_LIBRARIES( loops
  loops
  ode
)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# EXECUTABLES 
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

# Use GA library

ADD_EXECUTABLE( evolveChemotaxis ${NETWORKEVOLUTION_SOURCE_DIR}/chemotaxis/evolve_chemotaxis.c )
TARGET_LINK_LIBRARIES( evolveChemotaxis netga loops )

ADD_EXECUTABLE( bistabilize ${NETWORKEVOLUTION_SOURCE_DIR}/bistability/test_bistable.c )
TARGET_LINK_LIBRARIES( bistabilize gaBistabilize optim )

ADD_EXECUTABLE( evolveOscillator ${NETWORKEVOLUTION_SOURCE_DIR}/examples/ga_oscillator.c )
TARGET_LINK_LIBRARIES( evolveOscillator netga )

ADD_EXECUTABLE( evolveNoiseDamper ${NETWORKEVOLUTION_SOURCE_DIR}/examples/ga_reduceNoise.c )
TARGET_LINK_LIBRARIES( evolveNoiseDamper netga )

ADD_EXECUTABLE( evolveXOR ${NETWORKEVOLUTION_SOURCE_DIR}/examples/ga_logicGate.c )
TARGET_LINK_LIBRARIES( evolveXOR netga )

# Not related to GA

ADD_EXECUTABLE( findLoops ${NETWORKEVOLUTION_SOURCE_DIR}/loops/test_loops.c )
TARGET_LINK_LIBRARIES( findLoops loops )
