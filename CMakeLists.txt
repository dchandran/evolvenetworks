cmake_minimum_required(VERSION 2.0)

PROJECT( NETWORKEVOLUTION )

#---------------------------------------------------------------------------
# Python
#---------------------------------------------------------------------------
FIND_PACKAGE( PythonLibs )

IF( NOT PYTHONLIBS_FOUND )
  MESSAGE( SEND_ERROR "TinkerCell requires python to be installed" )
ENDIF( NOT PYTHONLIBS_FOUND )

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
	${PYTHON_INCLUDE_PATH}
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
	${NETWORKEVOLUTION_SOURCE_DIR}/GUI
)
LINK_DIRECTORIES(
	${NETWORKEVOLUTION_SOURCE_DIR}/lib
	${NETWORKEVOLUTION_SOURCE_DIR}/bin
	${NETWORKEVOLUTION_SOURCE_DIR}/GUI
	${NETWORKEVOLUTION_SOURCE_DIR}
)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# LIBRARIES
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

ADD_SUBDIRECTORY( external/lapack )

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

TARGET_LINK_LIBRARIES( ode ode  optim )

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


#---------------------------------------------------------------------------
# Evolve network library (using GA)
#---------------------------------------------------------------------------

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
# SWIG
#---------------------------------------------------------------------------

OPTION( SWIG_BINDINGS "make SWIG bindings" OFF )

IF (SWIG_BINDINGS)
    FIND_PACKAGE(SWIG)
    INCLUDE(${SWIG_USE_FILE})

    ADD_CUSTOM_COMMAND( 
	    OUTPUT ga_wrap.c
	    WORKING_DIRECTORY ${NETWORKEVOLUTION_SOURCE_DIR}/GA
        COMMAND `swig -python ga.i`
    )
	ADD_LIBRARY( pyga
      SHARED
	  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/ga_wrap.c
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
ENDIF(SWIG_BINDINGS)
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
TARGET_LINK_LIBRARIES( bistabilize gaBistabilize )

ADD_EXECUTABLE( evolveOscillator ${NETWORKEVOLUTION_SOURCE_DIR}/examples/ga_oscillator.c )
TARGET_LINK_LIBRARIES( evolveOscillator netga )

ADD_EXECUTABLE( evolveNoiseDamper ${NETWORKEVOLUTION_SOURCE_DIR}/examples/ga_reduceNoise.c )
TARGET_LINK_LIBRARIES( evolveNoiseDamper netga )

ADD_EXECUTABLE( evolveXOR ${NETWORKEVOLUTION_SOURCE_DIR}/examples/ga_logicGate.c )
TARGET_LINK_LIBRARIES( evolveXOR netga )

# Not related to GA

ADD_EXECUTABLE( findLoops ${NETWORKEVOLUTION_SOURCE_DIR}/loops/test_loops.c )
TARGET_LINK_LIBRARIES( findLoops loops )

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# GUI using QT
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

OPTION( BUILD_GUI_APP
    "Built the GUI for the evolution program. Qt Toolkit must be installed and configured properly." 
    OFF )
	
IF ( BUILD_GUI_APP )
  FIND_PACKAGE( Qt4 )
  IF( QT4_FOUND )

    INCLUDE_DIRECTORIES( BEFORE ${QT_INCLUDES} )
      IF( QT_USE_FILE )
	    SET( QT_USE_QTCORE TRUE )
	    SET( QT_USE_QTGUI TRUE )
        INCLUDE( ${QT_USE_FILE} )
      ELSE( QT_USE_FILE )
        SET( QT_LIBRARIES ${QT_QT_LIBRARY} )
      ENDIF( QT_USE_FILE )
      ADD_SUBDIRECTORY( GUI )
  ENDIF( QT4_FOUND )
ENDIF ( BUILD_GUI_APP )
