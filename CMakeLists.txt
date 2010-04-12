cmake_minimum_required(VERSION 2.0)

PROJECT( NETWORKEVOLUTION )
IF( COMMAND CMAKE_POLICY )
  CMAKE_POLICY( VERSION 2.4 )
  CMAKE_POLICY( SET CMP0005 NEW )
  CMAKE_POLICY( SET CMP0003 NEW )
  CMAKE_POLICY( SET CMP0002 OLD )
ENDIF(COMMAND CMAKE_POLICY)

SET(CMAKE_BUILD_TYPE "Release")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DQT_NO_DEBUG -DQT_XML_LIB -DQT_GUI_LIB -DQT_CORE_LIB -DQT_THREAD_SUPPORT")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --export-all-symbols")

IF (UNIX AND NOT APPLE)
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lrt")
	SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -lrt")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lrt")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
ENDIF(UNIX AND NOT APPLE)

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
	${NETWORKEVOLUTION_SOURCE_DIR}/external/lapack/INCLUDE
	${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/
	${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/cvode
	${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/nvec_ser
	${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/sundials
	${NETWORKEVOLUTION_SOURCE_DIR}/external/antimony/src
	${NETWORKEVOLUTION_SOURCE_DIR}/external/sbml
	${NETWORKEVOLUTION_SOURCE_DIR}/external/muparser
	${NETWORKEVOLUTION_SOURCE_DIR}/simulation
	
	${NETWORKEVOLUTION_SOURCE_DIR}/GA
	${NETWORKEVOLUTION_SOURCE_DIR}/galib
	${NETWORKEVOLUTION_SOURCE_DIR}/galib/ga
	${NETWORKEVOLUTION_SOURCE_DIR}/GA/blocks
	${NETWORKEVOLUTION_SOURCE_DIR}/parameters
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

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

ADD_SUBDIRECTORY( external )
ADD_SUBDIRECTORY( simulation )
ADD_SUBDIRECTORY( GA )
ADD_SUBDIRECTORY( galib )
ADD_SUBDIRECTORY( parameters )
ADD_SUBDIRECTORY( loops )

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
# SWIG
#---------------------------------------------------------------------------

OPTION( SWIG_BINDINGS "make SWIG bindings" OFF )

IF (SWIG_BINDINGS)
    FIND_PACKAGE(SWIG)
    INCLUDE(${SWIG_USE_FILE})

    ADD_CUSTOM_COMMAND( 
	    OUTPUT NotAFile
	    WORKING_DIRECTORY ${NETWORKEVOLUTION_SOURCE_DIR}/simulation
        COMMAND `swig -python sim.i`
		WORKING_DIRECTORY ${NETWORKEVOLUTION_SOURCE_DIR}/GA
        COMMAND `swig -python ga.i`
    )
	ADD_LIBRARY( _sim
      SHARED
	  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/sim_wrap.c
   )
   TARGET_LINK_LIBRARIES( _sim
      _sim
	  ode
	  ssa
	  ${PYTHON_LIBRARIES}
   )
	ADD_LIBRARY( _ga
      SHARED
	  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/ga_wrap.c
   )
   TARGET_LINK_LIBRARIES( _ga
      _ga
      ga
	  ${PYTHON_LIBRARIES}
   )
ENDIF(SWIG_BINDINGS)

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# EXECUTABLES 
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

#testing

ADD_EXECUTABLE( test_ode ${NETWORKEVOLUTION_SOURCE_DIR}/examples/test_cvode.c )
TARGET_LINK_LIBRARIES( test_ode ode )

ADD_EXECUTABLE( test_sbml_sim ${NETWORKEVOLUTION_SOURCE_DIR}/examples/optim_sbml.cpp ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/mtrand.c)
TARGET_LINK_LIBRARIES( test_sbml_sim sbml_sim galib)

# Use GA library

#ADD_EXECUTABLE( evolveChemotaxis ${NETWORKEVOLUTION_SOURCE_DIR}/chemotaxis/evolve_chemotaxis.c )
#TARGET_LINK_LIBRARIES( evolveChemotaxis ga loops )

ADD_EXECUTABLE( bistabilize ${NETWORKEVOLUTION_SOURCE_DIR}/parameters/test_bistable.c )
TARGET_LINK_LIBRARIES( bistabilize gaBistabilize )

ADD_EXECUTABLE( evolveOscillator ${NETWORKEVOLUTION_SOURCE_DIR}/examples/ga_oscillator.c )
TARGET_LINK_LIBRARIES( evolveOscillator ga )

ADD_EXECUTABLE( evolveNoiseDamper ${NETWORKEVOLUTION_SOURCE_DIR}/examples/ga_reduceNoise.c )
TARGET_LINK_LIBRARIES( evolveNoiseDamper ga )

ADD_EXECUTABLE( evolveXOR ${NETWORKEVOLUTION_SOURCE_DIR}/examples/ga_logicGate.c )
TARGET_LINK_LIBRARIES( evolveXOR ga )

# Not related to GA

#ADD_EXECUTABLE( findLoops ${NETWORKEVOLUTION_SOURCE_DIR}/loops/test_loops.c )
#TARGET_LINK_LIBRARIES( findLoops loops )

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
      ADD_SUBDIRECTORY( network_design )
  ENDIF( QT4_FOUND )
ENDIF ( BUILD_GUI_APP )
