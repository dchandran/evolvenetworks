#---------------------------------------------------------------------------
# Wrapper functions for Sundials CVODE
#---------------------------------------------------------------------------

FILE( GLOB CVODE_SRC ${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/cvode/*.c )
FILE( GLOB NVECSER_SRC ${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/nvec_ser/*.c )
FILE( GLOB SUNDIALS_SRC ${NETWORKEVOLUTION_SOURCE_DIR}/external/cvode_src/sundials/*.c )

ADD_DEFINITIONS(-DNO_FPRINTF_OUTPUT)
ADD_LIBRARY( ode
  STATIC
  ${CVODE_SRC}
  ${SUNDIALS_SRC}
  ${NVECSER_SRC}
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/eigen.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/cvodesim.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/external/optim/mat.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/external/optim/neldermead.c
)

TARGET_LINK_LIBRARIES( ode ode lapack )

SET_TARGET_PROPERTIES(ode PROPERTIES LINKER_LANGUAGE C)

#---------------------------------------------------------------------------
# stochastic simulation
#---------------------------------------------------------------------------

#Gillespie, Langevin, and custom multi-cell simulation

ADD_LIBRARY( ssa
  STATIC
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/ssa.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/cells_ssa.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/langevin.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/mtrand.c
)

SET_TARGET_PROPERTIES(ssa PROPERTIES LINKER_LANGUAGE C)

ADD_LIBRARY( sbml_sim
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/sbml_sim.cpp
)

TARGET_LINK_LIBRARIES( sbml_sim 
   sbml_sim 
   sbml 
   ode 
   ssa 
   muparser
)



