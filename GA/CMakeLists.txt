#---------------------------------------------------------------------------
# Genetic algorithm (GA)
#---------------------------------------------------------------------------

#modular GA

ADD_EXECUTABLE( antimonyToFunctions ${NETWORKEVOLUTION_SOURCE_DIR}/GA/antimonyToFunctions.c )
TARGET_LINK_LIBRARIES( antimonyToFunctions antimony )

ADD_LIBRARY( ga
  STATIC
  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/ga.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/simulation/mtrand.c
  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/blocks.c
)

TARGET_LINK_LIBRARIES( ga ga ode ssa)

ADD_EXECUTABLE( testblocks ${NETWORKEVOLUTION_SOURCE_DIR}/GA/testBlocks.c )
TARGET_LINK_LIBRARIES( testblocks ga )

