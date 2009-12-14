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

TARGET_LINK_LIBRARIES( ga )

#---------------------------------------------------------------------------
# Evolve network library (using GA)
#---------------------------------------------------------------------------

#ADD_LIBRARY( netga
#  SHARED
#  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/massActionNetwork.c
#  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/proteinInteractionNetwork.c
#  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/enzymeNetwork.c
#  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/geneRegulationNetwork.c
#  ${NETWORKEVOLUTION_SOURCE_DIR}/GA/reactionNetwork.c
#)

#TARGET_LINK_LIBRARIES( netga
#  netga
#  ga
#  ode
#  ssa
#)