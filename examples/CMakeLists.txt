#testing

ADD_EXECUTABLE( test_ode ${NETWORKEVOLUTION_SOURCE_DIR}/examples/test_cvode.c )
TARGET_LINK_LIBRARIES( test_ode ode )

ADD_EXECUTABLE( test_sbml ${NETWORKEVOLUTION_SOURCE_DIR}/examples/test_sbml.cpp )
TARGET_LINK_LIBRARIES( test_sbml sbml_sim )

ADD_EXECUTABLE( optim_sbml ${NETWORKEVOLUTION_SOURCE_DIR}/examples/optim_sbml.cpp)
TARGET_LINK_LIBRARIES( optim_sbml sbml_sim galib)

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

ADD_EXECUTABLE( test_modularga ${NETWORKEVOLUTION_SOURCE_DIR}/examples/test_modularga.cpp)
TARGET_LINK_LIBRARIES( test_modularga GAmodular)

