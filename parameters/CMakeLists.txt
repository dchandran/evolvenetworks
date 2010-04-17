#---------------------------------------------------------------------------
# Parameter related
#---------------------------------------------------------------------------

ADD_LIBRARY( gaBistabilize
  STATIC
  ${NETWORKEVOLUTION_SOURCE_DIR}/parameters/ga_bistable.c
)

TARGET_LINK_LIBRARIES( gaBistabilize
  gaBistabilize
  ga
  optim
)

ADD_EXECUTABLE( parameterInformation ${NETWORKEVOLUTION_SOURCE_DIR}/parameters/ParameterInformationContent.c )

TARGET_LINK_LIBRARIES( parameterInformation ode ssa )

ADD_EXECUTABLE( parameterMC ${NETWORKEVOLUTION_SOURCE_DIR}/parameters/ParameterMC.cpp )
TARGET_LINK_LIBRARIES( parameterMC ode ssa ${Boost_THREAD_LIBRARY} )

