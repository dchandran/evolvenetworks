#---------------------------------------------------------------------------
# Parameter related
#---------------------------------------------------------------------------

find_package( Boost 1.41 COMPONENTS threads)
set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)

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

include_directories ( ${Boost_INCLUDE_DIRS} )
ADD_EXECUTABLE( parameterMC ${NETWORKEVOLUTION_SOURCE_DIR}/parameters/ParameterMC.cpp )
TARGET_LINK_LIBRARIES( parameterMC ode ssa ${Boost_THREAD_LIBRARY} )

