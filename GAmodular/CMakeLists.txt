#---------------------------------------------------------------------------
# GA modular
#---------------------------------------------------------------------------

FILE(GLOB GAMODULAR_FILES ${NETWORKEVOLUTION_SOURCE_DIR}/GAmodular/*.cpp )

ADD_LIBRARY( GAmodular
  STATIC
  ${GAMODULAR_FILES}
)

TARGET_LINK_LIBRARIES( GAmodular
   GAmodular
   sbml_sim
   ${Boost_REGEX_LIBRARY}
#   ${Boost_THREADS_LIBRARY}
)
