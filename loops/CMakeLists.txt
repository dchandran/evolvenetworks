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
