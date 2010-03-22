#---------------------------------------------------------------------------
# GAlib from lancet.mit.edu/ga/
#---------------------------------------------------------------------------

FILE(GLOB GALIB_FILES ${NETWORKEVOLUTION_SOURCE_DIR}/galib/ga/*.cpp )

ADD_LIBRARY( galib
  STATIC
  ${GALIB_FILES}
)
