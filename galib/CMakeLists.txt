#---------------------------------------------------------------------------
# GAlib from lancet.mit.edu/ga/
#---------------------------------------------------------------------------

FILE(GLOB GALIB_FILES ${NETWORKEVOLUTION_SOURCE_DIR}/galib/ga/*.C )

ADD_LIBRARY( galib
  STATIC
  ${GALIB_FILES}
)
