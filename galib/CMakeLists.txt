#---------------------------------------------------------------------------
# GAlib from lancet.mit.edu/ga/
#---------------------------------------------------------------------------

ADD_LIBRARY( galib
  STATIC
  ${NETWORKEVOLUTION_SOURCE_DIR}/galib/ga/*.C
)
