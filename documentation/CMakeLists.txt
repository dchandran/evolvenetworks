CONFIGURE_FILE ( 
  ${NETWORKEVOLUTION_SOURCE_DIR}/documentation/Doxyfile.txt.in 
  ${NETWORKEVOLUTION_SOURCE_DIR}/documentation/Doxyfile.txt
) 

ADD_CUSTOM_TARGET(doc ALL 
  ${DOXYGEN_EXECUTABLE} 
  ${NETWORKEVOLUTION_BINARY_DIR}/documentation/Doxyfile.txt 
  DEPENDS 
  ${NETWORKEVOLUTION_BINARY_DIR}/documentation/Doxyfile.txt 
) 

