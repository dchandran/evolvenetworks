if(COMMAND cmake_policy)
  cmake_policy(SET CMP0005 NEW)
endif(COMMAND cmake_policy)

IF (WIN32)

	IF (MINGW)
		SET (LIBXML2_LIBRARIES ${NETWORKEVOLUTION_SOURCE_DIR}/sbml/win32/libxml2/lib/libxml2.dll.a)
	ELSE(MINGW)
		SET (LIBXML2_LIBRARIES ${NETWORKEVOLUTION_SOURCE_DIR}/sbml/win32/libxml2/lib/libxml2.lib)
	ENDIF(MINGW)
	SET (LIBXML2_INCLUDE_DIR ${NETWORKEVOLUTION_SOURCE_DIR}/sbml/win32/libxml2/include)
	
ELSE(WIN32)	

	FIND_PACKAGE(LIBXML2)
	
ENDIF(WIN32)

INCLUDE_DIRECTORIES( 
	BEFORE
	${LIBXML2_INCLUDE_DIR}
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/xml 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/util 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/validator 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/validator/constraints/
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/units 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/sbml
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/sbml/layout 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/math 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/common 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/compress 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/annotation
)

#ADD_DEFINITIONS(-DLIBSBML_STATIC)
#ADD_DEFINITIONS(-D__cplusplus)
ADD_DEFINITIONS(-DLIBSBML_EXPORTS)
ADD_DEFINITIONS(-DLIBLAX_EXPORTS)

ADD_DEFINITIONS(-DPACKAGE_NAME=\"libSBML\")
ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-DUSE_LIBXML)
ADD_DEFINITIONS(-DUSE_LAYOUT)

IF (WIN32)
	ADD_DEFINITIONS(-DWIN32)
	ADD_DEFINITIONS(-D_WIN32)
	ADD_DEFINITIONS(-D_WINDOWS)
	ADD_DEFINITIONS(-D_USRDLL)
	ADD_DEFINITIONS(-D_EXPORTS)
	ADD_DEFINITIONS(-DIN_LIBXML)
	IF (MINGW)
		ADD_DEFINITIONS(-DCYGWIN)
	ENDIF(MINGW)
ENDIF (WIN32)

FILE( GLOB SBML_SRC_FILES
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/xml/XML*.cpp 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/xml/XML*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/xml/LibXML*.cpp 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/xml/LibXML*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/util/*.cpp 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/util/*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/validator/*.cpp 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/validator/*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/validator/constraints/*.cpp
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/validator/constraints/*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/units/*.cpp 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/units/*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/math/*.cpp
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/math/*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/annotation/*.cpp 
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/annotation/*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/sbml/*.cpp
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/sbml/*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/sbml/layout/*.cpp
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/sbml/layout/*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/common/*.cpp
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/common/*.c
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/compress/CompressCommon.cpp
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/compress/InputDecompressor.cpp
	${NETWORKEVOLUTION_SOURCE_DIR}/sbml/compress/OutputCompressor.cpp
)

ADD_LIBRARY( sbml
  #SHARED
  STATIC
  ${SBML_SRC_FILES}
)

TARGET_LINK_LIBRARIES( sbml
  sbml
  ${LIBXML2_LIBRARIES}
)


IF ( WIN32 )
   FILE( GLOB ${NETWORKEVOLUTION_SOURCE_DIR}/sbml/win32/libxml2/lib/*.dll DESTINATION bin/../ )
   INSTALL(TARGETS DESTINATION bin/../)
#  INSTALL(TARGETS sbml DESTINATION bin/../)
#ELSE ( WIN32 )
  #INSTALL(TARGETS sbml LIBRARY DESTINATION bin/../)
ENDIF( WIN32 )


INCLUDE_DIRECTORIES(
	${NETWORKEVOLUTION_SOURCE_DIR}/ExternalCode/Antimony/src
)

SET( ANTIMONY_SRC
	src/Antimony.cpp
	src/antimony_api.cpp
	src/antimony.tab.cpp
	src/dnastrand.cpp
	src/event.cpp
	src/formula.cpp
	src/module.cpp
	src/reactantlist.cpp
	src/reaction.cpp
	src/registry.cpp
	src/sbmlx.cpp
	src/stringx.cpp
	src/typex.cpp
	src/userfunction.cpp
	src/variable.cpp
)
   
#ADD_DEFINITIONS(-DNSBML)
ADD_DEFINITIONS(-DNCELLML)
ADD_DEFINITIONS(-DLIB_EXPORTS)

ADD_LIBRARY ( antimony
  SHARED
  ${ANTIMONY_SRC}
)


TARGET_LINK_LIBRARIES( antimony
  antimony
  sbml
)

IF ( WIN32 )
  INSTALL(TARGETS antimony DESTINATION bin/../)
ELSE ( WIN32 )
  INSTALL(TARGETS antimony LIBRARY DESTINATION bin/../)
ENDIF( WIN32 )
