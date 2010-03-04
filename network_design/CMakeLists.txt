if(COMMAND cmake_policy) 
cmake_policy(SET CMP0003 NEW) 
endif(COMMAND cmake_policy) 

FILE( GLOB FFNDESIGN_HEADER_FILES *.h )
FILE( GLOB FFNDESIGN_SOURCE_FILES *.cpp *.c)

QT4_WRAP_CPP( FFNDESIGN_HEADER_MOC ${FFNDESIGN_HEADER_FILES} )

IF( APPLE AND BUILD_BUNDLE )
  SET( MACOSX_NETWORKEVOLUTION_BUNDLE_ICON_FILE
      ${NETWORKEVOLUTION_SOURCE_DIR}/GUI/networkevolve.icns
  )
  SET( MACOSX_BUNDLE_SHORT_VERSION_STRING
      ${NETWORKEVOLUTION_VERSION}
  )
  SET( MACOSX_BUNDLE_VERSION
      ${NETWORKEVOLUTION_VERSION}
  )
  SET( MACOSX_BUNDLE_LONG_VERSION_STRING
      Version ${NETWORKEVOLUTIONL_VERSION}
  )

  ADD_EXECUTABLE( FeedforwardGRNDesign 
    MACOSX_BUNDLE
    ${FFNDESIGN_SOURCE_FILES}
    ${FFNDESIGN_HEADER_MOC}
  )
  
  TARGET_LINK_LIBRARIES( FeedforwardGRNDesign
	${QT_LIBRARIES}
	antimony
  )

  ADD_CUSTOM_COMMAND( TARGET GenerateInputs POST_BUILD
      COMMAND mkdir ARGS -p
        ${EXECUTABLE_OUTPUT_PATH}/FeedforwardGRNDesign.app/Contents/Resources
        ${EXECUTABLE_OUTPUT_PATH}/FeedforwardGRNDesign.app/Contents/Frameworks
        ${EXECUTABLE_OUTPUT_PATH}/FeedforwardGRNDesign.app/Contents/MacOS
      COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtCore.framework
        ${EXECUTABLE_OUTPUT_PATH}/FeedforwardGRNDesign.app/Contents/Frameworks/
	  COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtGui.framework
        ${EXECUTABLE_OUTPUT_PATH}/FeedforwardGRNDesign.app/Contents/Frameworks/
	  COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtXml.framework
        ${EXECUTABLE_OUTPUT_PATH}/FeedforwardGRNDesign.app/Contents/Frameworks/
      COMMAND cp ARGS ${MACOSX_NETWORKEVOLUTION_BUNDLE_ICON_FILE}
        ${EXECUTABLE_OUTPUT_PATH}/FeedforwardGRNDesign.app/Contents/Resources
  )

ELSE( APPLE AND BUILD_BUNDLE )
  IF( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )
  	
	ADD_EXECUTABLE( FeedforwardGRNDesign 
	    ${FFNDESIGN_SOURCE_FILES}
		${FFNDESIGN_HEADER_MOC}
	)
  
	TARGET_LINK_LIBRARIES( FeedforwardGRNDesign
	  ${QT_LIBRARIES}
	  antimony
	)

  ELSE( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )
    IF( WIN32 )
	
	ADD_EXECUTABLE( FeedforwardGRNDesign WIN32
	    ${FFNDESIGN_SOURCE_FILES}
		${FFNDESIGN_HEADER_MOC}
	)
  
	TARGET_LINK_LIBRARIES( FeedforwardGRNDesign
	  ${QT_LIBRARIES}
	  antimony
	)

    ENDIF( WIN32 )
  ENDIF( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )

ENDIF( APPLE AND BUILD_BUNDLE )

#---------------------------------------------------------------------
#  INSTALL executable and required files and folders
#---------------------------------------------------------------------
   
INSTALL(TARGETS FeedforwardGRNDesign RUNTIME DESTINATION bin/../)
