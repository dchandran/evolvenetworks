if(COMMAND cmake_policy) 
cmake_policy(SET CMP0003 NEW) 
endif(COMMAND cmake_policy) 

#FILE( GLOB INPUTS_HEADER_FILES GenerateGAInputs.h CodeEditor.h SyntaxHighlighter.h )
#FILE( GLOB INPUTS_SOURCE_FILES GenerateGAInputs.cpp CodeEditor.cpp SyntaxHighlighter.cpp )

#FILE( GLOB RUN_HEADER_FILES RunEvolution.h CodeEditor.h SyntaxHighlighter.h )
#FILE( GLOB RUN_SOURCE_FILES RunEvolution.cpp CodeEditor.cpp SyntaxHighlighter.cpp )

FILE( GLOB MODULES_HEADER_FILES ModuleSelector.h )
FILE( GLOB MODULES_SOURCE_FILES ModuleSelector.cpp )

#QT4_WRAP_CPP( INPUTS_HEADER_MOC ${INPUTS_HEADER_FILES} )
#QT4_WRAP_CPP( RUN_HEADER_MOC ${RUN_HEADER_FILES} )
QT4_WRAP_CPP( MODULES_HEADER_MOC ${MODULES_HEADER_FILES} )

QT4_ADD_RESOURCES( NETWORK_EVOLVE_GUI_QRC networkevolve.qrc )

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

  ADD_EXECUTABLE( ModularEvolution 
    MACOSX_BUNDLE
    ${MODULES_SOURCE_FILES}
    ${MODULES_HEADER_MOC}
    ${NETWORK_EVOLVE_GUI_QRC}
  )
  
  TARGET_LINK_LIBRARIES( ModularEvolution
	${QT_LIBRARIES}
	ga
  )

  ADD_CUSTOM_COMMAND( TARGET GenerateInputs POST_BUILD
      COMMAND mkdir ARGS -p
        ${EXECUTABLE_OUTPUT_PATH}/ModularEvolution.app/Contents/Resources
        ${EXECUTABLE_OUTPUT_PATH}/ModularEvolution.app/Contents/Frameworks
        ${EXECUTABLE_OUTPUT_PATH}/ModularEvolution.app/Contents/MacOS
      COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtCore.framework
        ${EXECUTABLE_OUTPUT_PATH}/ModularEvolution.app/Contents/Frameworks/
	  COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtGui.framework
        ${EXECUTABLE_OUTPUT_PATH}/ModularEvolution.app/Contents/Frameworks/
	  COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtXml.framework
        ${EXECUTABLE_OUTPUT_PATH}/ModularEvolution.app/Contents/Frameworks/
      COMMAND cp ARGS ${MACOSX_NETWORKEVOLUTION_BUNDLE_ICON_FILE}
        ${EXECUTABLE_OUTPUT_PATH}/ModularEvolution.app/Contents/Resources
  )

ELSE( APPLE AND BUILD_BUNDLE )
  IF( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )
  	
	ADD_EXECUTABLE( ModularEvolution 
	  ${MODULES_SOURCE_FILES}
      ${MODULES_HEADER_MOC}
      ${NETWORK_EVOLVE_GUI_QRC}
	)
  
	TARGET_LINK_LIBRARIES( ModularEvolution
	  ${QT_LIBRARIES}
	  ga
	)

  ELSE( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )
    IF( WIN32 )
		  IF( MINGW )

			# resource compilation for mingw
  			ADD_CUSTOM_COMMAND( OUTPUT 
			${CMAKE_CURRENT_BINARY_DIR}/networkevolve_rc.o
				COMMAND windres.exe -I${CMAKE_CURRENT_SOURCE_DIR} 
					-i${CMAKE_CURRENT_SOURCE_DIR}/networkevolve.rc
					-o ${CMAKE_CURRENT_BINARY_DIR}/networkevolve_rc.o )		
  			SET( NETWORKEVOLVE_RC ${CMAKE_CURRENT_BINARY_DIR}/networkevolve_rc.o )

  		ELSE( MINGW )

	  		SET(NETWORKEVOLVE_RC networkevolve.rc)

		ENDIF( MINGW )
	
	ADD_EXECUTABLE( ModularEvolution WIN32
	  ${MODULES_SOURCE_FILES}
      ${MODULES_HEADER_MOC}
      ${NETWORK_EVOLVE_GUI_QRC}
	)
  
	TARGET_LINK_LIBRARIES( ModularEvolution
	  ${QT_LIBRARIES}
	  ga
	)

    ENDIF( WIN32 )
  ENDIF( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )

ENDIF( APPLE AND BUILD_BUNDLE )

#---------------------------------------------------------------------
#  INSTALL executable and required files and folders
#---------------------------------------------------------------------
   
INSTALL(TARGETS ModularEvolution RUNTIME DESTINATION bin/../)
   
INSTALL(FILES networkevolve.qss DESTINATION bin/../)
