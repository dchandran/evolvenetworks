FILE( GLOB HEADER_FILES *.h )

FILE( GLOB SOURCE_FILES *.cpp )

QT4_WRAP_CPP( HEADER_MOC ${HEADER_FILES} )

QT4_ADD_RESOURCES( NETWORK_EVOLVE_GUI_QRC networkevolve.qrc )

IF( APPLE AND BUILD_BUNDLE )
  SET( MACOSX_NETWORKEVOLUTION_BUNDLE_ICON_FILE
      ${NETWORKEVOLUTION_SOURCE_DIR}/Main/networkevolve.icns
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

  ADD_EXECUTABLE( NetworkEvolverGUI 
    MACOSX_BUNDLE
    ${SOURCE_FILES}
    ${HEADER_MOC}
    ${NETWORK_EVOLVE_GUI_QRC}
  )
  
  TARGET_LINK_LIBRARIES( NetworkEvolverGUI
	${QT_LIBRARIES}
  )

  ADD_CUSTOM_COMMAND( TARGET NetworkEvolverGUI POST_BUILD
      COMMAND mkdir ARGS -p
        ${EXECUTABLE_OUTPUT_PATH}/NetworkEvolverGUI.app/Contents/Resources
        ${EXECUTABLE_OUTPUT_PATH}/NetworkEvolverGUI.app/Contents/Frameworks
        ${EXECUTABLE_OUTPUT_PATH}/NetworkEvolverGUI.app/Contents/MacOS
      COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtXml.framework
        ${EXECUTABLE_OUTPUT_PATH}/NetworkEvolverGUI.app/Contents/Frameworks/
      COMMAND cp ARGS ${MACOSX_NETWORKEVOLUTION_BUNDLE_ICON_FILE}
        ${EXECUTABLE_OUTPUT_PATH}/NetworkEvolverGUI.app/Contents/Resources
  )

ELSE( APPLE AND BUILD_BUNDLE )
  IF( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )
    ADD_EXECUTABLE( NetworkEvolverGUI 
      ${SOURCE_FILES}
      ${HEADER_MOC}
      ${NETWORK_EVOLVE_GUI_QRC}
    )
	
	TARGET_LINK_LIBRARIES( NetworkEvolverGUI
		${QT_LIBRARIES}
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
	
    ADD_EXECUTABLE( NetworkEvolverGUI WIN32 
      ${SOURCE_FILES}
      ${HEADER_MOC}
      ${NETWORK_EVOLVE_GUI_QRC}
      ${NETWORKEVOLVE_RC}
    )
	
	TARGET_LINK_LIBRARIES( NetworkEvolverGUI
		${QT_LIBRARIES}
	)

    ENDIF( WIN32 )
  ENDIF( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )

ENDIF( APPLE AND BUILD_BUNDLE )

#---------------------------------------------------------------------
#  INSTALL NetworkEvolverGUI executable and required files and folders
#---------------------------------------------------------------------

INSTALL(TARGETS NetworkEvolverGUI 
   RUNTIME DESTINATION bin/../)
   
#INSTALL(FILES networkevolve.qss DESTINATION bin/../)