ADD_EXECUTABLE(
  SBWGA 
  sbw_ga_optim.cpp
)

TARGET_LINK_LIBRARIES( SBWGA sbml_sim galib SBWCore)


FIND_PACKAGE( Qt4 REQUIRED )
IF( QT4_FOUND )
    INCLUDE_DIRECTORIES( BEFORE ${QT_INCLUDES} )
    IF( QT_USE_FILE )
	  SET( QT_USE_QTCORE TRUE )
	  SET( QT_USE_QTGUI TRUE )
	  INCLUDE( ${QT_USE_FILE} )
    ELSE( QT_USE_FILE )
	  SET( QT_LIBRARIES ${QT_QT_LIBRARY} )
    ENDIF( QT_USE_FILE )
	
	QT4_WRAP_CPP( QWT_MOC ${QWT_HDRS} )

	FILE( GLOB HEADERS plotWidget.h sbwclient.h randomize_widget.h )
	FILE( GLOB SOURCES plotWidget.cpp sbwclient.cpp randomize_widget.cpp )	
    QT4_WRAP_CPP( HEADERS_MOC ${HEADERS} )

	IF( APPLE AND BUILD_BUNDLE )
	  ADD_EXECUTABLE(
	    JSimOptim
		MACOSX_BUNDLE
		${SOURCES}
		${HEADERS_MOC}
	  )
		
	  TARGET_LINK_LIBRARIES( JSimOptim  ${QT_LIBRARIES} ${SBWLIB} )
		
	  ADD_CUSTOM_COMMAND( TARGET JSimOptim POST_BUILD
		  COMMAND mkdir ARGS -p
			${EXECUTABLE_OUTPUT_PATH}/JSimOptim.app/Contents/Resources
			${EXECUTABLE_OUTPUT_PATH}/JSimOptim.app/Contents/Frameworks
			${EXECUTABLE_OUTPUT_PATH}/JSimOptim.app/Contents/MacOS
		  COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtCore.framework
			${EXECUTABLE_OUTPUT_PATH}/JSimOptim.app/Contents/Frameworks/
		  COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtGui.framework
			${EXECUTABLE_OUTPUT_PATH}/JSimOptim.app/Contents/Frameworks/
		  COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtXml.framework
			${EXECUTABLE_OUTPUT_PATH}/JSimOptim.app/Contents/Frameworks/
		  COMMAND cp ARGS -rf ${QT_LIBRARY_DIR}/QtOpenGL.framework
			${EXECUTABLE_OUTPUT_PATH}/JSimOptim.app/Contents/Frameworks/
		  COMMAND cp ARGS ${MACOSX_BUNDLE_ICON_FILE}
			${EXECUTABLE_OUTPUT_PATH}/JSimOptim.app/Contents/Resources
	  )
		
	ELSE( APPLE AND BUILD_BUNDLE )
	  IF( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )
		
		ADD_EXECUTABLE(
		  JSimOptim
		  ${SOURCES}
		  ${HEADERS_MOC}
		  ${QWT_SRC}
	      ${QWT_MOC}
		)
		
		TARGET_LINK_LIBRARIES(
		  JSimOptim
		  ${QT_LIBRARIES}
		  ${SBWLIB}
		)
		
	  ELSE( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )
		IF( WIN32 )
		
		ADD_EXECUTABLE(
		  JSimOptim 
		  ${SOURCES}
		  ${HEADERS_MOC}
		  ${QWT_SRC}
		  ${QWT_MOC}
		)
		
		TARGET_LINK_LIBRARIES( 
		  JSimOptim
		  ${QT_LIBRARIES}
		  ${SBWLIB}
		)
		
		ENDIF( WIN32 )
	  ENDIF( UNIX OR ( APPLE AND NOT BUILD_BUNDLE ) )

	ENDIF( APPLE AND BUILD_BUNDLE )

ENDIF( QT4_FOUND )
