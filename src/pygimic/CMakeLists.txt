include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if (CYTHON_FOUND)
    cython_add_module(connector
        connector.pxd
        connector.pyx
	)
    target_link_libraries( connector ${PYTHON_LIBRARIES} )
install(TARGETS connector
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION ${PYTHON_SITE_INSTALL_DIR}
    )

endif()

