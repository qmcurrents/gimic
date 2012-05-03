include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if (CYTHON_FOUND)
    cython_add_module(connector
        connector.pxd
        connector.pyx 
    )
install(TARGETS connector 
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    )

endif()
