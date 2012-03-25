add_executable(gimic.bin 
    gimic.f90
    integral.f90  
    cubeplot.f90   
    magnet.f90
    grid.f90	    
    jfield.f90	  
    divj_field.f90
    parallel.f90
    edens_field.f90
    )
target_link_libraries(gimic.bin gimic2)

# This is the "developer" driver
set (INST_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR})
set (GIMIC_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic.bin)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${PROJECT_BINARY_DIR}/gimic)
execute_process(COMMAND
    chmod 755 ${PROJECT_BINARY_DIR}/gimic OUTPUT_QUIET)

# This is the "installed" driver
set (INST_LIBDIR ${CMAKE_INSTALL_PREFIX}/${PYTHON_INST_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic)
execute_process(COMMAND
    chmod 755
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic
    OUTPUT_QUIET)

install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/getkw.py 
    DESTINATION ${PYTHON_INST_DIR}
    )

install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic 
    DESTINATION bin
    )

install(TARGETS gimic.bin
    DESTINATION bin
    )

if (ENABLE_CFOUR)
    find_package(CFour)
    if (CFOUR_FOUND)
        link_directories(${CFOUR_LIBRARY_DIR})
        add_executable(xcpdens 
            xcpdens.f90
            )
        target_link_libraries(xcpdens libr)

        install(TARGETS xcpdens
            DESTINATION bin
            )
    endif()
endif()

