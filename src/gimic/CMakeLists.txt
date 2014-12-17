add_executable(gimic.bin 
    gimic.f90
    integral.f90  
    cubeplot.f90   
    vtkplot.f90   
    magnet.f90
    grid.f90	    
    jfield.f90	  
    divj_field.f90
    parallel.f90
    edens_field.f90
    )
target_link_libraries(gimic.bin gimic2)

# This is the "developer" driver
set (_PYTHON_INSTDIR ${PYTHON_INSTDIR})
set (LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set (PYTHON_INSTDIR ${CMAKE_CURRENT_SOURCE_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${PROJECT_BINARY_DIR}/gimic)
execute_process(COMMAND
    chmod 755 ${PROJECT_BINARY_DIR}/gimic OUTPUT_QUIET)

# This is the "installed" driver
set (LIBDIR ${CMAKE_INSTALL_FULL_LIBDIR})
set (PYTHON_INSTDIR ${_PYTHON_INSTDIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic)
execute_process(COMMAND
    chmod 755
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic
    OUTPUT_QUIET)

install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/getkw.py 
    DESTINATION ${PYTHON_INSTDIR}
    )

install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic 
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

install(TARGETS gimic.bin
    DESTINATION ${CMAKE_INSTALL_BINDIR}
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
            DESTINATION ${CMAKE_INSTALL_BINDIR}
            )
    endif()
endif()

unset (_PYTHON_FULL_INSTDIR)
