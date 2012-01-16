include_directories(
    ${PROJECT_BINARY_DIR}/src/libgetkw/Fortran
    )
link_directories(
    ${PROJECT_BINARY_DIR}/external/lib
    )
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
target_link_libraries(gimic.bin 
    gimic2
    fgetkw
    )

# This is the "developer" driver
set (PYTHON_INST_LIBDIR ${CMAKE_SOURCE_DIR}/src/libgetkw/Python) 
set (GIMIC_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic.bin)
if (ENABLE_MRCPP)
    set (GIMLET_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimlet.bin)
else()
    set (GIMLET_EXECUTABLE None)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${PROJECT_BINARY_DIR}/gimic)
execute_process(COMMAND
    chmod 755 ${PROJECT_BINARY_DIR}/gimic OUTPUT_QUIET)

# This is the "installed" driver
set (PYTHON_INST_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib/${PYTHON_INST_DIR})
set (GIMIC_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/gimic.bin)
if (ENABLE_MRCPP)
    set (GIMLET_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/gimlet.bin)
else()
    set (GIMLET_EXECUTABLE None)
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic)
execute_process(COMMAND
    chmod 755
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic
    OUTPUT_QUIET)

install(PROGRAMS ${CMAKE_SOURCE_DIR}/libgetkw/Python/getkw.py 
    DESTINATION lib/${PYTHON_INST_DIR}
    CONFIGURATIONS Release
    )

install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic 
    DESTINATION bin
    CONFIGURATIONS Release
    )

install(TARGETS gimic.bin
    DESTINATION bin
    CONFIGURATIONS Release
    )

if (ENABLE_CFOUR)
    link_directories(${CFOUR_LIBRARY_DIR})
    add_executable(xcpdens 
        xcpdens.f90
        )
    target_link_libraries(xcpdens lib)

    install(PROGRAMS xcpdens
        DESTINATION bin
        CONFIGURATIONS Release
        )
endif()

