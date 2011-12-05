add_subdirectory(libgimic)

add_executable(gimic.bin 
    edens.f90
    gimic.f90
    gopenmol.f90  
    integral.f90  
    magnet.f90
    cubeplot.f90   
    grid.f90	    
    jfield.f90	  
    divj_field.f90
    )
target_link_libraries(gimic.bin gimic)

# This is the "developer" driver
set (INST_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR})
set (GIMIC_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic.bin)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${PROJECT_BINARY_DIR}/gimic)
execute_process(COMMAND
    chmod 755 ${PROJECT_BINARY_DIR}/gimic OUTPUT_QUIET)

# This is the "installed" driver
set (INST_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib/gimic)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic)
execute_process(COMMAND
    chmod 755
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic
    OUTPUT_QUIET)

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

add_executable(gomdiff 
    gomdiff.f90
    )

