include(LocalBuildTree)
include_directories (${PROJECT_SOURCE_DIR}/src)

add_custom_target(tests ALL 
    DEPENDS test.bin
    )

# Configure test.bin
add_executable(test.bin ${PROJECT_SOURCE_DIR}/src/gimic.f90)
target_link_libraries(test.bin gimic)

# Configure driver
set (GIMIC_EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/test.bin)
set (INST_LIBDIR ${PROJECT_SOURCE_DIR}/src)

configure_file(${CMAKE_SOURCE_DIR}/src/gimic.in
    ${CMAKE_CURRENT_BINARY_DIR}/gimic-test @ONLY
    )
execute_process(COMMAND 
    chmod 755 ${CMAKE_CURRENT_BINARY_DIR}/gimic-test OUTPUT_QUIET) 

# Configure input
if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/gimic.inp) 
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.inp.in
        ${CMAKE_CURRENT_BINARY_DIR}/gimic.inp @ONLY
        )
endif()

# Configure script
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/test.py @ONLY
    )
execute_process(
    COMMAND chmod 755 ${CMAKE_CURRENT_BINARY_DIR}/test.py 
    OUTPUT_QUIET
    ) 

# Add tests
add_test(NAME Vectors
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test.py vectors
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

add_test(NAME BondIntegral
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test.py bond_integral
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

