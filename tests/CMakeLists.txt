#include(LocalBuildTree)
include_directories (${PROJECT_SOURCE_DIR}/src)

# Configure driver
set (GIMIC_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic.bin)
set (INST_LIBDIR ${PROJECT_SOURCE_DIR}/src/gimic)

configure_file(${CMAKE_SOURCE_DIR}/src/gimic/gimic.in
    ${CMAKE_CURRENT_BINARY_DIR}/gimic @ONLY
    )
execute_process(COMMAND 
    chmod 755 ${CMAKE_CURRENT_BINARY_DIR}/gimic OUTPUT_QUIET) 

# Configure input
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/vectors.inp.in
    ${CMAKE_CURRENT_BINARY_DIR}/vectors.inp @ONLY
    )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bondint.inp.in
    ${CMAKE_CURRENT_BINARY_DIR}/bondint.inp @ONLY
    )

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

