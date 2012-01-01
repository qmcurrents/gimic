#include(LocalBuildTree)
include_directories (${PROJECT_SOURCE_DIR}/src)

if (MPI_FOUND)
    set(MPI_NUMPROC 3)
    set(MPIEXEC_PREFLAGS "")
    set (MPIRUN ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUMPROC}
        ${MPIEXEC_PREFLAGS})
endif()

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
    #WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test.py vectors
    )

add_test(NAME BondIntegral
    #WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test.py bond_integral
    )

add_test(NAME VectorsMPI
    #WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${MPIRUN} ${CMAKE_CURRENT_BINARY_DIR}/test.py vectors
    )

add_test(NAME BondIntegralMPI
    #WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${MPIRUN} ${CMAKE_CURRENT_BINARY_DIR}/test.py bond_integral
    )
