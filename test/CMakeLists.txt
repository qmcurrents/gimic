include_directories (${PROJECT_SOURCE_DIR}/src)

if (MPI_FOUND)
    set(MPI_NUMPROC 3)
    set(MPIEXEC_PREFLAGS "")
    set (MPIRUN ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUMPROC}
        ${MPIEXEC_PREFLAGS})
endif()

# Configure driver
set (GIMIC_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic.bin)
set (LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set (PYTHON_INSTDIR ${PROJECT_SOURCE_DIR}/src/gimic)

configure_file(${CMAKE_SOURCE_DIR}/src/gimic.in
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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/PyGimicTest.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/PyGimicTest.py @ONLY
    )
execute_process(
    COMMAND chmod 755 ${CMAKE_CURRENT_BINARY_DIR}/test.py
    OUTPUT_QUIET
    )

# Add tests
add_test(NAME Vectors
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test.py vectors
    )

add_test(NAME BondIntegral
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test.py bondint
    )

foreach(_name
    benzene/integration
    benzene/2d
    benzene/3d
    c4h4/integration
    )
    add_test(
        NAME ${_name}
        COMMAND ${PROJECT_SOURCE_DIR}/test/${_name}/test
                --binary-dir=${PROJECT_BINARY_DIR}/bin
                --work-dir=${PROJECT_BINARY_DIR}/test/${_name}
                --verbose
        )
endforeach()
