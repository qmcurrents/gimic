# Configure driver
set (GIMIC_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic.bin)
set (LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set (PYTHON_INSTDIR ${PROJECT_SOURCE_DIR}/src/gimic)

configure_file(${CMAKE_SOURCE_DIR}/src/gimic.in
    ${CMAKE_CURRENT_BINARY_DIR}/gimic @ONLY
    )
execute_process(COMMAND
    chmod 755 ${CMAKE_CURRENT_BINARY_DIR}/gimic OUTPUT_QUIET)

foreach(_name
    benzene/integration
    benzene/bond-integration
    benzene/vectors
    benzene/2d
    benzene/3d
    c4h4/integration
    )
    add_test(
        NAME ${_name}
        COMMAND ./test
                --binary-dir=${PROJECT_BINARY_DIR}/bin
                --work-dir=${PROJECT_BINARY_DIR}/test/${_name}
                --verbose
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/${_name}
        )
endforeach()
