set (PROJECT_VERSION 
    ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})

find_package(PythonInterp)
if (PYTHONINTERP_FOUND)
    configure_file(
        ${CMAKE_SOURCE_DIR}/cmake/bump-version.in
        ${CMAKE_SOURCE_DIR}/cmake/bump-version
        )
    execute_process(COMMAND chmod 755
        ${CMAKE_SOURCE_DIR}/cmake/bump-version OUTPUT_QUIET)
endif()

