set (PROJECT_VERSION_MAJOR 1)
set (PROJECT_VERSION_MINOR 5)
set (PROJECT_VERSION_PATCH 1)

set (PROJECT_VERSION 
    ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})

find_package(PythonInterp)
if (PYTHONINTERP_FOUND) 
    configure_file(
        ${PROJECT_SOURCE_DIR}/cmake/bump-version.in 
        ${PROJECT_SOURCE_DIR}/cmake/bump-version
        )
    execute_process(COMMAND chmod 755 
        ${PROJECT_SOURCE_DIR}/cmake/bump-version OUTPUT_QUIET) 
endif()

