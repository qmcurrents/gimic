set (PROJECT_VERSION_MAJOR 0)
set (PROJECT_VERSION_MINOR 1)
set (PROJECT_VERSION_PATCH 0)

set (PROJECT_VERSION 
    ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})

find_package(PythonInterp)
if (PYTHONINTERP-FOUND) 
    configure_file(
        ${CMAKE_SOURCE_DIR}/cmake/bump-version.in 
        ${CMAKE_SOURCE_DIR}/cmake/bump-version
        )
    execute_process(COMMAND chmod 755 
        ${CMAKE_SOURCE_DIR}/cmake/bump-version OUTPUT_QUIET) 
endif()

