add_subdirectory(libgimic)
add_subdirectory(pygimic)
add_subdirectory(fgimic)

# This is the "developer" driver
set (INST_LIBDIR ${CMAKE_BINARY_DIR}/lib)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${PROJECT_BINARY_DIR}/gimic)
execute_process(COMMAND
    chmod 755 ${PROJECT_BINARY_DIR}/gimic OUTPUT_QUIET)


if (ENABLE_MRCPP)
    find_package(Git)
    include(ExternalProject)

    ExternalProject_Add(libmwrepr
        DOWNLOAD_COMMAND git submodule update
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libmwrepr
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libmwrepr
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/libmwrepr
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/libmwrepr 
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        #-DPROJECT_RUNTIME_OUTPUT_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        )

    ExternalProject_Add_Step(libmwrepr git-update
        COMMAND git submodule init
        COMMAND git submodule update
        ALWAYS 1
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
endif()
