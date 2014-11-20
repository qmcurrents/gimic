include(ExternalProject)
add_subdirectory(libgimic)
add_subdirectory(pygimic)
add_subdirectory(fgimic)

# This is the "developer" driver
set (LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set (_PYTHON_INSTDIR ${PYTHON_INSTDIR})
set (PYTHON_INSTDIR ${CMAKE_CURRENT_SOURCE_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${PROJECT_BINARY_DIR}/gimic)
execute_process(COMMAND
    chmod 755 ${PROJECT_BINARY_DIR}/gimic OUTPUT_QUIET)

# This is the "installed" driver
set (PYTHON_INSTDIR ${_PYTHON_INSTDIR})
set (LIBDIR ${CMAKE_INSTALL_FULL_LIBDIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gimic.in
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic)
execute_process(COMMAND
    chmod 755 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic OUTPUT_QUIET)

ExternalProject_Add(libgetkw
    DOWNLOAD_COMMAND git submodule update
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libgetkw
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libgetkw
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/libgetkw
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    #-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    #-DPROJECT_RUNTIME_OUTPUT_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    )

ExternalProject_Add_Step(libgetkw git-update
    COMMAND git submodule init
    COMMAND git submodule update
    ALWAYS 1
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

if (ENABLE_MRCPP)
    add_subdirectory(gimlet)
    find_package(Git)

    ExternalProject_Add(libmwrepr
        DOWNLOAD_COMMAND git submodule update
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libmwrepr
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libmwrepr
        INSTALL_DIR ${CMAKE_BINARY_DIR}/external
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external
        -DCMAKE_BUILD_TYPE=Release
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

install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimic
    DESTINATION ${CMAKE_INSTALL_FULL_BINDIR}
    )

install(DIRECTORY pygimic
    DESTINATION ${PYTHON_INSTDIR}
    FILES_MATCHING PATTERN "*.py"
    )

install(DIRECTORY london
    DESTINATION ${PYTHON_INSTDIR}
    FILES_MATCHING PATTERN "*.py"
    )
