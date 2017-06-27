# hack for Travis because we have a strange dependency on numpy header files
# remove this once we get rid of Cython
if(DEFINED ENV{TRAVIS_TEST})
    include_directories(/home/travis/virtualenv/python2.7.9/lib/python2.7/site-packages/numpy/core/include/)
endif()

include_directories(${PROJECT_SOURCE_DIR}/venv/lib/python2.7/site-packages/numpy/core/include/)

include_directories(${CMAKE_BINARY_DIR})
include(ExternalProject)
add_subdirectory(libgimic)
add_subdirectory(pygimic)
add_subdirectory(fgimic)

configure_file(gimic.in gimic)

if (ENABLE_MRCPP)
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

    add_subdirectory(gimlet)
    find_package(Git)

    ExternalProject_Add(libmwrepr
        DOWNLOAD_COMMAND git clone git@repo.ctcc.no:mrcpp/libmwrepr.git
        DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        UPDATE_COMMAND git pull
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

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/gimic
    DESTINATION bin
    )

install(FILES fgimic/getkw.py
    DESTINATION ${PYTHON_SITE_INSTALL_DIR}
    )

install(DIRECTORY pygimic
    DESTINATION ${PYTHON_SITE_INSTALL_DIR}
    FILES_MATCHING PATTERN "*.py"
    )

install(DIRECTORY london
    DESTINATION ${PYTHON_SITE_INSTALL_DIR}
    FILES_MATCHING PATTERN "*.py"
    )
