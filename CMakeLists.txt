cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/custom ${CMAKE_MODULE_PATH})

project(GIMIC)
enable_language(Fortran)

option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_OPENMP         "Enable OpenMP parallelization" OFF)
option(ENABLE_CFOUR          "Enable the CFOUR interface" OFF)

unset (CMAKE_LIBRARY_ARCHITECTURE)
include(GNUInstallDirs)
include(ConfigProjectVersion)
include(ConfigSafeGuards)
include(ConfigGitRevision)
include(ConfigCompilerFlags)

find_package(Python REQUIRED)
set(PYTHON_SITE_INSTALL_DIR 
    lib/python${PYTHON_VERSION}/site-packages/gimic)
include(UseCython)

set(LIBS)

if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DHAVE_MPI)
        include_directories(${MPI_INCLUDE_PATH})
    endif()
endif()

if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        add_definitions(-DHAVE_OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
        include_directories(${OpenMP_INCLUDE_PATH})
    endif()
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

configure_file (
    ${PROJECT_SOURCE_DIR}/config.h.in
    ${PROJECT_BINARY_DIR}/config.h
    )

include(ConfigTesting)

add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(test)

include(ConfigPackaging) 
