cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project(GIMIC)
enable_language(Fortran)

option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization" ON)
option(ENABLE_CFOUR          "Enable the CFOUR interface" OFF)
option(ENABLE_TESTS          "Enable the test suite" ON)

include(ConfigProjectVersion)
include(ConfigSafeGuards)

include(ConfigFortranCompilerFlags)
include(ConfigCCompilerFlags)

find_package(PythonInterp REQUIRED)
include(UseCython)
include(ConfigureScript)
include(FortranNameManglingScheme)

set(LIBS)
if (BLAS_FOUND)
    add_definitions(-DHAVE_BLAS)
endif()

if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DHAVE_MPI)
        include_directories(${MPI_INCLUDE_PATH})
    endif()
endif()

if(ENABLE_OMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        add_definitions(-DHAVE_OMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
        include_directories(${OpenMP_INCLUDE_PATH})
    endif()
endif()

find_package(PythonInterp REQUIRED)

include(CommonBuildTree)
include(ConfigGitRevision)

configure_file (
    ${PROJECT_SOURCE_DIR}/config.h.in
    ${PROJECT_BINARY_DIR}/config.h
    )

add_subdirectory(src)

if (ENABLE_TESTS)
    set (PROJECT_REPOSITORY "git@repo.ctcc.no:gimic.git")
    include(ConfigTesting)
    add_subdirectory(tests)
endif()

if (EXISTS ${CMAKE_SOURCE_DIR}/cmake/ConfigPackaging.cmake)
    include(ConfigPackaging) 
endif()

