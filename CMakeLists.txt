cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

project(GIMIC)
enable_language(Fortran)

set (PROJECT_VERSION_MAJOR 2)
set (PROJECT_VERSION_MINOR 0)
set (PROJECT_VERSION_PATCH 0)

option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_OPENMP         "Enable OpenMP parallelization" OFF)
option(ENABLE_CFOUR          "Enable the CFOUR interface" OFF)
option(ENABLE_TESTS          "Enable the test suite" ON)
option(ENABLE_MRCPP          "Enable the MRCPP wavelet librbary" OFF)

include(GNUInstallDirs)
include(ConfigProjectVersion)
include(ConfigSafeGuards)
include(UseGlobalBuildTree)
include(ConfigGitRevision)
include(ConfigureScript)
include(ConfigCompilerFlags)

find_package(Python REQUIRED)
set (PYTHON_INST_DIR 
    lib/python${PYTHON_VERSION}/site-packages
    )
include(UseCython)

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

if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        add_definitions(-DHAVE_OPENMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
        include_directories(${OpenMP_INCLUDE_PATH})
    endif()
endif()

configure_file (
    ${PROJECT_SOURCE_DIR}/config.h.in
    ${PROJECT_BINARY_DIR}/config.h
    )

add_subdirectory(src)
add_subdirectory(tools)

if (ENABLE_TESTS)
    set (PROJECT_REPOSITORY "git@repo.ctcc.no:gimic.git")
    include(ConfigTesting)
    add_subdirectory(test)
endif()

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "GIMIC - Gauge-Including Magnetically Induced Currents program")
set(CPACK_PACKAGE_VENDOR "University of Helsinki")
set(CPACK_PACKAGE_CONTACT "Jonas Juselius <jonas.juselius@uit.no>")
set(CPACK_PACKAGE_EXECUTABLES "gimic" "The GIMIC program")
set(CPACK_PACKAGE_URL "http://repo.ctcc.no/projects/gimic")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/INSTALL")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
#set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_CURRENT_SOURCE_DIR}/INFO")
set(CPACK_GENERATOR TGZ DEB RPM)
include(ConfigPackaging) 

