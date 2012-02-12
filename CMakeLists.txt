cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(GIMIC)
enable_language(Fortran)

set (PROJECT_VERSION_MAJOR 1)
set (PROJECT_VERSION_MINOR 5)
set (PROJECT_VERSION_PATCH 1)

option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_OPENMP          "Enable OpenMP parallelization" OFF)
option(ENABLE_CFOUR          "Enable the CFOUR interface" OFF)
option(ENABLE_TESTS          "Enable the test suite" ON)
option(ENABLE_MRCPP          "Enable the MRCPP wavelet librbary" OFF)

include(GNUInstallDirs)
include(ConfigProjectVersion)
include(ConfigSafeGuards)
include(ConfigCompilerFlags)

find_package(PythonInterp REQUIRED)
include(UseCython)
include(ConfigureScript)

include(FortranCInterface)
if (FortranCInterface_GLOBAL_FOUND) 
    FortranCInterface_HEADER(fcmangle.h
        MACRO_NAMESPACE "FC_"
        SYMBOL_NAMESPACE "FC_"
        )
else()
    message(WARNING "Failed to detect Fortran name mangling scheme")
endif()

if (PYTHONINTERP_FOUND) 
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} --version
        ERROR_VARIABLE PYTHON_VERSION
        )
    string(REGEX REPLACE "Python *([0-9][.][0-9])[.].*" \\1
        PYTHON_VERSION ${PYTHON_VERSION}
        )
    set (PYTHON_INST_DIR 
        python${PYTHON_VERSION}/site-packages
        )
endif()

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

find_package(PythonInterp REQUIRED)

include(CommonBuildTree)
include(ConfigGitRevision)

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

