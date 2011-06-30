cmake_minimum_required(VERSION 2.6)
project(GIMIC)
enable_language(Fortran)

set(GIMIC_VERSION_MAJOR 1)
set(GIMIC_VERSION_MINOR 2)
set(GIMIC_VERSION_PATCH 3)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    )

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
		Release
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(    NOT cmake_build_type_tolower STREQUAL "debug"
		AND NOT cmake_build_type_tolower STREQUAL "release"
		AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
	message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\".
	Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()

# uncomment this to not rebuild if rules (compiler flags) change
# this has advantages and disadvantages
# set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization" OFF)
#option(ENABLE_CUDA           "Enable CUDA GPU acceleration"  OFF)
option(ENABLE_CFOUR          "Enable the CFOUR interface" OFF)

include(Archs)
include(Compilers)

set(LIBS)

include(Math)
if (MATH_FOUND)
	add_definitions(-DHAVE_BLAS)
endif()

if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
		add_definitions(-DHAVE_MPI)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})
    endif()
endif()

if(ENABLE_OMP)
    find_package(OpenMP)
    if(OpenMP_FOUND)
		add_definitions(-DHAVE_OMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_COMPILE_FLAGS}")
        include_directories(${OpenMP_INCLUDE_PATH})
    endif()
endif()

if (ENABLE_CUDA)
	find_package(CUDA)
endif ()

find_package(PythonInterp REQUIRED)

include_directories(
    ${PROJECT_BINARY_DIR}
    ${PROJECT_BINARY_DIR}/modules
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

# This is the "developer" driver
set (INST_LIBDIR ${CMAKE_SOURCE_DIR})
set (GIMIC_EXECUTABLE ${CMAKE_BINARY_DIR}/gimic.bin)
configure_file(${CMAKE_SOURCE_DIR}/gimic.in
	${CMAKE_BINARY_DIR}/gimic)
execute_process(COMMAND
	chmod 755 ${CMAKE_BINARY_DIR}/gimic OUTPUT_QUIET)

# This is the "installed" driver
set (INST_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib/gimic)
configure_file(${CMAKE_SOURCE_DIR}/gimic.in
	${CMAKE_BINARY_DIR}/bin/gimic)
execute_process(COMMAND
	chmod 755
	${CMAKE_BINARY_DIR}/bin/gimic
	OUTPUT_QUIET)

configure_file(${CMAKE_SOURCE_DIR}/config.h.in
	${CMAKE_BINARY_DIR}/config.h)

add_executable(gimic.bin
	basis.f90
	bfeval.f90
	cao2sao.f90
	caos.f90
	cubeplot.f90
	d2fdrdb.f90
	dbop.f90
	dens.f90
	dfdb.f90
	dfdr.f90
	divj.f90
	edens.f90
	factorial.f90
	gaussint.f90
	getkw.f90
	gimic.f90
	globals.f90
	gopenmol.f90
	grid.f90
	gtodefs.f90
	integral.f90
	intgrl.f90
	jfield.f90
	jtensor.f90
	kinds.f90
	lip.f90
	magnet.f90
	parallel.f90
	reorder.f90
	teletype.f90
	tensor.f90
	timer.f90
	)

install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/gimic 
	DESTINATION bin
	CONFIGURATIONS Release
	)

install(TARGETS gimic.bin
	DESTINATION bin
	CONFIGURATIONS Release
	)

if (ENABLE_CFOUR)
	link_directories(${CFOUR_LIBRARY_DIR})
	add_executable(xcpdens 
		xcpdens.f90
		)
	target_link_libraries(xcpdens lib)

	install(PROGRAMS xcpdens
		DESTINATION bin
		CONFIGURATIONS Release
		)
endif()

add_executable(gomdiff 
    gomdiff.f90
    )

# Packaging. Currently tgz, but zip, rmp, deb, dmg, nsis, etc. possible.
include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
	"GIMIC -- Gauge-Including Magnetically Induced Currents")
set(CPACK_PACKAGE_VENDOR "njj")
set(CPACK_PACKAGE_DESCRIPTION_FILE
	"${CMAKE_CURRENT_SOURCE_DIR}/INSTALL")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${GIMIC_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${GIMIC_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${GIMIC_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY
	"gimic/${GIMIC_VERSION_MAJOR}.${GIMIC_VERSION_MINOR}.${GIMIC_VERSION_PATCH}")
set(CPACK_STRIP_FILES FALSE)
set(CPACK_SOURCE_STRIP_FILES FALSE)
set(CPACK_PACKAGE_EXECUTABLES "gimic" "GIMIC driver executable")
set(CPACK_GENERATOR TGZ)
set(CPACK_SOURCE_GENERATOR TGZ)
set(CPACK_SOURCE_IGNORE_FILES "/CVS/;/.svn/;/.bzr/;/.hg/;/.git.*/;.swp$;.#;/#")
set(CPACK_SOURCE_IGNORE_FILES
	"/(Build|build|BUILD)-?.*/;/(Release|release|RELEASE)-?.*/;/(Debug|debug|DEBUG)-?.*/;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)
