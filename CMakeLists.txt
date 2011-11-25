cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project(GIMIC)
enable_language(Fortran)

option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization" OFF)
option(ENABLE_CFOUR          "Enable the CFOUR interface" OFF)

include(ConfigProjectVersion)
include(ConfigSafeGuards)
include(ConfigFortranCompilerFlags)
include(ConfigCCompilerFlags)

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
    if(OpenMP_FOUND)
		add_definitions(-DHAVE_OMP)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_COMPILE_FLAGS}")
        include_directories(${OpenMP_INCLUDE_PATH})
    endif()
endif()

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

if (ENABLE_TESTS)
	set (PROJECT_REPOSITORY "git@repo.ctcc.no:gimic.git")
	include(ConfigTesting)
	add_subdirectory(tests)
endif()

include(ConfigGitRevision)

if (EXISTS ${CMAKE_SOURCE_DIR}/cmake/ConfigPackaging.cmake)
	include(ConfigPackaging) 
endif()

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

configure_file (
	${PROJECT_SOURCE_DIR}/config.h.in
	${PROJECT_BINARY_DIR}/config.h
	)

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

