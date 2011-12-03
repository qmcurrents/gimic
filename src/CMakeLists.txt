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

add_library(gimic STATIC
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

add_executable(gimic.bin gimic.f90)
target_link_libraries(gimic.bin gimic)

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


