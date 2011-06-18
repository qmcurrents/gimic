set(MATH_FOUND  FALSE)
set(MATH_LANG   "Fortran")

# user sets
set(USERDEFINED_LIBS
	"${USERDEFINED_LIBS}"
	CACHE STRING
	"User set math libraries"
	FORCE
	)
if(NOT "${USERDEFINED_LIBS}" STREQUAL "")
	set(MATH_LIBS
		"${USERDEFINED_LIBS}"
		CACHE STRING
		"User set math libraries"
		FORCE
		)
	message("-- User set math libraries: ${MATH_LIBS}")
	set(MATH_FOUND TRUE)
endif()

# try to find the best library using environment variables
if(NOT MATH_FOUND)
	find_package(BLAS)
	find_package(LAPACK)
	if(BLAS_FOUND AND LAPACK_FOUND)
		set(MATH_LIBS
			${BLAS_LIBRARIES}
			${LAPACK_LIBRARIES}
			)
		set(MATH_FOUND TRUE)
	endif()
endif()

if(MATH_FOUND)
	set(LIBS
		${LIBS}
		${MATH_LIBS}
		)
else()
	message("-- No external math libraries found")
endif()

