if (CMAKE_C_COMPILER_WORKS)
    include(${PROJECT_SOURCE_DIR}/cmake/compilers/CFlags.cmake)
endif()

if (CMAKE_CXX_COMPILER_WORKS)
    include(${PROJECT_SOURCE_DIR}/cmake/compilers/CXXFlags.cmake)
endif()

if (CMAKE_Fortran_COMPILER_WORKS)
    include(${PROJECT_SOURCE_DIR}/cmake/compilers/FortranFlags.cmake)
endif()
