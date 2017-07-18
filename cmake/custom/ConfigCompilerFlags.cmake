if (CMAKE_C_COMPILER_WORKS)
    include(${PROJECT_SOURCE_DIR}/cmake/custom/compilers/CFlags.cmake)
endif()

if (CMAKE_CXX_COMPILER_WORKS)
    include(${PROJECT_SOURCE_DIR}/cmake/custom/compilers/CXXFlags.cmake)
endif()

if (CMAKE_Fortran_COMPILER_WORKS)
    include(${PROJECT_SOURCE_DIR}/cmake/custom/compilers/FortranFlags.cmake)
endif()
