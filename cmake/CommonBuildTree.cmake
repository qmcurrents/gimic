
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
        ${CMAKE_BINARY_DIR}/bin
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
        ${CMAKE_BINARY_DIR}/lib
    )

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
        ${CMAKE_BINARY_DIR}/lib
    )

if (CMAKE_Fortran_COMPILER_WORKS)
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${PROJECT_BINARY_DIR}/modules
        )
    include_directories(${PROJECT_BINARY_DIR}/modules)
endif()

link_directories(${CMAKE_BINARY_DIR}/lib)
include_directories(${PROJECT_BINARY_DIR})
