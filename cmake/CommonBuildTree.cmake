if (PROJECT_RUNTIME_OUTPUT_DIR)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
        ${PROJECT_RUNTIME_OUTPUT_DIR}
        )
else()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
        ${CMAKE_BINARY_DIR}/bin
        )
endif()

if (PROJECT_ARCHIVE_OUTPUT_DIR)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
        ${PROJECT_ARCHIVE_OUTPUT_DIR}
        )
else()
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
        ${CMAKE_BINARY_DIR}/lib
        )
endif()

if (PROJECT_LIBRARY_OUTPUT_DIR)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
        ${PROJECT_LIBRARY_OUTPUT_DIR}
        )
else()
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
        ${CMAKE_BINARY_DIR}/lib
        )
endif()

if (CMAKE_Fortran_COMPILER_WORKS)
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${PROJECT_BINARY_DIR}/modules
        )
    include_directories(${PROJECT_BINARY_DIR}/modules)
endif()

link_directories(${CMAKE_BINARY_DIR}/lib)
include_directories(${PROJECT_BINARY_DIR})
