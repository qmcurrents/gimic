find_package(Eigen3)
set(Boost_USE_MULTITHREADED OFF)
find_package(Boost REQUIRED serialization)

if (ENABLE_MPI)
        find_package(MPI)
        if (MPI_FOUND)
                find_package(Boost REQUIRED mpi)
                set(HAVE_MPI TRUE)
                include_directories(${MPI_INCLUDE_PATH})
        endif()
endif ()


include_directories(
    ${PROJECT_SOURCE_DIR}/src/libmwrepr/src
    ${PROJECT_SOURCE_DIR}/src/libgimic
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    )
link_directories(
    ${PROJECT_BINARY_DIR}/src/libmwrepr/lib
    ${Boost_LIBRARY_DIRS}
    )
add_executable(gimlet.bin 
    gimlet.cpp
    MREnv.cpp
    MREnv.h
    GimicFunc.h
    CurrentFunc.cpp
    CurrentFunc.h
    ModulusFunc.cpp
    ModulusFunc.h
    )
target_link_libraries(gimlet.bin 
    gimic2
    mwrepr 
    ${BLAS_LIBRARIES}
    ${Boost_LIBRARIES})

# This is the "developer" driver
set (INST_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR})
set (GIMIC_EXECUTABLE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimlet.bin)

configure_file(${PROJECT_SOURCE_DIR}/src/gimic/gimic.in
    ${PROJECT_BINARY_DIR}/gimlet)
execute_process(COMMAND
    chmod 755 ${PROJECT_BINARY_DIR}/gimlet OUTPUT_QUIET)

# This is the "installed" driver
set (INST_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib/${PYTHON_INST_DIR})
configure_file(${PROJECT_SOURCE_DIR}/src/gimic/gimic.in
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimlet)
execute_process(COMMAND
    chmod 755
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimlet
    OUTPUT_QUIET)

install(TARGETS gimlet.bin
    RUNTIME DESTINATION bin
    CONFIGURATIONS Release
    )
install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gimlet 
    DESTINATION bin
    CONFIGURATIONS Release
    )
