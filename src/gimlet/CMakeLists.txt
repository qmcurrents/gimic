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
    ${PROJECT_SOURCE_DIR}/src/libgetkw/C++
    ${PROJECT_SOURCE_DIR}/src/libgimic
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    )
link_directories(
    ${PROJECT_BINARY_DIR}/external/lib
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
    getkw
    ${BLAS_LIBRARIES}
    ${Boost_LIBRARIES})

install(TARGETS gimlet.bin
    RUNTIME DESTINATION bin
    )
