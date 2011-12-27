find_package(Eigen3)

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
add_executable(gimlet.bin gimlet.cpp)
target_link_libraries(gimlet.bin 
    gimic2
    mwrepr 
    ${BLAS_LIBRARIES}
    ${Boost_LIBRARIES})

install(TARGETS gimlet.bin
    RUNTIME DESTINATION bin
    CONFIGURATIONS Release
    )
