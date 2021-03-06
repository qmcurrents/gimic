
add_executable(gimic.bin
    gimic.F90
    getkw.f90
    integral.f90
    vtkplot.f90
    magnet.f90
    grid.f90
    jfield.f90
    parallel.F90
    )
target_link_libraries(gimic.bin gimic2)

if(ENABLE_BLAS)
    target_link_libraries(gimic.bin ${BLAS_LIBRARIES})
endif()

if(ENABLE_LAPACK)
    target_link_libraries(gimic.bin ${LAPACK_LIBRARIES})
endif()

if(ENABLE_MKL_FLAG)
    set_target_properties(gimic.bin PROPERTIES LINK_FLAGS "-mkl=parallel")
endif()

# without this line, linking with Intel stack fails
set_target_properties(gimic.bin PROPERTIES LINKER_LANGUAGE Fortran)

configure_file(gimic.py.in gimic.py)

install(FILES __init__.py ${CMAKE_CURRENT_BINARY_DIR}/gimic.py getkw.py
    DESTINATION ${PYTHON_SITE_INSTALL_DIR}/fgimic
    )

install(TARGETS gimic.bin
    DESTINATION bin
    )

if (ENABLE_CFOUR)
    find_package(CFour)
    if (CFOUR_FOUND)
        link_directories(${CFOUR_LIBRARY_DIR})
        add_executable(xcpdens
            xcpdens.F90
            )
        target_link_libraries(xcpdens libr)

        install(TARGETS xcpdens
            DESTINATION bin
            )
    endif()
endif()

