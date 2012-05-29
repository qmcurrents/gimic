include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/pygimic)

set(sources
    globals.f90
    getkw.f90  
    basis.f90
    bfeval.f90
    cao2sao.f90
    caos.f90
    dens.f90
    divj.f90
    edens.f90
    factorial.f90
    gaussint.f90
    gtodefs.f90
    intgrl.f90
    jtensor.f90
    kinds.f90
    lip.f90
    reorder.f90
    settings.f90
    teletype.f90
    tensor.f90
    timer.f90
    gausspoints.f90
    gausspoints.c
    gausspoints.h
    gimic_interface.f90
    gimic_interface.h
    GimicInterface.cpp
    GimicInterface.h
    )


add_library(gimic2 SHARED
    ${sources}
    )

if (CYTHON_FOUND)
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/gimic.pxd
        ${CMAKE_CURRENT_SOURCE_DIR}/gimic.pyx
        PROPERTIES CYTHON_IS_CXX TRUE 
        )

    cython_add_module(gimic 
        gimic.pxd 
        gimic.pyx 
        )
    target_link_libraries(gimic 
        gimic2
        )
    cython_add_module(gengauss
        gengauss.pxd 
        gengauss.pyx 
        )
    target_link_libraries(gengauss
        gimic2
        )

    install(TARGETS gimic gengauss
        ARCHIVE DESTINATION ${PYTHON_INST_DIR}
        LIBRARY DESTINATION ${PYTHON_INST_DIR}
        )
endif()

install(TARGETS gimic2 
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    )

install(FILES GimicInterface.h gimic_interface.h gimic.pxd
    DESTINATION include
    )

