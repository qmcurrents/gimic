include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(gimic2 SHARED
    globals.f90
    getkw.f90  
    basis.f90
    bfeval.f90
    cao2sao.f90
    caos.f90
    d2fdrdb.f90
    dbop.f90
    dens.f90
    dfdb.f90
    dfdr.f90
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
    gimic_cif.f90
    gimic_cif.h
    GimicInterface.cpp
    GimicInterface.h
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
install(TARGETS gimic 
    ARCHIVE DESTINATION lib/${PYTHON_INST_DIR}
    LIBRARY DESTINATION lib
    CONFIGURATIONS Release
    )
endif()

install(TARGETS gimic2 
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    CONFIGURATIONS Release
    )

