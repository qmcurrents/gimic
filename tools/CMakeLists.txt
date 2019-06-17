
if (NOT DEFINED CFOUR_ROOT)
    if (EXISTS $ENV{CFOUR_ROOT})
        set(CFOUR_ROOT $ENV{CFOUR_ROOT})
    else()
        set(CFOUR_ROOT /usr/local/cfour)
    endif()
endif()

configure_file(MOL2mol.sh.in MOL2mol.sh)

configure_file(qgimic.sh.in qgimic.sh)
configure_file(grd2node.sh.in grd2node.sh)

configure_file(quesys.sh.in quesys.sh) 
configure_file(turbo2gimic.py.in turbo2gimic.py)
configure_file(sortorbitals.py.in sortorbitals.py)
configure_file(unstructured_grid_gen.py.in unstructured_grid_gen.py)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/PyGimicTest.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/PyGimicTest.py @ONLY
    )

install(PROGRAMS 
    ${CMAKE_CURRENT_BINARY_DIR}/MOL2mol.sh
    ${CMAKE_CURRENT_BINARY_DIR}/qgimic.sh
    ${CMAKE_CURRENT_BINARY_DIR}/grd2node.sh
    ${CMAKE_CURRENT_BINARY_DIR}/quesys.sh
    ${CMAKE_CURRENT_BINARY_DIR}/turbo2gimic.py
    ${CMAKE_CURRENT_BINARY_DIR}/sortorbitals.py
    ${CMAKE_CURRENT_BINARY_DIR}/unstructured_grid_gen.py
    DESTINATION bin
    )

install(DIRECTORY QCTools 
#   DESTINATION ${PYTHON_SITE_INSTALL_DIR}
    DESTINATION bin 
    FILES_MATCHING PATTERN *.py
    )
