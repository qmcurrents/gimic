
if (NOT DEFINED CFOUR_ROOT)
    if (EXISTS $ENV{CFOUR_ROOT})
        set(CFOUR_ROOT $ENV{CFOUR_ROOT})
    else()
        set(CFOUR_ROOT /usr/local/cfour)
    endif()
endif()

configure_file(MOL2mol.sh.in MOL2mol.sh)

configure_file(qgimic.sh.in qgimic.sh)

configure_file(quesys.sh.in quesys.sh) 
configure_file(turbo2gimic.py.in turbo2gimic.py)

install(PROGRAMS 
    ${CMAKE_CURRENT_BINARY_DIR}/MOL2mol.sh
    ${CMAKE_CURRENT_BINARY_DIR}/qgimic.sh
    ${CMAKE_CURRENT_BINARY_DIR}/quesys.sh
    ${CMAKE_CURRENT_BINARY_DIR}/turbo2gimic.py
    DESTINATION bin
    )

install(DIRECTORY QCTools 
    DESTINATION ${PYTHON_SITE_INSTALL_DIR}
    FILES_MATCHING PATTERN *.py
    )
