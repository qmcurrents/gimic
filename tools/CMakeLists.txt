include(ConfigureScript)

if (NOT DEFINED CFOUR_ROOT)
    if (EXISTS $ENV{CFOUR_ROOT})
        set(CFOUR_ROOT $ENV{CFOUR_ROOT})
    else()
        set(CFOUR_ROOT /usr/local/cfour)
    endif()
endif()

configure_script(${CMAKE_CURRENT_SOURCE_DIR}/MOL2mol.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/MOL2mol.sh
    )

configure_script(${CMAKE_CURRENT_SOURCE_DIR}/qgimic.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/qgimic.sh
    )

configure_script(${CMAKE_CURRENT_SOURCE_DIR}/quesys.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/quesys.sh
    )

configure_script(${CMAKE_CURRENT_SOURCE_DIR}/turbo2gimic.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/turbo2gimic.py
    )

install(PROGRAMS 
    ${CMAKE_CURRENT_BINARY_DIR}/MOL2mol.sh
    ${CMAKE_CURRENT_BINARY_DIR}/qgimic.sh
    ${CMAKE_CURRENT_BINARY_DIR}/quesys.sh
    ${CMAKE_CURRENT_BINARY_DIR}/turbo2gimic.py
    DESTINATION bin
    )

install(FILES 
    Atom.py
    Const.py
    Elements.py
    QCTools.py
    DESTINATION ${PYTHON_INST_DIR}
    )
