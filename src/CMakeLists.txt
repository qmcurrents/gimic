# hack for Travis because we have a strange dependency on numpy header files
# remove this once we get rid of Cython
if(DEFINED ENV{TRAVIS_TEST})
    include_directories(/home/travis/virtualenv/python2.7.9/lib/python2.7/site-packages/numpy/core/include/)
endif()

include_directories(${PROJECT_SOURCE_DIR}/venv/lib/python2.7/site-packages/numpy/core/include/)

include_directories(${CMAKE_BINARY_DIR})
include(ExternalProject)
add_subdirectory(libgimic)
add_subdirectory(pygimic)
add_subdirectory(fgimic)

configure_file(gimic.in gimic)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/gimic
    DESTINATION bin
    )

install(FILES fgimic/getkw.py
    DESTINATION ${PYTHON_SITE_INSTALL_DIR}
    )

install(DIRECTORY pygimic
    DESTINATION ${PYTHON_SITE_INSTALL_DIR}
    FILES_MATCHING PATTERN "*.py"
    )

install(DIRECTORY london
    DESTINATION ${PYTHON_SITE_INSTALL_DIR}
    FILES_MATCHING PATTERN "*.py"
    )
