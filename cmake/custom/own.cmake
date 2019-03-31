set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})

include(GNUInstallDirs)

include(version)
include(git)

set(PYTHON_SITE_INSTALL_DIR
    lib/python${PYTHON_VERSION}/site-packages/gimic)

include(cython)

configure_file(
    ${PROJECT_SOURCE_DIR}/config.h.in
    ${PROJECT_BINARY_DIR}/config.h
    )

add_subdirectory(tools)

enable_testing()
add_subdirectory(test)
