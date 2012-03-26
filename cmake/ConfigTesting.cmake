find_program(VALGRIND_EXECUTABLE valgrind)
if (VALGRIND_EXECUTABLE)
    set(MEMORYCHECK_COMMAND ${VALGRIND_EXECUTABLE})
    set(CTEST_MEMORYCHECK_COMMAND ${VALGRIND_EXECUTABLE})
    set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full")
    set(MEMORYCHECK_SUPPRESSIONS_FILE
        ${CMAKE_BINARY_DIR}/valgrind-suppressions.txt)
endif()
mark_as_advanced(VALGRIND_EXECUTABLE)

if (EXISTS ${CMAKE_SOURCE_DIR}/CTestConfig.cmake)
    include(CTest)
    if (EXISTS ${CMAKE_SOURCE_DIR}/cdash)
        set (DASHBOARD_DIR ${CMAKE_SOURCE_DIR}/cdash)
        add_subdirectory(cdash)
    endif()
endif()

if (EXISTS ${CMAKE_SOURCE_DIR}/CTestCustom.cmake.in)
    configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake.in
        ${CMAKE_BINARY_DIR}/CTestCustom.cmake @ONLY)
endif()

enable_testing()

set(BUILDNAME ${BUILDNAME}
    CACHE STRING "Name of build on the dashboard"
    )

