find_program(VALGRIND_EXECUTABLE valgrind)
if (VALGRIND_EXECUTABLE)
    set(MEMORYCHECK_COMMAND ${VALGRIND_EXECUTABLE})
    set(CTEST_MEMORYCHECK_COMMAND ${VALGRIND_EXECUTABLE})
    set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full")
    set(MEMORYCHECK_SUPPRESSIONS_FILE
        ${CMAKE_BINARY_DIR}/valgrind-suppressions.txt)
endif()
mark_as_advanced(VALGRIND_EXECUTABLE)

enable_testing()
