include_directories (${PROJECT_SOURCE_DIR}/src)

add_custom_target(tests ALL 
    DEPENDS test.bin
    )

add_executable(test.bin test.f90)
target_link_libraries(test.bin gimic)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/test.sh @ONLY
    )
execute_process(COMMAND 
    chmod 755 ${CMAKE_CURRENT_BINARY_DIR}/test.sh OUTPUT_QUIET) 

configure_file(${CMAKE_SOURCE_DIR}/src/gimic.in
    ${CMAKE_CURRENT_BINARY_DIR}/gimic @ONLY
    )
execute_process(COMMAND 
    chmod 755 ${CMAKE_CURRENT_BINARY_DIR}/gimic OUTPUT_QUIET) 

add_test(NAME ExampleTests 
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test.sh 
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

