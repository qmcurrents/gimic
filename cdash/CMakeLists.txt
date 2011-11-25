#
# Generate the valgrind suppressions file.
#

file(GLOB suppression_files ${CMAKE_SOURCE_DIR}/cdash/*-suppressions.txt)

file(REMOVE ${CMAKE_BINARY_DIR}/valgrind-suppressions.txt)

foreach(fname ${suppression_files})
	file(READ ${fname} vgs)
	file(APPEND ${CMAKE_BINARY_DIR}/valgrind-suppressions.txt ${vgs})
endforeach()

configure_file(${CMAKE_SOURCE_DIR}/cdash/cdash.py.in
	${CMAKE_SOURCE_DIR}/cdash/cdash.py @ONLY
	)
execute_process(COMMAND 
	chmod 755 ${CMAKE_SOURCE_DIR}/cdash/cdash.py OUTPUT_QUIET) 

configure_file(${CMAKE_SOURCE_DIR}/cdash/cdash.run.in
	${CMAKE_SOURCE_DIR}/cdash/cdash.run @ONLY
	)
execute_process(COMMAND 
	chmod 755 ${CMAKE_SOURCE_DIR}/cdash/cdash.run OUTPUT_QUIET) 

