message("-- Detected CMAKE_SYSTEM_NAME is \"${CMAKE_SYSTEM_NAME}\"")

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
	#if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i686")
	#    add_definitions(-D32BIT)
	#elseif(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
	#else()
	#endif()
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
	#add_definitions(-DSYS_DARWIN)
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
endif()

