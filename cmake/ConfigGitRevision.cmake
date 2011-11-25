set (GIT_REVISION)
find_package(Git)
if (GIT_FOUND)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} rev-list --abbrev-commit --max-count=1 HEAD
		OUTPUT_VARIABLE GIT_REVISION
		)
	string(STRIP
		${GIT_REVISION}
		GIT_REVISION)
endif()

