function(gimic_test _name)
  add_test(
    NAME
      ${_name}
    COMMAND
      ./test --binary-dir=${PROJECT_BINARY_DIR}/bin
             --work-dir=${PROJECT_BINARY_DIR}/test/${_name}
             --verbose
    WORKING_DIRECTORY
      ${PROJECT_SOURCE_DIR}/test/${_name}
    )
endfunction()


gimic_test(benzene/integration-gauss)
gimic_test(benzene/integration-lobatto)
gimic_test(benzene/vectors)
gimic_test(benzene/2d)
gimic_test(benzene/3d)
gimic_test(benzene/keyword-magnet)
gimic_test(benzene/3d-keyword-magnet)
gimic_test(benzene/2d-keyword-magnet)
gimic_test(benzene/keyword-spacing)
gimic_test(benzene/keyword-rotation)
gimic_test(benzene/keyword-rotation_origin)
gimic_test(benzene/keyword-radius)
gimic_test(benzene/int-grid-bond-even)
gimic_test(benzene/int-cdens)
gimic_test(benzene/diamag-off)
gimic_test(benzene/paramag-off)
gimic_test(benzene/giao-test)
gimic_test(c4h4/integration)
gimic_test(c4h4/read-grid)
gimic_test(open-shell/3d)
gimic_test(open-shell/integration)
gimic_test(benzene/skip-jmod-integration)
gimic_test(benzene/magnetizability)


option(ENABLE_VALGRIND_TESTS "Enable Valgrind tests" OFF)
if(ENABLE_VALGRIND_TESTS)
  add_test(
    NAME
      benzene/valgrind-integral
    COMMAND
      ./test ${PROJECT_BINARY_DIR}/bin/gimic
    WORKING_DIRECTORY
      ${PROJECT_SOURCE_DIR}/test/benzene/valgrind-integral
    )
  add_test(
    NAME
      benzene/valgrind-2d
    COMMAND
      ./test ${PROJECT_BINARY_DIR}/bin/gimic
    WORKING_DIRECTORY
      ${PROJECT_SOURCE_DIR}/test/benzene/valgrind-2d
    )
endif()
