foreach(_name
    benzene/integration
    benzene/bond-integration
    benzene/vectors
    benzene/2d
    benzene/3d
    benzene/keyword-magnet
    benzene/3d-keyword-magnet
    benzene/2d-keyword-magnet
    benzene/keyword-spacing
    benzene/keyword-rotation
    benzene/keyword-rotation_origin
    benzene/keyword-radius
    benzene/int-grid-bond-even
    benzene/int-cdens
    benzene/diamag-off
    benzene/paramag-off
    benzene/giao-test   
    c4h4/integration
    c4h4/read-grid
    open-shell/3d
    open-shell/integration
    benzene/valgrind-integral
    benzene/valgrind-2d
    )
    add_test(
        NAME ${_name}
        COMMAND ./test
                --binary-dir=${PROJECT_BINARY_DIR}/bin
                --work-dir=${PROJECT_BINARY_DIR}/test/${_name}
                --verbose
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test/${_name}
        )
endforeach()
