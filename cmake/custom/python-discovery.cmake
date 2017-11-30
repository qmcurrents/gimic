# Find Python and NumPy
find_package( PythonLibs 2.7 REQUIRED )
include_directories( ${PYTHON_INCLUDE_DIRS} )
find_package( NumPy REQUIRED )
include_directories( ${NUMPY_INCLUDE_DIRS} )
