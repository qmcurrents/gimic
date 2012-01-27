# - Find the CFOUR includes and library
#
# This module defines
#  CFOUR_LIBRARY_DIR, where to find CFOUR libraries
#  CFOUR_LIBRARIES is set to liblibr.a if found. 
#  CFOUR_FOUND, If false, do not try to use CFOUR.
# also defined, but not for general use are
# None of the above will be defined unles CFOUR can be found.
# 

#=============================================================================
# Copyright 2012 Jonas Juselius <jonas.juselius@uit.no>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

if (CFOUR_LIBRARY_DIR)
	set(CFOUR_FIND_QUIETLY TRUE)
endif ()

find_path(CFOUR_LIBRARY_DIR liblibr.a
  PATHS ${CFOUR_ROOT} $ENV{CFOUR_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
)
find_library(CFOUR_LIBRARIES libr 
    HINTS ${CFOUR_LIBRARY_DIR}
    NO_DEFAULT_PATH
    )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CFOUR DEFAULT_MSG CFOUR_LIBRARIES)

mark_as_advanced(CFOUR_LIBRARY_DIR CFOUR_LIBRARIES)
