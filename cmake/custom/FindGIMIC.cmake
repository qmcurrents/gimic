# - Find the GIMIC includes and library
#
# This module defines
#  GIMIC_INCLUDE_DIRS, where to find GimicInterface.h, etc.
#  GIMIC_LIBRARIES, the libraries to link against to use GIMIC.
#  GIMIC_DEFINITIONS - You should add_definitons(${GIMIC_DEFINITIONS}) before
#  compiling
#  GIMIC_FOUND, If false, do not try to use GIMIC.
# also defined, but not for general use are
# None of the above will be defined unles GIMIC can be found.
#

#=============================================================================
# Copyright 2010 Jonas Juselius <jonas.juselius@uit.no>
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

if (GIMIC_INCLUDE_DIRS AND GIMIC_LIBRARIES)
    set(GIMIC_FIND_QUIETLY TRUE)
endif ()

find_path(GIMIC_INCLUDE_DIRS
  NAMES GimicInterface.h
  PATHS ${GIMIC_ROOT} $ENV{GIMIC_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
)
find_path(GIMIC_INCLUDE_DIRS GimicInterface.h)

find_path(GIMIC_LIBRARIES gimic2
  PATHS ${GIMIC_ROOT} $ENV{GIMIC_ROOT}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
)
find_library(GIMIC_LIBRARIES gimic2)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GIMIC DEFAULT_MSG
                                  GIMIC_INCLUDE_DIR GIMIC_LIBRARIES)

mark_as_advanced(GIMIC_INCLUDE_DIR GIMIC_LIBRARIES)
