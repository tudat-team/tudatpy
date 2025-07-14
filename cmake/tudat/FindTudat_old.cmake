 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #
 #    References
 #      FindEigen3.cmake (2-clause BSD license)

macro(_tudat_check_version)
  file(READ "${TUDAT_INCLUDE_DIR}/Tudat/tudatVersion.h" _tudat_header)

  string(REGEX MATCH "define[ \t]+TUDAT_VERSION_MAJOR[ \t]+([0-9]+)" _tudat_major_version_match "${_tudat_header}")
  set(TUDAT_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+TUDAT_VERSION_MINOR[ \t]+([0-9]+)" _tudat_minor_version_match "${_tudat_header}")
  set(TUDAT_MINOR_VERSION "${CMAKE_MATCH_1}")

  set(TUDAT_VERSION ${TUDAT_MAJOR_VERSION}.${TUDAT_MINOR_VERSION})
  if(${TUDAT_VERSION} VERSION_LESS ${Tudat_FIND_VERSION})
    set(TUDAT_VERSION_OK FALSE)
  else(${TUDAT_VERSION} VERSION_LESS ${Tudat_FIND_VERSION})
    set(TUDAT_VERSION_OK TRUE)
  endif(${TUDAT_VERSION} VERSION_LESS ${Tudat_FIND_VERSION})

  if(NOT TUDAT_VERSION_OK)

    message(STATUS "tudat version ${TUDAT_VERSION} found in ${TUDAT_INCLUDE_DIR}, "
                   "but at least version ${Tudat_FIND_VERSION} is required")
  endif(NOT TUDAT_VERSION_OK)

  set(TUDAT_LIBRARIES "tudat")
  link_directories(${TUDAT_LIBRARIES_DIR})
endmacro(_tudat_check_version)

if (TUDAT_INCLUDE_DIR)

  # in cache already
  _tudat_check_version()
  set(TUDAT_FOUND ${TUDAT_VERSION_OK})

else (TUDAT_INCLUDE_DIR)

  find_path(TUDAT_BASE_PATH NAMES tudatVersion.h
      PATHS
      ${PROJECT_SOURCE_DIR}/External
      ${PROJECT_SOURCE_DIR}/../../tudat/trunk
      ${PROJECT_SOURCE_DIR}/../../../tudat/trunk
      ${PROJECT_SOURCE_DIR}/../../../../tudat/trunk
      ${PROJECT_SOURCE_DIR}/../../tudat
      ${PROJECT_SOURCE_DIR}/../../../tudat
      ${PROJECT_SOURCE_DIR}/../../../../tudat
      ${CMAKE_INSTALL_PREFIX}/include
      PATH_SUFFIXES Tudat
    )
  set(TUDAT_INCLUDE_DIR ${TUDAT_BASE_PATH}/..)
  set(TUDAT_LIBRARIES_DIR ${TUDAT_BASE_PATH}/../lib)

  if(TUDAT_INCLUDE_DIR)
    _tudat_check_version()
  endif(TUDAT_INCLUDE_DIR)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Tudat DEFAULT_MSG TUDAT_INCLUDE_DIR TUDAT_VERSION_OK)

  mark_as_advanced(TUDAT_INCLUDE_DIR)

  if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" TUDAT_RELATIVE_PROJECT_PATH ${TUDAT_BASE_PATH})
    message(STATUS "Relative path to Tudat found: ${TUDAT_RELATIVE_PROJECT_PATH}")
    add_definitions(-DTUDAT_RELATIVE_PROJECT_PATH="${TUDAT_RELATIVE_PROJECT_PATH}")
  endif()

endif(TUDAT_INCLUDE_DIR)
