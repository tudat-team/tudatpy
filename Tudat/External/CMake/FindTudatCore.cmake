# - Try to find Tudat Core library
# This module supports requiring a minimum version, e.g. you can do
#   find_package(TudatCore 3.1.2)
# to require version 3.1.2 or newer of Tudat Core.
#
# Once done this will define
#
#  TUDAT_CORE_FOUND - system has tudat core lib with correct version
#  TUDAT_CORE_INCLUDE_DIR - the tudat code include directory
#
# This file is baded on FindEigen3.cmake

# 
# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Copyright (c) 2012 Bryan Tong Minh <b.tongminh@student.tudelft.nl>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.

macro(_tudat_core_check_version)
	MESSAGE( STATUS "Checking for TudatCore in:         " ${TUDAT_CORE_BASE_PATH} )
  file(READ "${TUDAT_CORE_BASE_PATH}/tudatCoreVersion.h" _tudat_core_header)

  string(REGEX MATCH "define[ \t]+TUDAT_CORE_VERSION_MAJOR[ \t]+([0-9]+)" _tudat_core_major_version_match "${_tudat_core_header}")
  set(TUDAT_CORE_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+TUDAT_CORE_VERSION_MINOR[ \t]+([0-9]+)" _tudat_core_minor_version_match "${_tudat_core_header}")
  set(TUDAT_CORE_MINOR_VERSION "${CMAKE_MATCH_1}")

  set(TUDAT_CORE_VERSION ${TUDAT_CORE_MAJOR_VERSION}.${TUDAT_CORE_MINOR_VERSION})
  if(${TUDAT_CORE_VERSION} VERSION_LESS ${TudatCore_FIND_VERSION})
    set(TUDAT_CORE_VERSION_OK FALSE)
  else(${TUDAT_CORE_VERSION} VERSION_LESS ${TudatCore_FIND_VERSION})
    set(TUDAT_CORE_VERSION_OK TRUE)
  endif(${TUDAT_CORE_VERSION} VERSION_LESS ${TudatCore_FIND_VERSION})

  if(NOT TUDAT_CORE_VERSION_OK)

    message(WARNING "Tudat core version ${TUDAT_CORE_VERSION} found in ${TUDAT_CORE_INCLUDE_DIR}, "
                    "but at least version ${TudatCore_FIND_VERSION} is required")
  endif(NOT TUDAT_CORE_VERSION_OK)

  set(TUDAT_CORE_LIBRARIES "tudat_core")
  set(TUDAT_CORE_INCLUDE_DIR ${TUDAT_CORE_BASE_PATH}/..)
  set(TUDAT_CORE_LIBRARIES_DIR ${TUDAT_CORE_BASE_PATH}/../lib)
  link_directories(${TUDAT_CORE_LIBRARIES_DIR})
endmacro(_tudat_core_check_version)

if (TUDAT_CORE_BASE_PATH)

  # in cache already
  _tudat_core_check_version( )
  set(TUDAT_CORE_FOUND ${TUDAT_CORE_VERSION_OK})

else (TUDAT_CORE_BASE_PATH)

  find_path(TUDAT_CORE_BASE_PATH NAMES tudatCoreVersion.h
      PATHS
      ${PROJECT_SOURCE_DIR}/External
      ${PROJECT_SOURCE_DIR}/../../../tudatCore/trunk
      ${PROJECT_SOURCE_DIR}/../../core/trunk
      ${CMAKE_INSTALL_PREFIX}/include
      PATH_SUFFIXES TudatCore
    )

  if(TUDAT_CORE_BASE_PATH)
    _tudat_core_check_version( )
  endif(TUDAT_CORE_BASE_PATH)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(TudatCore DEFAULT_MSG TUDAT_CORE_BASE_PATH TUDAT_CORE_VERSION_OK)

  mark_as_advanced(TUDAT_CORE_BASE_PATH)

endif(TUDAT_CORE_BASE_PATH)
