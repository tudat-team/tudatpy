 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #
 #    References
 #      FindEigen3.cmake.
 #
 #    Notes
 #      This script tries to find the Tudat Core library. This module supports requiring a minimum
 #      version, e.g. you can do find_package(TudatCore 3.1.2) to require version 3.1.2 or newer of
 #      Tudat Core.
 #
 #      Once done, this will define:
 #
 #          TUDAT_FOUND - system has Tudat Core lib with correct version;
 #          TUDAT_INCLUDE_DIR - the Tudat Core include directory.
 #
 #      Original copyright statements (from FindEigen3.cmake:
 #          Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
 #          Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
 #          Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
 #
 #      FindEigen3.cmake states that redistribution and use is allowed according to the terms of
 #      the 2-clause BSD license.
 #

macro(_tudat_core_check_version)
	MESSAGE( STATUS "Checking for TudatCore in:         " ${TUDAT_BASE_PATH} )
  file(READ "${TUDAT_BASE_PATH}/tudatCoreVersion.h" _tudat_core_header)

  string(REGEX MATCH "define[ \t]+TUDAT_VERSION_MAJOR[ \t]+([0-9]+)" _tudat_core_major_version_match "${_tudat_core_header}")
  set(TUDAT_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+TUDAT_VERSION_MINOR[ \t]+([0-9]+)" _tudat_core_minor_version_match "${_tudat_core_header}")
  set(TUDAT_MINOR_VERSION "${CMAKE_MATCH_1}")

  set(TUDAT_VERSION ${TUDAT_MAJOR_VERSION}.${TUDAT_MINOR_VERSION})
  if(${TUDAT_VERSION} VERSION_LESS ${TudatCore_FIND_VERSION})
    set(TUDAT_VERSION_OK FALSE)
  else(${TUDAT_VERSION} VERSION_LESS ${TudatCore_FIND_VERSION})
    set(TUDAT_VERSION_OK TRUE)
  endif(${TUDAT_VERSION} VERSION_LESS ${TudatCore_FIND_VERSION})

  if(NOT TUDAT_VERSION_OK)

    message(WARNING "Tudat core version ${TUDAT_VERSION} found in ${TUDAT_INCLUDE_DIR}, "
                    "but at least version ${TudatCore_FIND_VERSION} is required")
  endif(NOT TUDAT_VERSION_OK)

  set(TUDAT_LIBRARIES "tudat_core")
  set(TUDAT_INCLUDE_DIR ${TUDAT_BASE_PATH}/..)
  set(TUDAT_LIBRARIES_DIR ${TUDAT_BASE_PATH}/../lib)
  link_directories(${TUDAT_LIBRARIES_DIR})
endmacro(_tudat_core_check_version)

if (TUDAT_BASE_PATH)

  # in cache already
  _tudat_core_check_version( )
  set(TUDAT_FOUND ${TUDAT_VERSION_OK})

else (TUDAT_BASE_PATH)

  find_path(TUDAT_BASE_PATH NAMES tudatCoreVersion.h
      PATHS
      ${PROJECT_SOURCE_DIR}/External
      ${PROJECT_SOURCE_DIR}/../../../tudatCore/trunk
      ${PROJECT_SOURCE_DIR}/../../tudatCore/trunk
      ${PROJECT_SOURCE_DIR}/../../tudatCore
      ${CMAKE_INSTALL_PREFIX}/include
      PATH_SUFFIXES TudatCore
    )

  if(TUDAT_BASE_PATH)
    _tudat_core_check_version( )
  endif(TUDAT_BASE_PATH)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(TudatCore DEFAULT_MSG TUDAT_BASE_PATH TUDAT_VERSION_OK)

  mark_as_advanced(TUDAT_BASE_PATH)

endif(TUDAT_BASE_PATH)
