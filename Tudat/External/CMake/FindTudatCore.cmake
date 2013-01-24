 #    Copyright (c) 2010-2013, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      12xxxx    B. Tong Minh      File created based on FindEigen3.cmake.
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
 #          TUDAT_CORE_FOUND - system has Tudat Core lib with correct version;
 #          TUDAT_CORE_INCLUDE_DIR - the Tudat Core include directory.
 #
 #      Original copyright statements (from FindEigen3.cmake:
 #          Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
 #          Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
 #          Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
 #
 #      FindEigen3.cmake states that redistribution and use is allowed according to the terms of
 #      the 2-clause BSD license.

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
      ${PROJECT_SOURCE_DIR}/../../tudatCore/trunk
      ${PROJECT_SOURCE_DIR}/../tudatCore/
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
