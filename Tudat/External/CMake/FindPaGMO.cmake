 #    Copyright (c) 2010-2014, Delft University of Technology
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
 #      13xxxx    K. Kumar          Original file
 #      140610    J. Geul           Implemented version finder for main.cpp
 #
 #    References
 #      FindEigen3.cmake.
 #
 #    Notes
 #
 #	This script tries to find the PaGMO library. This module supports requiring a minimum 
 #    	version, e.g. you can do version, e.g. you can do find_package(PaGMO 1.1.2) to require
 #    	version 1.1.2 or newer of PaGMO.
 #
 #	Sets the following variables:
 #          PAGMO_INCLUDE_DIR    - Source directory to include
 #          PAGMO_LIBRARYDIR     - Path to build static library
 #          PAGMO_VERSION_OK     - True of PaGMO found.
 #	     
 #   	    PAGMO_VERSION        - PaGMO version found e.g. 1.1.5
 #   	    PAGMO_VERSION_MAJOR  - PaGMO major version found e.g. 1
 #   	    PAGMO_VERSION_MINOR  - PaGMO minor version found e.g. 1
 #   	    PAGMO_VERSION_PATCH  - PaGMO patch version found e.g. 5
 #
 #      Original file copied from FindEigen3.cmake
 #          Original copyright statements	
 #              Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
 #              Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
 #              Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
 #
 #          FindEigen3.cmake states that redistribution and use is allowed according to the 
 #          terms of the 2-clause BSD license.
 #

macro(_pagmo_check_version)
  message(STATUS "Checking for PaGMO in: " ${PAGMO_BASE_PATH})      
  file(READ "${PAGMO_BASE_PATH}/main.cpp" _pagmo_header)
  STRING(REGEX REPLACE "^.*PaGMO ([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1" PAGMO_VERSION_MAJOR "${_pagmo_header}")
  STRING(REGEX REPLACE "^.*PaGMO [0-9]+\\.([0-9]+)\\.[0-9]+.*" "\\1" PAGMO_VERSION_MINOR "${_pagmo_header}")
  STRING(REGEX REPLACE "^.*PaGMO [0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" PAGMO_VERSION_PATCH "${_pagmo_header}")

  set(PAGMO_VERSION ${PAGMO_VERSION_MAJOR}.${PAGMO_VERSION_MINOR}.${PAGMO_VERSION_PATCH})

  # Only check version if a required is set.
  if(PaGMO_FIND_VERSION)
    if(${PAGMO_VERSION} VERSION_LESS ${PaGMO_FIND_VERSION})
      set(PAGMO_VERSION_OK FALSE)
    else(${PAGMO_VERSION} VERSION_LESS ${PaGMO_FIND_VERSION})
      set(PAGMO_VERSION_OK TRUE)
    endif(${PAGMO_VERSION} VERSION_LESS ${PaGMO_FIND_VERSION})
  else(PaGMO_FIND_VERSION)
      set(PAGMO_VERSION_OK TRUE)
  endif(PaGMO_FIND_VERSION)

  if(NOT PAGMO_VERSION_OK)
    message(STATUS "PaGMO version ${PAGMO_VERSION} found in ${PAGMO_INCLUDE_DIR}, "
                   "but at least version ${PaGMO_FIND_VERSION} is required!")
  endif(NOT PAGMO_VERSION_OK)

  set(PAGMO_LIBRARY "pagmo")
  set(PAGMO_INCLUDE_DIR ${PAGMO_BASE_PATH}/../)
  set(PAGMO_LIBRARYDIR ${PAGMO_BASE_PATH}/build/src/)
  link_directories(${PAGMO_LIBRARYDIR})
endmacro(_pagmo_check_version)

if(PAGMO_BASE_PATH)

  _pagmo_check_version()
  set(PAGMO_FOUND ${PAGMO_VERSION_OK})

else (PAGMO_BASE_PATH)
  find_path(PAGMO_BASE_PATH NAMES main.cpp
      PATHS
      ${PROJECT_SOURCE_DIR}/../pagmo
      ${PROJECT_SOURCE_DIR}/../../pagmo
      ${PROJECT_SOURCE_DIR}/../../../pagmo
      ${PROJECT_SOURCE_DIR}/../../../../pagmo
      ${PROJECT_SOURCE_DIR}/../../../../../pagmo
    )

  if(PAGMO_BASE_PATH)
    _pagmo_check_version()
  endif(PAGMO_BASE_PATH)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(PaGMO DEFAULT_MSG PAGMO_BASE_PATH PAGMO_VERSION_OK)

  mark_as_advanced(PAGMO_BASE_PATH)

endif(PAGMO_BASE_PATH)
