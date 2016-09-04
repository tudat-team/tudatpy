 #    Copyright (c) 2010-2015, Delft University of Technology
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
 #      151120    J. Geul           File created based on FindSpice.cmake (originally FindEigen3.cmake)
 #
 #
 #    References
 #      FindEigen3.cmake.
 #
 #    Notes
 #      This script tries to find NRLMSISE00 library.
 #

# If the path has not been set previously or manually, try to autodetect the path
if(NOT NRLMSISE00_BASE_PATH)
    find_path(NRLMSISE00_BASE_PATH NAMES nrlmsise-00.h
        PATHS
            ${PROJECT_SOURCE_DIR}/External
	    ${PROJECT_SOURCE_DIR}/..
	    ${PROJECT_SOURCE_DIR}/../..
	    ${PROJECT_SOURCE_DIR}/../../..
	    ${PROJECT_SOURCE_DIR}/../../../..
	    ${PROJECT_SOURCE_DIR}/../../../../..
	    ${PROJECT_SOURCE_DIR}/../../../../../..
        PATH_SUFFIXES nrlmsise-00
)
endif(NOT NRLMSISE00_BASE_PATH)

# If the path is still not set, then autodetection went wrong
if(NOT NRLMSISE00_BASE_PATH)

  # Throw a warning and disable NRLMSISE00
  message(WARNING "WARNING: NRLMSISE00 not found! USE_NRLMSISE00 flag has been disabled.")
  SET(USE_NRLMSISE00 false)

else(NOT NRLMSISE00_BASE_PATH)

  # Good, path has been set/found, now set important variables and find libraries.
  set(NRLMSISE00_BASE_PATH ${NRLMSISE00_BASE_PATH})
  set(NRLMSISE00_INCLUDE_DIR ${NRLMSISE00_BASE_PATH})
  set(NRLMSISE00_LIBRARIES_DIR ${NRLMSISE00_BASE_PATH}/lib)

  find_library(NRLMSISE00_LIBRARIES
	NAMES libnrlmsise00.a libnrlmsise00.lib
	PATHS ${NRLMSISE00_LIBRARIES_DIR})

  # Force NRLMSISE00 libraries, to be used when nrlmsise00 and other libraries are simultaneously compiled.
  if(NOT NRLMSISE00_LIBRARIES)
    set(NRLMSISE00_LIBRARIES "nrlmsise00")
  endif( )

  # Let user know which NRLMSISE00 library was found.
  message(STATUS "NRLMSISE00_LIBRARIES: ${NRLMSISE00_LIBRARIES}")

  link_directories(${NRLMSISE00_LIBRARIES_DIR})

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(NRLMSISE00 DEFAULT_MSG NRLMSISE00_INCLUDE_DIR)

  mark_as_advanced(NRLMSISE00_INCLUDE_DIR)

endif(NOT NRLMSISE00_BASE_PATH)
