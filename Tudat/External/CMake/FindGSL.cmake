#    Copyright (c) 2010-2016, Delft University of Technology
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
 #      160322    R.Hoogendoorn     File created based on FindSpice.cmake (originally FindEigen3.cmake)
 #
 #
 #    References
 #      FindEigen3.cmake.
 #
 #    Notes
 #      This script tries to find GSL library.
 #


message(STATUS ${GSL_BASE_PATH})

# If the path has not been set previously or manually, try to autodetect the path
if(NOT GSL_BASE_PATH)
    find_path(GSL_BASE_PATH NAMES gsl_version.h
        PATHS
            ${PROJECT_SOURCE_DIR}/External
            ${PROJECT_SOURCE_DIR}/..
            ${PROJECT_SOURCE_DIR}/../..
            ${PROJECT_SOURCE_DIR}/../../..
            ${PROJECT_SOURCE_DIR}/../../../..
            ${PROJECT_SOURCE_DIR}/../../../../..
            ${PROJECT_SOURCE_DIR}/../../../../../..
            ${PROJECT_SOURCE_DIR}/../../../../../../..
        PATH_SUFFIXES "gsl"
)
endif(NOT GSL_BASE_PATH)

# If the path is still not set, then autodetection went wrong
if(NOT GSL_BASE_PATH)

  # Throw a warning and disable GSL
  message(WARNING "WARNING: GSL not found! USE_GSL flag has been disabled.")
  SET(USE_GSL false)

else(NOT GSL_BASE_PATH)

  # Good, path has been set/found, now set important variables and find libraries.
  set(GSL_INCLUDE_DIR "${GSL_BASE_PATH}/include")
  set(GSL_LIBRARIES_DIR "${GSL_BASE_PATH}/lib")
  set(GSL_LIBRARY_DIR "${GSL_BASE_PATH}/lib")
  set(GSL_CBLAS_LIBRARY_DIR "${GSL_BASE_PATH}/lib")

  find_library(GSL_LIBRARY
	NAMES libgsl.a
	PATHS ${GSL_LIBRARIES_DIR})

  find_library(GSL_CBLAS_LIBRARY
	NAMES libgslcblas.a
	PATHS ${GSL_LIBRARIES_DIR})

  set( GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} )

  # Force GSL libraries, to be used when GSL and other libraries are simultaneously compiled.
  if(NOT GSL_LIBRARY)
    set(GSL_LIBRARY "gsl")
  endif( )

  if(NOT GSL_CBLAS_LIBRARY)
    set(GSL_CBLAS_LIBRARY "gslcblas")
  endif( )

  # Let user know which GSL library was found.
  message(STATUS "GSL_LIBRARY: ${GSL_LIBRARY}")
  message(STATUS "GSL_CBLAS_LIBRARY: ${GSL_CBLAS_LIBRARY}")

  link_directories(${GSL_LIBRARIES_DIR})

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GSL DEFAULT_MSG GSL_INCLUDE_DIR)

  mark_as_advanced(GSL_INCLUDE_DIR)

endif(NOT GSL_BASE_PATH)
