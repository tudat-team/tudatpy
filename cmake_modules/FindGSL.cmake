 #    Copyright (c) 2010-2018, Delft University of Technology
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
