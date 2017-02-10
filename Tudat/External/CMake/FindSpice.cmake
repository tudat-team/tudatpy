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
 #      FindEigen3.cmake licensed under the 2-clause BSD license.

# If the path has not been set previously or manually, try to autodetect the path
if(NOT SPICE_BASE_PATH)
    find_path(SPICE_BASE_PATH NAMES SpiceUsr.h
        PATHS
            ${PROJECT_SOURCE_DIR}/External
	    ${PROJECT_SOURCE_DIR}/..
	    ${PROJECT_SOURCE_DIR}/../..
	    ${PROJECT_SOURCE_DIR}/../../..
	    ${PROJECT_SOURCE_DIR}/../../../..
        PATH_SUFFIXES cspice/include
)
endif(NOT SPICE_BASE_PATH)

# If the path is still not set, then autodetection went wrong
if(NOT SPICE_BASE_PATH)

  # Throw a warning and disable SPICE
  message(WARNING "WARNING: SPICE not found! USE_CSPICE flag has been disabled.")
  SET(USE_CSPICE false)

else(NOT SPICE_BASE_PATH)

  # Good, path has been set/found, now set important variables and find libraries.
  set(SPICE_BASE_PATH ${SPICE_BASE_PATH}/..)
  set(SPICE_INCLUDE_DIR ${SPICE_BASE_PATH}/..)
  set(SPICE_LIBRARIES_DIR ${SPICE_BASE_PATH}/lib)

  find_library(SPICE_LIBRARIES
	NAMES libcspice.a libcspice.lib cspice.a cspice.lib
	PATHS ${SPICE_LIBRARIES_DIR})

  # Force SPICE libraries, to be used when spice and other libraries are simultaneously compiled.
  if(NOT SPICE_LIBRARIES)
    set(SPICE_LIBRARIES "cspice")
  endif( )

  # Let user know which SPICE library was found.
  message(STATUS "SPICE_LIBRARIES: ${SPICE_LIBRARIES}")

  link_directories(${SPICE_LIBRARIES_DIR})

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(SPICE DEFAULT_MSG SPICE_INCLUDE_DIR)

  mark_as_advanced(SPICE_INCLUDE_DIR)

endif(NOT SPICE_BASE_PATH)
