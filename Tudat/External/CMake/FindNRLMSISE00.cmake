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
 #      FindEigen3.cmake (2-clause BSD license)

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
