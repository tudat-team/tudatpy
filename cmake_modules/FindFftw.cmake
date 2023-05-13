# - Try to find Fftw library
#
# This file is based on FindEigen3.cmake.
#
# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Copyright (c) 2012 Bryan Tong Minh <b.tongminh@student.tudelft.nl>
# Copyright (c) 2012 Paul Musegaas <p.musegaas@student.tudelft.nl>
# Copyright (c) 2012 Dominic Dirkx <d.dirkx@tudelft.nl>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.

if(NOT FFTW_BASE_PATH)
find_path(FFTW_BASE_PATH NAMES fftw.pc
  PATHS
      ${TUDAT_BASE_PATH}/External
      ${TUDAT_BASE_PATH}/../../../RequiredLibraries
      ${TUDAT_BASE_PATH}/../../RequiredLibraries
      ${TUDAT_BASE_PATH}/../RequiredLibraries
  PATH_SUFFIXES fftw
)
endif()

if(NOT FFTW_BASE_PATH)
  message(STATUS "WARNING: Fftw not found! Make sure Fftw is installed")
else()
    MESSAGE( STATUS "Fftw found in: " ${FFTW_BASE_PATH})
endif( )

set(FFTW_INCLUDE_DIR ${FFTW_BASE_PATH}/api)
set(FFTW_LIBRARIES_DIR ${FFTW_BASE_PATH}/.libs)
set(FFTW_LIBRARIES "fftw3")
link_directories(${FFTW_LIBRARIES_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_INCLUDE_DIR)

mark_as_advanced(FFTW_INCLUDE_DIR)
