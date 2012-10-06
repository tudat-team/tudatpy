# - Try to find Spice library
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

find_path(SPICE_BASE_PATH NAMES SpiceUsr.h
  PATHS
      ${PROJECT_SOURCE_DIR}/External
      ${PROJECT_SOURCE_DIR}/../../..
      ${PROJECT_SOURCE_DIR}/../..
      ${PROJECT_SOURCE_DIR}/..
  PATH_SUFFIXES cspice/include
)

if(NOT SPICE_BASE_PATH)
  message(STATUS "WARNING: SPICE not found! Make sure SPICE is installed or set USE_CSPICE to 'false'.")
endif( )

set(SPICE_BASE_PATH ${SPICE_BASE_PATH}/..)
set(SPICE_INCLUDE_DIR ${SPICE_BASE_PATH}/..)
set(SPICE_LIBRARIES_DIR ${SPICE_BASE_PATH}/lib)
set(SPICE_LIBRARIES "cspice")
link_directories(${SPICE_LIBRARIES_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPICE DEFAULT_MSG SPICE_INCLUDE_DIR)

mark_as_advanced(SPICE_INCLUDE_DIR)
