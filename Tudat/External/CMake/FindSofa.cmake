# - Try to find Sofa library
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

if(NOT SOFA_BASE_PATH)
find_path(SOFA_BASE_PATH NAMES sofa.h
  PATHS
      ${TUDAT_BASE_PATH}/External
      ${TUDAT_BASE_PATH}/../../..
      ${TUDAT_BASE_PATH}/../..
      ${TUDAT_BASE_PATH}/..
  PATH_SUFFIXES sofa/src
)
endif()

if(NOT SOFA_BASE_PATH)
  message(STATUS "WARNING: SOFA not found! Make sure SOFA is installed")
else()
    MESSAGE( STATUS "Sofa found in: " ${SOFA_BASE_PATH} )

endif( )


set(SOFA_BASE_PATH ${SOFA_BASE_PATH}/..)
set(SOFA_INCLUDE_DIR ${SOFA_BASE_PATH}/..)
set(SOFA_LIBRARIES_DIR ${SOFA_BASE_PATH}/lib)
set(SOFA_LIBRARIES "sofa")
link_directories(${SOFA_LIBRARIES_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SOFA DEFAULT_MSG SOFA_INCLUDE_DIR)

mark_as_advanced(SOFA_INCLUDE_DIR)
