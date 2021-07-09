#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# Specify minimum CMake version required to compile Tudat.
cmake_minimum_required(VERSION 2.8.10)

# Clean CMakeList Cache
if (EXISTS "${CMAKE_BINARY_DIR}/CMakeCache.txt")
    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeCache.txt")
endif ()

# CMake configuration
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")

# Specific Tudat project name.
project(TudatBundle)

# As both project share the same API docstrings, the solution is to add an external step to generate documented
#  versions of the sources. In this case, tudat-bundle/tudat-docstrings/main.py will need to be executed, if
#  documented versions of the source is desired. This step is OPTIONAL, so if development is not related to
#  docstrings present in the releases, then no need to worry.
# Tudat project.
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.tudat-documented)
    message(STATUS "Using DOCUMENTED tudat")
    add_subdirectory(".tudat-documented")
else ()
    message(STATUS "Using UNDOCUMENTED tudat")
    add_subdirectory(tudat)
endif ()
# TudatPy project.
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.tudatpy-documented)
    message(STATUS "Using DOCUMENTED tudatpy")
    add_subdirectory(".tudatpy-documented")
else ()
    message(STATUS "Using UNDOCUMENTED tudatpy")
    add_subdirectory(tudatpy)
endif ()