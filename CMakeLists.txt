#    Copyright (c) 2010-2018, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Specify minimum CMake version required.
cmake_minimum_required(VERSION 2.6)

# Specify project name.
project(TudatPy)

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(tudatLinkLibraries)


# Find Boost libraries on local system.
find_package(Boost 1.55.0
        COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)


FIND_PACKAGE(PythonInterp)

if (PYTHONINTERP_FOUND)
    if (UNIX AND NOT APPLE)
        if (PYTHON_VERSION_MAJOR EQUAL 3)
            FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_SUFFIX})
            FIND_PACKAGE(PythonInterp 3)
            FIND_PACKAGE(PythonLibs 3 REQUIRED)
        else ()
            FIND_PACKAGE(Boost COMPONENTS python)
            FIND_PACKAGE(PythonInterp)
            FIND_PACKAGE(PythonLibs REQUIRED)
        endif ()
    else ()
        if (PYTHON_VERSION_MAJOR EQUAL 3)
            FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
            FIND_PACKAGE(PythonInterp 3)
            FIND_PACKAGE(PythonLibs 3 REQUIRED)
        else ()
            FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
            FIND_PACKAGE(PythonInterp)
            FIND_PACKAGE(PythonLibs REQUIRED)
        endif ()
    endif ()
else ()
    message("Python not found")
endif ()

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")

ENABLE_TESTING()
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} ${PYTHON_LIBRARIES}) # Deprecated but so convenient!

# Add source files.
#file(GLOB_RECURSE SIMULATION_ENVIRONMENT_SETUP_HEADERS ${ENVIRONMENTSETUPDIR} ABSOLUTE ${CODEROOT} *.h)
#file(GLOB_RECURSE SIMULATION_ENVIRONMENT_SETUP_SOURCES ${ENVIRONMENTSETUPDIR} ABSOLUTE ${CODEROOT} *.cpp)
#

# Add static libraries.
#add_library(tudatPy tudat_environment_setup)
#setup_tudat_library_target(tudat_environment_setup "${SRCROOT}${ENVIRONMENTSETUPDIR}")
# Define the directory with the source code.
#set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}")

PYTHON_ADD_MODULE(simulation_setup simulation_setup.cpp)
#FILE(COPY simulation_setup.py DESTINATION .)
#ADD_TEST(NAME SimulationSetup COMMAND ${PYTHON_EXECUTABLE} classes.py)
