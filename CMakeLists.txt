# \file CMakeLists.txt
 #    This is the main CMake file for the Tudat project.
 #
 #    Path              : /
 #    Version           : 4
 #    Check status      : Checked
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Author            : K. Kumar
 #    Affiliation       : Delft University of Technology
 #    E-mail address    : K.Kumar@tudelft.nl
 #
 #    Checker           : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Date created      : 20 August, 2011
 #    Last modified     : 27 October, 2011
 #
 #    References
 #
 #    Notes
 #
 #    Copyright (c) 2010-2011 Delft University of Technology.
 #
 #    This software is protected by national and international copyright.
 #    Any unauthorized use, reproduction or modification is unlawful and
 #    will be prosecuted. Commercial and non-private application of the
 #    software in any form is strictly prohibited unless otherwise granted
 #    by the authors.
 #
 #    The code is provided without any warranty; without even the implied
 #    warranty of merchantibility or fitness for a particular purpose.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      111025    K. Kumar          Adapted file to work with Revision 194.
 #      111027    K. Kumar          Updated to auto-search for Eigen3 and Boost; manual setting
 #                                  examples included as comments.
 #                                  Added suppression of warnings from Boost.
 #/

# Specify minimum CMake version required to compile Tudat.
cmake_minimum_required( VERSION 2.6 )

# Specific Tudat project name.
project( Tudat )

# Set module path to local CMake scripts.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/External/Cmake)

# Set full-path to Tudat trunk directory.
# This should be used in case the FILE token does not given the full-path in the getRootPath()
# function in basicFunctions.cpp.
#add_definitions(-DTUDAT_CUSTOM_ROOT_PATH="C:\My\Custom\Path\To\Tudat\trunk")

# Set root-directory for code to current source directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Deine the basic input (trunk) and output (lib & include) directories.
set(SRCROOT "${CODEROOT}/trunk")

# Set CMake build-type to release.
set(CMAKE_BUILD_TYPE Release)

# Set testing options based on platform.
enable_testing()

if (WIN32)
  if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
    set(LIBROOT "${CODEROOT}/lib_msvc")
    set(INCLUDEROOT "${CODEROOT}/include")
    set(BINROOT "${CODEROOT}/bin_msvc")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(LIBROOT "${CODEROOT}/lib_w32")
    set(INCLUDEROOT "${CODEROOT}/include")
    set(BINROOT "${CODEROOT}/bin_w32")
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  set(LIBROOT "${CODEROOT}/lib")
  set(INCLUDEROOT "${CODEROOT}/include")
  set(BINROOT "${CODEROOT}/bin")
endif()

# Set the global macros for setting up targets.
macro(setup_tudat_executable_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
endmacro(setup_tudat_executable_target)

macro(setup_tudat_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
  install(TARGETS ${target_name} LIBRARY DESTINATION "${BINROOT}" ARCHIVE DESTINATION "${LIBROOT}")
endmacro(setup_tudat_library_target)

macro(setup_custom_test_program target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_custom_test_program)

macro(setup_headers HEADER_FILES HEADER_PATH)
    foreach(CURRENT_HEADER_FILE ${HEADER_FILES})
      install(FILES "${CURRENT_HEADER_FILE}" DESTINATION "${INCLUDEROOT}${HEADER_PATH}")
    endforeach(CURRENT_HEADER_FILE)
  endmacro(setup_headers)

# Set the main sub-directories
set(APPLICATIONSDIR "/Applications")
set(ASTRODYNAMICSDIR "/Astrodynamics")
set(BASICSDIR "/Basics")
set(EXTERNALDIR "/External")
set(INPUTDIR "/Input")
set(MATHEMATICSDIR "/Mathematics")
set(OUTPUTDIR "/Output")

# Include the top-level directories.
include_directories(AFTER
  "${SRCROOT}"
)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)
# If Eigen3 is not found, use the following to manually set the path to the Eigen3 directory.
# A default location is given as an example.
#set(EIGEN3_INCLUDE_DIR "${SRCROOT}/External/eigen3")

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif()

# Configure Boost libraries.
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)

# Find Boost libraries on local system.
find_package(Boost 1.45.0 COMPONENTS thread date_time REQUIRED)
# If Boost is not found, use the following to manually set the path to the Boost include and
# library directories.
#SET(BOOST_INCLUDEDIR "C:\My\Custom\Path\To\Boost\Include\Directory")
#SET(BOOST_LIBRARYDIR "C:\My\Custom\Path\To\Boost\Library\Directory")

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif()

# Set sub-directories.
set(SUBDIRS ${APPLICATIONSDIR} ${ASTRODYNAMICSDIR} ${BASICSDIR} ${INPUTDIR} ${MATHEMATICSDIR} ${OUTPUTDIR})

# Add sub-directories to CMake process.
foreach(CURRENT_SUBDIR ${SUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_SUBDIR}")
endforeach(CURRENT_SUBDIR)

# End of file.
