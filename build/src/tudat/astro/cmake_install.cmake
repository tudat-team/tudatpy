# Install script for directory: /home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/src/tudat/astro

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/dominic/miniconda3/envs/tudat-bundle-clean")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/system_models/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/aerodynamics/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/basic_astro/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/earth_orientation/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/electromagnetism/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/ephemerides/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/gravitation/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/ground_stations/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/low_thrust/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/mission_segments/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/observation_models/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/orbit_determination/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/propagators/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/propulsion/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/reference_frames/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/dominic/Tudat/tudat-bundle/tudat-bundle/tudatpy/build/src/tudat/astro/relativity/cmake_install.cmake")
endif()

