#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#
#    Notes
#        This add_boost.cmake heavily relies on code taken from the MaidSafe
#        project, in specific add_boost.cmake and utils.cmake.
#
#        Original MaidSafe copyright distributed under the MaidSafe
#        Contributor Agreement, version 1.0:
#        http://www.maidsafe.net/licenses

#
# Hepler function(s)
#
# Gets the path to the temp directory using the same method as Boost.Filesystem:
# http://www.boost.org/doc/libs/release/libs/filesystem/doc/reference.html#temp_directory_path
function(ms_get_temp_dir)
  if(TempDir)
    return()
  elseif(WIN32)
    file(TO_CMAKE_PATH "$ENV{TEMP}" WindowsTempDir)
    set(Temp "${WindowsTempDir}")
  else()
    foreach(Var TMPDIR TMP TEMP TEMPDIR)
      if(IS_DIRECTORY "$ENV{${Var}}")
        set(Temp $ENV{${Var}})
        break()
      endif()
    endforeach()
    if(NOT TempDir AND IS_DIRECTORY "/tmp")
      set(Temp /tmp)
    endif()
  endif()
  set(TempDir "${Temp}" CACHE INTERNAL "Path to temp directory")
endfunction()

#
# Compiler detection
#
# Allow forcing of building with GNU or Clang (ignore detected)
# NB: do not combine USE_XXX and detection cases in one multi-arg if
# statement, we want overriding before detection.
if( USE_CLANG )
  message(STATUS "BOOST: Using clang.")
  set( BOOST_BUILD_CLANG ON)
elseif( USE_GNU )
  message(STATUS "BOOST: Using gnu.")
  set( BOOST_BUILD_GNU ON)
elseif( "${CMAKE_CXX_COMPILER_ID}" MATCHES "^(Apple)?Clang$" )
  message(STATUS "BOOST: Using clang.")
  set( BOOST_BUILD_CLANG ON)
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  message(STATUS "BOOST: Using gnu.")
  set( BOOST_BUILD_GNU   ON)
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" )
  message(STATUS "BOOST: Using msvc.")
  set( BOOST_BUILD_MSVC ON)
endif()

#
# Boost Source
#
# Create build folder name derived from version
string(REGEX REPLACE "beta\\.([0-9])$" "beta\\1" BoostFolderName ${BoostVersion})
string(REPLACE "." "_" BoostFolderName ${BoostFolderName})
set(BoostFolderName boost_${BoostFolderName})

# Create directory in top-level source dir.
set(BOOST_INCLUDEDIR "${CMAKE_CURRENT_SOURCE_DIR}/boost")
set(BOOST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/boost")
set(BOOST_LIBRARYDIR "${BOOST_INCLUDEDIR}/stage/lib")
set(Boost_NO_SYSTEM_PATHS ON)

set(BoostCacheDir   "${BOOST_INCLUDEDIR}/build")
file(MAKE_DIRECTORY "${BOOST_INCLUDEDIR}")
file(MAKE_DIRECTORY "${BoostCacheDir}")

# Force using static libraries as the build libraries are not installed to the system
# or the libs dir added to the path.
set(Boost_USE_STATIC_LIBS ON)

#
# Check if C POSIX library is present
#
include(CheckIncludeFiles)
# usage: CHECK_INCLUDE_FILES (<header> <RESULT_VARIABLE> )
unset(HAVE_UNISTD_H CACHE)
check_include_files(unistd.h HAVE_UNISTD_H)
if(NOT HAVE_UNISTD_H)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        execute_process(COMMAND xcode-select;-p OUTPUT_VARIABLE Output ERROR_VARIABLE Error)
        message(FATAL_ERROR "C POSIX libraries not found on system. From command line execute:\nxcode-select --install.\n${Output}\n${Error}\n")
    else()
        message(FATAL_ERROR "C POSIX libraries not found on system. Please (re)install development environment.")
    endif()
endif()

#
# Check if local Boost is not already present
#
# Prevent find_package from throwing an error with Boost_NO_SYSTEM_PATHS ON
if(EXISTS "${BOOST_INCLUDEDIR}/boost/version.hpp")
# Set variables to true, if boost is already on system these stay true!
set(BoostComponentsFound ON)
set(BoostComponentsDir ON)
foreach(Component ${BoostComponents})
  find_package(Boost ${BoostVersion} COMPONENTS ${Component} QUIET)

  # Convert component name to upper case
  string(TOUPPER "${Component}" ComponentUpper)

  # Variable variable names! Second-level unwrapping. 
  # message(STATUS "${ComponentUpper}: ${Boost_${ComponentUpper}_FOUND} - ${Boost_${ComponentUpper}_LIBRARY}")
  if(NOT ${Boost_${ComponentUpper}_FOUND})
    set(BoostComponentsFound OFF)
	message(STATUS "Boost ${Component} not found on system. We will build it then. Please ignore boost warnings.")
  endif()

  # Check if the library is located in the local library directory 
  string(FIND "${Boost_${ComponentUpper}_LIBRARY}" "${BOOST_INCLUDEDIR}" BoostComponentsDirYes)
  if(BoostComponentsDirYes LESS 0)
    set(BoostComponentsDir OFF)
  endif()
  
  # Need to unset these too, otherwise other find_package calls willl not update them. Also some are put in current scope and cache
  unset("Boost_${ComponentUpper}_FOUND")
  unset("Boost_${ComponentUpper}_LIBRARY")
  unset("Boost_${ComponentUpper}_FOUND" CACHE)
  unset("Boost_${ComponentUpper}_LIBRARY" CACHE)
  unset("Boost_${ComponentUpper}_LIBRARIES")
  unset("Boost_${ComponentUpper}_LIBRARY_DEBUG" CACHE)
  unset("Boost_${ComponentUpper}_LIBRARY_RELEASE" CACHE)
  
  # Exit the for loop if a single component fails
  if(NOT ${BoostComponentsDir})
    break()
  endif()
  if(NOT ${BoostComponentsFound})
    break()
  endif()

endforeach()

# Unset all variable from find_package(Boost), preventing future usages of this macro becoming lazy.
# As before some variables are also cached, so need double cleaning
unset(Boost_FOUND)
unset(Boost_INCLUDE_DIRS)
unset(Boost_LIBRARY_DIRS)
unset(Boost_LIBRARY_DIRS CACHE)
unset(Boost_LIBRARIES)
unset(Boost_VERSION)
unset(Boost_VERSION CACHE)
unset(Boost_LIB_VERSION)
unset(Boost_LIB_VERSION CACHE)
unset(Boost_MAJOR_VERSION)
unset(Boost_MINOR_VERSION)
unset(Boost_SUBMINOR_VERSION)

endif()

# Check if all components were found and if their location is local and not on the system.
if(${BoostComponentsFound} AND ${BoostComponentsDir})
  message(STATUS "Boost was already build on system")
  return()
endif()


# Set up the full path to the source directory
set(BoostSourceDir "${BoostFolderName}_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}")
if(HAVE_LIBC++)
  set(BoostSourceDir "${BoostSourceDir}_LibCXX")
endif()
if(HAVE_LIBC++ABI)
  set(BoostSourceDir "${BoostSourceDir}_LibCXXABI")
endif()
if(CMAKE_CL_64)
  set(BoostSourceDir "${BoostSourceDir}_Win64")
endif()
string(REPLACE "." "_" BoostSourceDir ${BoostSourceDir})
set(BoostSourceDir "${BoostCacheDir}/${BoostSourceDir}")

# Check the full path to the source directory is not too long for Windows.  File paths must be less
# than MAX_PATH which is 260.  The current longest relative path Boost tries to create is:
# Build\boost\bin.v2\libs\program_options\build\fd41f4c7d882e24faa6837508d6e5384\libboost_program_options-vc120-mt-gd-1_55.lib.rsp
# which along with a leading separator is 129 chars in length.  This gives a maximum path available
# for 'BoostSourceDir' as 130 chars.
if(WIN32)
  get_filename_component(BoostSourceDirName "${BoostSourceDir}" NAME)
  string(LENGTH "/${BoostSourceDirName}" BoostSourceDirNameLengthWithSeparator)
  math(EXPR AvailableLength 130-${BoostSourceDirNameLengthWithSeparator})
  string(LENGTH "${BoostSourceDir}" BoostSourceDirLength)
  if(BoostSourceDirLength GREATER "130")
    set(Msg "\n\nThe path to boost's source is too long to handle all the files which will ")
    set(Msg "${Msg}be created when boost is built.  To avoid this, set the CMake variable ")
    set(Msg "${Msg}USE_BOOST_CACHE to ON and set the variable BOOST_CACHE_DIR to a path ")
    set(Msg "${Msg}which is at most ${AvailableLength} characters long.  For example:\n")
    set(Msg "${Msg}  mkdir C:\\maidsafe_boost\n")
    set(Msg "${Msg}  cmake . -DUSE_BOOST_CACHE=ON -DBOOST_CACHE_DIR=C:\\maidsafe_boost\n\n")
    message(FATAL_ERROR "${Msg}")
  endif()
endif()

# Download boost if required
set(ZipFilePath "${BoostCacheDir}/${BoostFolderName}.tar.bz2")
if(NOT EXISTS ${ZipFilePath})
  message(STATUS "Downloading boost ${BoostVersion} to ${BoostCacheDir}")
endif()
file(DOWNLOAD http://sourceforge.net/projects/boost/files/boost/${BoostVersion}/${BoostFolderName}.tar.bz2/download
     ${ZipFilePath}
     STATUS Status
     SHOW_PROGRESS
     EXPECTED_HASH SHA1=${BoostSHA1}
     )

# Extract boost if required
string(FIND "${Status}" "returning early" Found)
if(Found LESS "0" OR NOT IS_DIRECTORY "${BoostSourceDir}")
  set(BoostExtractFolder "${BoostCacheDir}/boost_unzip")
  file(REMOVE_RECURSE ${BoostExtractFolder})
  file(MAKE_DIRECTORY ${BoostExtractFolder})
  file(COPY ${ZipFilePath} DESTINATION ${BoostExtractFolder})
  message(STATUS "Extracting boost ${BoostVersion} to ${BoostExtractFolder}")
  if(WIN32)
    message(STATUS "On Windows this can take up to several (tens) of minutes.")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz ${BoostFolderName}.tar.bz2
                  WORKING_DIRECTORY ${BoostExtractFolder}
                  RESULT_VARIABLE Result
                  )
  if(NOT Result EQUAL "0")
    message(FATAL_ERROR "Failed extracting boost ${BoostVersion} to ${BoostExtractFolder}")
  endif()
  file(REMOVE ${BoostExtractFolder}/${BoostFolderName}.tar.bz2)

  # Get the path to the extracted folder
  file(GLOB ExtractedDir "${BoostExtractFolder}/*")
  list(LENGTH ExtractedDir n)
  if(NOT n EQUAL "1" OR NOT IS_DIRECTORY ${ExtractedDir})
    message(FATAL_ERROR "Failed extracting boost ${BoostVersion} to ${BoostExtractFolder}")
  endif()
  file(RENAME ${ExtractedDir} ${BoostSourceDir})
  file(REMOVE_RECURSE ${BoostExtractFolder})
endif()

# Set platform depenedent executables
if(WIN32)
  set(b2Bootstrap ".\\bootstrap.bat")
  set(b2Args ".\\b2.exe")
  set(b2ArgsToolsetPrefix "")
else()
  set(b2Bootstrap "./bootstrap.sh")
  set(b2Args "./b2")
  set(b2ArgsToolsetPrefix "--with-toolset=")
endif()

#
# Custom config
#
set(BoostCmakeConfig "${BoostSourceDir}/cmake-config.jam")
file(WRITE "${BoostCmakeConfig}" "# Generated automatically by: add_boost.cmake\n\n")

#
# Python magic
#
# If python library is build, we need to point to the interpreter, the version and directories
if(";${BoostComponents}" MATCHES ";python")

  # Use either automatically detected or user set variables
  find_package(PythonInterp)
  find_package(PythonLibs)
  if(NOT BoostPythonExe AND NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "Please set the BoostPythonExe variable to your Python executable.")
  elseif(NOT BoostPythonExe)
    set(BoostPythonExe "${PYTHON_EXECUTABLE}")
  endif()
  if(NOT BoostPythonVer AND NOT PYTHON_VERSION_STRING)
    message(FATAL_ERROR "Please set the BoostPythonVer variable to your Python version (major.minor).")
  elseif(NOT BoostPythonVer)
    set(BoostPythonVer "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
  endif()
  if(NOT BoostPythonInc AND NOT PYTHON_INCLUDE_DIRS)
    message(FATAL_ERROR "Please set the BoostPythonInc variable to your Python include directory.")
  elseif(NOT BoostPythonInc)
    set(BoostPythonInc "${PYTHON_INCLUDE_DIRS}")
  endif()
  if(NOT BoostPythonLib AND NOT PYTHON_LIBRARIES)
    message(FATAL_ERROR "Please set the BoostPythonLib variable to your Python libs directory.")
  elseif(NOT BoostPythonLib)
    list(LENGTH ${PYTHON_LIBRARIES} PYTHON_LIBRARIES_COUNT)
    if(${PYTHON_LIBRARIES_COUNT} GREATER 0)
      list(GET ${PYTHON_LIBRARIES} 0 BoostPythonLib)
    else()
     set(BoostPythonLib ${PYTHON_LIBRARIES})
    endif()
    get_filename_component(BoostPythonLib ${BoostPythonLib} PATH)
  endif()

  # Add this to user configuration (used by b2)
  file(APPEND "${BoostCmakeConfig}"
    "using python\n    : ${BoostPythonVer}\n    : \"${BoostPythonExe}\"\n    : \"${BoostPythonInc}\"\n    : \"${BoostPythonLib}\" ;\n")
  # Need to add extra bootstrap arguments on Linux and OS X, bootstrap
  # detection of Python is flawed (e.g. wrong library directories
  # etcetera), moreover, this ends up in project-config.jam and takes
  # priority over whatever we set in our custom user-config.jam.
  if(NOT WIN32)
    list(APPEND b2Bootstrap "--with-python=purposefully_failing_python_detection")
  endif()
endif()

#
# Bootstrap
#
# Build b2 (bjam) if required
unset(b2Path CACHE)
find_program(b2Path NAMES b2 PATHS ${BoostSourceDir} NO_DEFAULT_PATH)
if(NOT b2Path)
  if( BOOST_BUILD_MSVC )
    list(APPEND b2Bootstrap "${b2ArgsToolsetPrefix}msvc")
  elseif( BOOST_BUILD_CLANG )
    list(APPEND b2Bootstrap "${b2ArgsToolsetPrefix}clang")
  elseif( BOOST_BUILD_GNU )
    list(APPEND b2Bootstrap "${b2ArgsToolsetPrefix}gcc")
  endif()
  message(STATUS "Building b2 (bjam)")
  message(STATUS "  ${b2Bootstrap}")
  execute_process(COMMAND ${b2Bootstrap} WORKING_DIRECTORY ${BoostSourceDir}
                  RESULT_VARIABLE Result OUTPUT_VARIABLE Output ERROR_VARIABLE Error)
  file(WRITE "${BoostSourceDir}/build_bootstrap.log" "-- ARGS\n${b2Bootstrap}\n\n-- RESULT\n${Result}\n\n-- OUTPUT\n${Output}\n\n-- ERROR\n${Error}")
  if(NOT Result EQUAL "0")
    message(FATAL_ERROR "Failed running bootstrap:\n${Output}\n${Error}\n")
  endif()
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BoostSourceDir}/Build)

#
# Patching
#
# Apply patched files
# if(NOT BoostVersion STREQUAL "1.57.0")
#  message(FATAL_ERROR "Remove patched files from the source tree and delete corresponding 'configure_file' commands in this 'add_boost' cmake file.")
# endif()
# configure_file(patches/boost_1_57/boost/thread/pthread/thread_data.hpp ${BoostSourceDir}/boost/thread/pthread/thread_data.hpp COPYONLY)

# Expose BoostSourceDir to parent scope
set(BoostSourceDir ${BoostSourceDir})

#
# B2 and compiler configuration
#
# Set up general b2 (bjam) command line arguments
list(APPEND b2Args --user-config=cmake-config.jam link=static,shared threading=multi --build-dir=Build stage)

if(CMAKE_BUILD_TYPE STREQUAL "ReleaseNoInline")
  list(APPEND b2cxxflags "${RELEASENOINLINE_FLAGS}")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "DebugLibStdcxx")
  list(APPEND b2Args define=_GLIBCXX_DEBUG)
endif()
list(APPEND b2cxxflags "-fPIC")

# Apply ::hypot not declared patch for MinGW32
# http://stackoverflow.com/questions/10660524/error-building-boost-1-49-0-with-gcc-4-7-0
if(";${BoostComponents}" MATCHES ";python" AND MINGW)
  list(APPEND b2cxxflags "-include cmath")
endif()

# Set up platform-specific b2 (bjam) command line arguments
if( BOOST_BUILD_MSVC )
  if(MSVC11)
    list(APPEND b2Args toolset=msvc-11.0)
  elseif(MSVC12)
    list(APPEND b2Args toolset=msvc-12.0)
  elseif(MSVC14)
    list(APPEND b2Args toolset=msvc-14.0)
  endif()
  list(APPEND b2Args define=_BIND_TO_CURRENT_MFC_VERSION=1 define=_BIND_TO_CURRENT_CRT_VERSION=1 --layout=versioned)
  if(TargetArchitecture STREQUAL "x86_64")
    list(APPEND b2Args address-model=64)
  endif()
  set(b2toolset "msvc")
  set(b2cxxexe  "")
else()
  # Apply MinGW32 fix for GCC 4.8.1
  # http://stackoverflow.com/questions/29450016/o1-2-3-with-std-c1y-11-98-if-cmath-is-included-im-getting-error-hypo
  # http://stackoverflow.com/questions/21826649/boost-test-on-windows-with-mingw-compiler-error-putenv-not-declared
  if(MINGW AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
    list(APPEND b2CXX -D__NO_INLINE__ -Dputenv=)
  endif()
  if( BOOST_BUILD_CLANG )
    list(APPEND b2Args architecture=combined address-model=32_64)
    if(HAVE_LIBC++)
      list(APPEND b2Args linkflags=-stdlib=libc++)
	  list(APPEND b2cxxflags  -stdlib=libc++)
    endif()
    set(b2toolset "clang")
  elseif( BOOST_BUILD_GNU )
    list(APPEND b2Args -sNO_BZIP2=1)
    set(b2toolset "gcc")
  endif()
  # Need to configure the toolset based on CMAKE_CXX_COMPILER
  list(APPEND b2Args variant=release --layout=tagged toolset=${b2toolset})
  list(APPEND b2cxxflags -std=c++11)
  set(b2cxxexe ${CMAKE_CXX_COMPILER})
endif()
string (REPLACE ";" " " b2cxxflagsstr "${b2cxxflags}")
file(APPEND "${BoostCmakeConfig}"
  "using ${b2toolset}\n    :\n    : \"${b2cxxexe}\"\n    : <cxxflags>\"${b2cxxflagsstr}\" ;\n")

# Build only the necessary components
foreach(Component ${BoostComponents})
  if("${Component}" STREQUAL "unit_test_framework")
    set(Component "test")
  elseif("${Component}" STREQUAL "python3")
    # Check if python is also build, otherwise we still need it
    if(";${BoostComponents};" MATCHES ";python;")
      continue()
    else()
      set(Component "python")
    endif()
  endif()
  list(APPEND b2Args "--with-${Component}")
endforeach()

#
# B2 run
#
# Start boost build
message(STATUS "Build boost (note that this may take a while, please sit back)")
message(STATUS "  ${b2Args}")
execute_process(COMMAND ${b2Args} WORKING_DIRECTORY ${BoostSourceDir}
                  RESULT_VARIABLE Result OUTPUT_VARIABLE Output ERROR_VARIABLE Error)
file(WRITE "${BoostSourceDir}/build_b2.log"
  "-- ARGS\n${b2Args}\n\n-- RESULT\n${Result}\n\n-- OUTPUT\n${Output}\n\n-- ERROR\n${Error}")

#
# Copy libraries and source
#
# Set FindBoost.cmake helpers and copy source and libraries
file(COPY "${BoostSourceDir}/stage" DESTINATION "${BOOST_INCLUDEDIR}")
file(COPY "${BoostSourceDir}/boost" DESTINATION "${BOOST_INCLUDEDIR}")

# Clean
# file(REMOVE "${CMAKE_CURRENT_SOURCE_DIR}/${BoostFolderName}.tar.bz2")
