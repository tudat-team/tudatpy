
# This add_boost.cmake heavily relies on code taken from the MaidSafe project, in specific:
# add_boost.cmake:
#   https://github.com/maidsafe-archive/MaidSafe/blob/next/cmake_modules/utils.cmake
# utils.cmake:
#   https://github.com/maidsafe-archive/MaidSafe/blob/next/cmake_modules/utils.cmake
#
# Original MaidSafe copyright and license below.

#==================================================================================================#
#                                                                                                  #
#  Copyright 2013 MaidSafe.net limited                                                             #
#                                                                                                  #
#  This MaidSafe Software is licensed to you under (1) the MaidSafe.net Commercial License,        #
#  version 1.0 or later, or (2) The General Public License (GPL), version 3, depending on which    #
#  licence you accepted on initial access to the Software (the "Licences").                        #
#                                                                                                  #
#  By contributing code to the MaidSafe Software, or to this project generally, you agree to be    #
#  bound by the terms of the MaidSafe Contributor Agreement, version 1.0, found in the root        #
#  directory of this project at LICENSE, COPYING and CONTRIBUTOR respectively and also available   #
#  at: http://www.maidsafe.net/licenses                                                            #
#                                                                                                  #
#  Unless required by applicable law or agreed to in writing, the MaidSafe Software distributed    #
#  under the GPL Licence is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF   #
#  ANY KIND, either express or implied.                                                            #
#                                                                                                  #
#  See the Licences for the specific language governing permissions and limitations relating to    #
#  use of the MaidSafe Software.                                                                   #
#                                                                                                  #
#==================================================================================================#

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


set(BoostVersion 1.57.0)
set(BoostSHA1 e151557ae47afd1b43dc3fac46f8b04a8fe51c12)

# Create build folder name derived from version
string(REGEX REPLACE "beta\\.([0-9])$" "beta\\1" BoostFolderName ${BoostVersion})
string(REPLACE "." "_" BoostFolderName ${BoostFolderName})
set(BoostFolderName boost_${BoostFolderName})

# If user wants to use a cache copy of Boost, get the path to this location.
if(USE_BOOST_CACHE)
  if(BOOST_CACHE_DIR)
    file(TO_CMAKE_PATH "${BOOST_CACHE_DIR}" BoostCacheDir)
  elseif(WIN32)
    ms_get_temp_dir()
    set(BoostCacheDir "${TempDir}")
  elseif(APPLE)
    set(BoostCacheDir "$ENV{HOME}/Library/Caches")
  else()
    set(BoostCacheDir "$ENV{HOME}/.cache")
  endif()
endif()

# If the cache directory doesn't exist, fall back to use the build root.
if(NOT IS_DIRECTORY "${BoostCacheDir}")
  if(BOOST_CACHE_DIR)
    set(Message "\nThe directory \"${BOOST_CACHE_DIR}\" provided in BOOST_CACHE_DIR doesn't exist.")
    set(Message "${Message}  Falling back to default path at \"${CMAKE_BINARY_DIR}/MaidSafe\"\n")
    message(WARNING "${Message}")
  endif()
  set(BoostCacheDir ${CMAKE_BINARY_DIR})
else()
  if(NOT USE_BOOST_CACHE AND NOT BOOST_CACHE_DIR)
    set(BoostCacheDir "${BoostCacheDir}/MaidSafe")
  endif()
  file(MAKE_DIRECTORY "${BoostCacheDir}")
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

# Build b2 (bjam) if required
unset(b2Path CACHE)
find_program(b2Path NAMES b2 PATHS ${BoostSourceDir} NO_DEFAULT_PATH)
if(NOT b2Path)
  if(MSVC)
    list(APPEND b2Bootstrap "${b2ArgsToolsetPrefix}msvc")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    list(APPEND b2Bootstrap "${b2ArgsToolsetPrefix}clang")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    list(APPEND b2Bootstrap "${b2ArgsToolsetPrefix}gcc")
  endif()
  message(STATUS "Building b2 (bjam)")
  message(STATUS "  ${b2Bootstrap}")
  execute_process(COMMAND ${b2Bootstrap} WORKING_DIRECTORY ${BoostSourceDir}
                  RESULT_VARIABLE Result OUTPUT_VARIABLE Output ERROR_VARIABLE Error)
  if(NOT Result EQUAL "0")
    message(FATAL_ERROR "Failed running bootstrap:\n${Output}\n${Error}\n")
  endif()
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BoostSourceDir}/Build)

# Apply patched files
# if(NOT BoostVersion STREQUAL "1.57.0")
#  message(FATAL_ERROR "Remove patched files from the source tree and delete corresponding 'configure_file' commands in this 'add_boost' CMake file.")
# endif()
# configure_file(patches/boost_1_57/boost/thread/pthread/thread_data.hpp ${BoostSourceDir}/boost/thread/pthread/thread_data.hpp COPYONLY)

# Expose BoostSourceDir to parent scope
set(BoostSourceDir ${BoostSourceDir})

# Set up general b2 (bjam) command line arguments
list(APPEND b2Args link=static threading=multi runtime-link=shared --build-dir=Build stage -d+2 --hash --ignore-site-config)

if(CMAKE_BUILD_TYPE STREQUAL "ReleaseNoInline")
  list(APPEND b2Args "cxxflags=${RELEASENOINLINE_FLAGS}")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "DebugLibStdcxx")
  list(APPEND b2Args define=_GLIBCXX_DEBUG)
endif()

# Set up platform-specific b2 (bjam) command line arguments
if(MSVC)
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
else()
  # Need to configure the toolset based on CMAKE_CXX_COMPILER
  list(APPEND b2Args variant=release cxxflags=-fPIC cxxflags=-std=c++11 --layout=tagged)
  # Apply MinGW32 fix for GCC 4.8.1
  # http://stackoverflow.com/questions/29450016/o1-2-3-with-std-c1y-11-98-if-cmath-is-included-im-getting-error-hypo
  # http://stackoverflow.com/questions/21826649/boost-test-on-windows-with-mingw-compiler-error-putenv-not-declared
  if(MINGW AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
    list(APPEND b2Args cxxflags=-D__NO_INLINE__ cxxflags=-Dputenv=)
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    file(WRITE "${BoostSourceDir}/tools/build/src/user-config.jam" "using clang : : ${CMAKE_CXX_COMPILER} ;\n")
    list(APPEND b2Args toolset=clang architecture=combined address-model=32_64)
    if(HAVE_LIBC++)
      list(APPEND b2Args cxxflags=-stdlib=libc++ linkflags=-stdlib=libc++)
    endif()
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    file(WRITE "${BoostSourceDir}/tools/build/src/user-config.jam" "using gcc : : ${CMAKE_CXX_COMPILER} ;\n")
    list(APPEND b2Args toolset=gcc -sNO_BZIP2=1)
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
    list(APPEND b2Args toolset=clang)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    list(APPEND b2Args toolset=gcc)
  endif()
endif()

# Build only the necessary components
foreach(Component ${BoostComponents})
  if("${Component}" STREQUAL "unit_test_framework")
    set(Component "test")
  endif()
  list(APPEND b2Args "--with-${Component}")
endforeach()

# Start boost build
message(STATUS "Build boost (note that this may take a while, please sit back)")
message(STATUS "  ${b2Args}")
execute_process(COMMAND ${b2Args} WORKING_DIRECTORY ${BoostSourceDir}
                  RESULT_VARIABLE Result OUTPUT_VARIABLE Output ERROR_VARIABLE Error)

                  
# Create the directory and copy source and libraries
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/boost")
file(COPY "${BoostSourceDir}/stage" DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/boost/")
file(COPY "${BoostSourceDir}/boost" DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/boost/")

SET(Boost_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/boost/")

# Clean
# file(REMOVE "${CMAKE_CURRENT_SOURCE_DIR}/${BoostFolderName}.tar.bz2")
