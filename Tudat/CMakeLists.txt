 #    Copyright (c) 2010-2012 Delft University of Technology.
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
 #                                  examples included as comments. Added suppression of warnings
 #                                  from Boost.
 #      111028    K. Kumar          Corrected Boost variables and changed build-type to Debug.
 #      111103    K. Kumar          Updated CMake directory name in module path.
 #      111117    K. Kumar          Added system and filesystem Boost components.
 #      170212    S. Billemont      Added support for UserSettings.txt file.
 #
 #    References
 #

# Specify minimum CMake version required to compile Tudat.
cmake_minimum_required( VERSION 2.6 )

# Specific Tudat project name.
project( Tudat )

# Set module path to local CMake scripts.
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/External/CMake)

# Set root-directory for code to current source directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# Define the basic input (trunk) and output (lib & include) directories.
set(SRCROOT "${CODEROOT}/Tudat")

# Set default build-type.
set(CMAKE_BUILD_TYPE Debug)

# Set testing options based on platform.
enable_testing( )

if (WIN32)
  if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Ox /W3 /FC -D_SCL_SECURE_NO_WARNINGS")
    # Because we are using static boost libraries, with static runtime, we need to force MSVC to
    # also use static runtime: (from http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace)
    foreach(flag_var
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      # Find all dynamic runtime (MD) references and replace with static (MT)
      if(${flag_var} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
      endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
    set(LIBROOT "${CODEROOT}/lib")
    set(INCLUDEROOT "${CODEROOT}/include")
    set(BINROOT "${CODEROOT}/bin")
  else( )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(LIBROOT "${CODEROOT}/lib")
    set(INCLUDEROOT "${CODEROOT}/include")
    set(BINROOT "${CODEROOT}/bin")
  endif( )
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(LIBROOT "${CODEROOT}/lib")
    set(INCLUDEROOT "${CODEROOT}/include")
    set(BINROOT "${CODEROOT}/bin")
endif( )

# Defines the install target to create a distribution of Tudat.
if(NOT TUDAT_DISTRIBUTION_PATH)
    set(TUDAT_DISTRIBUTION_PATH "${CODEROOT}/dist/packages")
endif(NOT TUDAT_DISTRIBUTION_PATH)

# ------------------------------------------------------------------------------------
# Don't define any new variables with SET that the user might change beyond this point.
#
# Load user specific settings. These can overwrite all previous settings.
if (NOT TUDAT_USER_SETTINGS AND NOT DEFINED ENV{TUDAT_USER_SETTINGS})
  # Load from default paths, least important first.
  INCLUDE( "../UserSettings.txt" OPTIONAL)
  INCLUDE( "UserSettings.txt" OPTIONAL)
else()
	MESSAGE(STATUS "Loading preset")
  # Load the user specified file
  if ( TUDAT_USER_SETTINGS )
    INCLUDE( $TUDAT_USER_SETTINGS OPTIONAL)
  else()  
    INCLUDE( $ENV{TUDAT_USER_SETTINGS} OPTIONAL)
  endif()
endif()

# Install the Tudat files.
install(DIRECTORY "${SRCROOT}/"
        DESTINATION "${TUDAT_DISTRIBUTION_PATH}/Tudat"
        FILES_MATCHING
            PATTERN "*.*"
            PATTERN "CMakeLists.txt.user" EXCLUDE
) 

# Set the global macros for setting up targets.
macro(setup_tudat_library_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
  install(TARGETS ${target_name} LIBRARY DESTINATION "${BINROOT}" ARCHIVE DESTINATION "${TUDAT_DISTRIBUTION_PATH}/lib")
endmacro(setup_tudat_library_target)

macro(setup_custom_test_program target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
  get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
  add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_custom_test_program)

macro(setup_headers HEADER_FILES HEADER_PATH)
    foreach(CURRENT_HEADER_FILE ${HEADER_FILES})
      install(FILES "${CURRENT_HEADER_FILE}" DESTINATION "${HEADER_PATH}")
    endforeach(CURRENT_HEADER_FILE)
endmacro(setup_headers)

# Set the main sub-directories.
set(ASTRODYNAMICSDIR "/Astrodynamics")
set(EXTERNALDIR "/External")
set(INPUTOUTPUTDIR "/InputOutput")
set(MATHEMATICSDIR "/Mathematics")

# Include the top-level directories.
include_directories(AFTER
  "${CODEROOT}"
)

# Load the environment file from the core distribution if specified.
if(TUDAT_CORE_ENVIRONMENT_PATH)
    message(STATUS "Loading Tudat Core distribution from ${TUDAT_CORE_ENVIRONMENT_PATH}/TudatCoreEnvironment.cmake")
    include("${TUDAT_CORE_ENVIRONMENT_PATH}/TudatCoreEnvironment.cmake")
endif(TUDAT_CORE_ENVIRONMENT_PATH)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif( )

# Configure Boost libraries.
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)

# Find Boost libraries on local system.
find_package(Boost 1.45.0 COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif( )

# Find Tudat Core library on local system.
find_package(TudatCore 1.1 REQUIRED)

# Include TudatCore directories.
# Set CMake flag to suppress TudatCore warnings (platform-dependent solution).
if(NOT APPLE)
  include_directories(SYSTEM AFTER "${TUDAT_CORE_INCLUDE_DIR}")
else( )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_CORE_INCLUDE_DIR}\"")
endif( )

# Set sub-directories.
set(SUBDIRS ${APPLICATIONSDIR} ${ASTRODYNAMICSDIR} ${INPUTDIR} ${INPUTOUTPUTDIR} ${MATHEMATICSDIR})

# Add sub-directories to CMake process.
foreach(CURRENT_SUBDIR ${SUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_SUBDIR}")
endforeach(CURRENT_SUBDIR)
