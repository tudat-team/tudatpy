cmake_minimum_required(VERSION 2.8.3)

project(tudatpy)

# Add FindNumpy.cmake in boost_numpy_eigen.
list(APPEND CMAKE_MODULE_PATH
        "${CMAKE_CURRENT_SOURCE_DIR}/external/boost_numpy_eigen/cmake")


set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}-Wl,--no-undefined" )

# TudatPy Versioning.
set(tudatpy_VERSION_MAJOR 0)
set(tudatpy_VERSION_MINOR 0)
set(tudatpy_VERSION_PATCH 1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set the main sub-directories.
set(CONSTANTSDIR "/constants")
set(SIMULATIONDIR "/simulation")
set(SPICEDIR "/spice")
set(ENVIRONMENTDIR "/simulation/environment")

# Set library directories.
set(BOOST_EIGEN_NUMPY_SRCDIR "${CMAKE_CURRENT_SOURCE_DIR}/external/boost_numpy_eigen/src")

# `PyLong_Type' can not be used when making a shared object; recompile with -fPIC
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Load UserSettings.txt
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    set(BUILD_STYLE "standalone")
    include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
else ()
    set(BUILD_STYLE "part of ${CMAKE_PROJECT_NAME}")
    include("${CMAKE_CURRENT_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
    include("${CMAKE_SOURCE_DIR}/UserSettings.txt" OPTIONAL)
    STRING(REGEX REPLACE ${CMAKE_SOURCE_DIR} "" RELATIVE_PROJECT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
    set(RELATIVE_PROJECT_PATH "${RELATIVE_PROJECT_PATH}" CACHE STRING "Relative path wrt to project for function")
    # message(STATUS "Relative path (wrt to project): ${RELATIVE_PROJECT_PATH}")
endif ()

# Set CMake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_BUILD_TYPE Release)
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE Debug)
endif ()

message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

# Add local module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules")
message(STATUS "CMake Module path(s): ${CMAKE_MODULE_PATH}")

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(compiler)

# Define the directory with the source code.
set(SRCROOT "${CMAKE_CURRENT_SOURCE_DIR}/tudatpy")

# Define the code root directory.
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# Set testing options based on platform.
enable_testing()

# Set lib and bin directories where static libraries and unit tests are built.
if (NOT LIBROOT)
    set(LIBROOT "${CODEROOT}/lib")
endif ()
if (NOT BINROOT)
    set(BINROOT "${CODEROOT}/bin")
endif ()

# Set the global macros for setting up targets.
macro(setup_tudatpy_library_target target_name CUSTOM_OUTPUT_PATH)
    set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
    set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
endmacro(setup_tudatpy_library_target)

## Set the global macros for setting up targets.
#macro(setup_executable_target target_name CUSTOM_OUTPUT_PATH)
#    set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/applications")
#    install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}/applications")
#endmacro(setup_executable_target)
#
#macro(setup_library_target target_name CUSTOM_OUTPUT_PATH)
#    set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
#    set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
#endmacro(setup_library_target)
#
#macro(setup_unit_test_target target_name CUSTOM_OUTPUT_PATH)
#    set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
#    get_property(CUSTOM_TEST_PROGRAM_NAME TARGET ${target_name} PROPERTY OUTPUT_NAME)
#    add_test("${target_name}" "${BINROOT}/unit_tests/${target_name}")
#endmacro(setup_unit_test_target)

# Include the top-level directories.
include_directories(AFTER
        "${CODEROOT}"
        )

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
if (NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif ()

# Configure Boost libraries.
if (NOT Boost_USE_STATIC_LIBS)
    set(Boost_USE_STATIC_LIBS ON)
endif ()
if (NOT Boost_USE_MULTITHREADED)
    set(Boost_USE_MULTITHREADED ON)
endif ()
if (NOT Boost_USE_STATIC_RUNTIME)
    set(Boost_USE_STATIC_RUNTIME ON)
endif ()

# Find Boost libraries on local system.
find_package(Boost 1.55.0
        COMPONENTS thread date_time system unit_test_framework filesystem regex python3 REQUIRED)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if (NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif ()

# Find Tudat library on local system.
find_package(Tudat 2.0 REQUIRED)

# Include Tudat directories.
# Set CMake flag to suppress Tudat warnings (platform-dependent solution).
if (NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${TUDAT_INCLUDE_DIR}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${TUDAT_INCLUDE_DIR}\"")
endif ()

# Find CSPICE library on local system.
find_package(Spice)

# Include CSpice directories.
if (NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SPICE_INCLUDE_DIR}\"")
endif ()

if (NOT USE_NRLMSISE00)
    message(STATUS "NRLMSISE-00 disabled!")
    add_definitions(-DUSE_NRLMSISE00=0)
else ()
    message(STATUS "NRLMSISE-00 enabled!")
    add_definitions(-DUSE_NRLMSISE00=1)
    # Find USE_NRLMSISE00 library on local system.
    find_package(NRLMSISE00)

    # Include NRLMSISE00 directories.
    if (NOT APPLE OR APPLE_INCLUDE_FORCE)
        include_directories(SYSTEM AFTER "${NRLMSISE00_INCLUDE_DIR}")
    else ()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${NRLMSISE00_INCLUDE_DIR}\"")
    endif ()
endif ()

# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(tudatLinkLibraries)

# TODO: Allow for choice between python 3.6 and python 2.7 for compilation.
FIND_PACKAGE(PythonInterp 3 REQUIRED)

if (PYTHONINTERP_FOUND)
    if (UNIX AND NOT APPLE)
        if (PYTHON_VERSION_MAJOR EQUAL 3)
            # NOTE: As far as is known, building boost.python3 provides the lib files
            # in the form of python3 and not python3.6 as is suggested by original
            # code.
            FIND_PACKAGE(Boost 1.61 COMPONENTS filesystem system regex date_time unit_test_framework thread python${PYTHON_VERSION_MAJOR})
            FIND_PACKAGE(PythonInterp 3)
            FIND_PACKAGE(PythonLibs 3 REQUIRED)
        else ()
            FIND_PACKAGE(Boost 1.61 COMPONENTS filesystem system regex date_time unit_test_framework thread python)
            FIND_PACKAGE(PythonInterp)
            FIND_PACKAGE(PythonLibs REQUIRED)
        endif ()
    else ()
        if (PYTHON_VERSION_MAJOR EQUAL 3)
            FIND_PACKAGE(Boost 1.61 COMPONENTS filesystem system regex date_time unit_test_framework thread python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
            FIND_PACKAGE(PythonInterp 3)
            FIND_PACKAGE(PythonLibs 3 REQUIRED)
        else ()
            FIND_PACKAGE(Boost 1.61 COMPONENTS filesystem system regex date_time unit_test_framework thread python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
            FIND_PACKAGE(PythonInterp)
            FIND_PACKAGE(PythonLibs REQUIRED)
        endif ()
    endif ()
else ()
    message("Python not found")
endif ()

# Add boost_numpy_eigen submodule.
set(POSITION_INDEPENDENT_CODE True)

# Find NumPy
find_package(NumPy ${BOOST_PYTHON_NUMPY_EIGEN_PYTHON_VERSION} REQUIRED)

enable_testing()
include_directories(
        ${Boost_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${PYTHON_INCLUDE_DIRS}
        ${BOOST_NUMPY_EIGEN_SOURCE}
        ${PYTHON_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        ${NUMPY_INCLUDE_DIRS}
        )


add_library(boost_numpy_eigen SHARED
        include/eigen_numpy.cc
        include/eigen_numpy.h)
setup_tudatpy_library_target(boost_numpy_eigen "${BOOST_EIGEN_NUMPY_SRCDIR}")

target_link_libraries(boost_numpy_eigen ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

target_compile_options(boost_numpy_eigen PUBLIC PRIVATE
        -Wignored-attributes
        -Wdeprecated-declarations
        -Wold-style-cast)

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions
#        -Wignored-attributes
#        -Wdeprecated-declarations
#        -Wold-style-cast")

message("${boost_numpy_eigen}")
link_libraries(
        ${Boost_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${PYTHON_LIBRARIES}
        ${tudat_basic_mathematics}
        ${TUDAT_PROPAGATION_LIBRARIES}
        ${tudat_environment_setup}
        ${LIBROOT}/libboost_numpy_eigen.so) # Deprecated but so convenient!

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")
#message(STATUS "boost_numpy_eigen = ${boost_numpy_eigen}")
message(STATUS "NUMPY_INCLUDE_DIRS = ${NUMPY_INCLUDE_DIRS}")
message(STATUS "TUDAT_PROPAGATION_LIBRARIES = ${TUDAT_PROPAGATION_LIBRARIES}")
message(STATUS "tudat_environment_setup = ${tudat_environment_setup}")

# TODO: Setup the CMakelist.txt such that tudatPy is installed to the detected dist-packages of python library.

## Simulation setup tests.
add_test(NAME to_python_conversion COMMAND ${PYTHON_EXECUTABLE} -m unittest -v test/test_to_python_conversion.py)
#add_test(NAME from_python_conversion COMMAND ${PYTHON_EXECUTABLE} -m unittest -v test/test_from_python_conversion.py)
add_test(NAME constants COMMAND ${PYTHON_EXECUTABLE} -m unittest -v test/test_constants.py)
add_test(NAME environment_setup COMMAND ${PYTHON_EXECUTABLE} -m unittest -v test/test_environment_setup.py)

add_subdirectory(test)
add_subdirectory(tudatpy/${CONSTANTSDIR})
add_subdirectory(tudatpy/${SIMULATIONDIR})
add_subdirectory(tudatpy/${SPICEDIR})

