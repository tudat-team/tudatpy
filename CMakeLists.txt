cmake_minimum_required(VERSION 2.8.3)

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

# Get version+label from file.
file(STRINGS "version" tudatpy_VERSION_LABEL)

# Get version from version_label argument using regex (x.y.z-rc -> x.y.z).
string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)" tudatpy_VERSION "${tudatpy_VERSION_LABEL}")

# Project declaration.
project(tudatpy VERSION "${tudatpy_VERSION}" LANGUAGES CXX C)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
else()
#    add_definitions(-DVERSION_INFO "${tudatpy_VERSION_LABEL}")
endif()

include_directories(include)
# Add cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")
message(STATUS "tudatpy version: ${tudatpy_VERSION}")

if (_ENABLE_EXTENDED_ALIGNED_STORAGE)
    add_definitions("-D_ENABLE_EXTENDED_ALIGNED_STORAGE")
endif ()

if (MINGW)
    # https://github.com/tudat-team/tudatpy-feedstock/issues/3
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_hypot=hypot")
    # https://stackoverflow.com/questions/16596876/object-file-has-too-many-sections
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
    # https://sourceforge.net/p/mingw-w64/mailman/message/33154210/
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--large-address-aware")
endif ()

if (MSVC)
    # https://docs.microsoft.com/en-us/cpp/build/reference/bigobj-increase-number-of-sections-in-dot-obj-file?view=vs-2019
    # This is equivalent to the big-obj for MinGW, this should stop error C1128.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif ()

message(STATUS "=============================================================")
message(STATUS "******************** BUILD CONFIGURATION ********************")
message(STATUS "=============================================================")
message(STATUS "PROJECT_NAME ........................................ ${PROJECT_NAME}")
message(STATUS "CMAKE_BUILD_TYPE .................................... ${CMAKE_BUILD_TYPE}")
message(STATUS "BUILD_STYLE ......................................... ${BUILD_STYLE}")

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_STATIC_RUNTIME OFF)

if (NOT TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    add_definitions(-DTUDAT_BUILD_WITH_ESTIMATION_TOOLS=0)
else ()
    add_definitions(-DTUDAT_BUILD_WITH_ESTIMATION_TOOLS=1)
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Find Boost libraries on local system.
 find_package(Boost 1.72.0
        COMPONENTS thread date_time system unit_test_framework filesystem regex serialization REQUIRED)

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")

if (TUDAT_INSTALL)
    add_definitions(-DTUDAT_INSTALL=1)
endif ()

if (TUDAT_CONDA_INSTALL)
    add_definitions(-DTUDAT_CONDA_INSTALL=1)
endif ()

# Core Tudat dependency.
# Sofa dependency if in build settings.
if (NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(Tudat_DIR ${CMAKE_BINARY_DIR}/tudat)
endif ()
set (_TUDATPY_MIN_TUDAT_VERSION 2.9)
find_package(Tudat REQUIRED)
if(${Tudat_VERSION} VERSION_LESS ${_TUDATPY_MIN_TUDAT_VERSION})
    message(FATAL_ERROR "The minimum tudat version required by tudatpy is ${_TUDATPY_MIN_TUDAT_VERSION}, but version ${Tudat_VERSION} was found instead.")
endif()

message(STATUS "=============================================================")
message(STATUS "***************** TUDAT LINK LIBRARIES *********************")
message(STATUS "=============================================================")
message(STATUS Tudat_VERSION ........................................ "${Tudat_VERSION}")
message(STATUS Tudat_INCLUDE_DIRS ................................... "${Tudat_INCLUDE_DIRS}")
message(STATUS Tudat_LIBRARY_DIRS ................................... "${Tudat_LIBRARY_DIRS}")
message(STATUS Tudat_DATA_DIRS ...................................... "${Tudat_DATA_DIRS}")
message(STATUS "Tudat_PROPAGATION_LIBRARIES  ........................ [LISTED BELOW] ")
foreach(line IN LISTS Tudat_PROPAGATION_LIBRARIES)
    message(STATUS "                                                    - ${line}")
endforeach()
message(STATUS "Tudat_ESTIMATION_LIBRARIES  ......................... [LISTED BELOW] ")
foreach(line IN LISTS Tudat_ESTIMATION_LIBRARIES)
    message(STATUS "                                                    - ${line}")
endforeach()
message(STATUS "Tudat_INCLUDE_DIRS  ................................. [LISTED BELOW] ")
foreach(line IN LISTS Tudat_INCLUDE_DIRS)
    message(STATUS "                                                    - ${line}")
endforeach()

set(POSITION_INDEPENDENT_CODE True)

# python.
include(YACMAPythonSetup)

# python version check.
if (${PYTHON_VERSION_MAJOR} LESS 3 OR (${PYTHON_VERSION_MAJOR} EQUAL 3 AND ${PYTHON_VERSION_MINOR} LESS 4))
    message(FATAL_ERROR "Minimum supported python version is 3.4.")
endif ()

# pybind11.
find_package(pybind11 REQUIRED)

# Configure the sphinx config file.
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py" @ONLY)

# Add the module directory.
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/docs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_subdirectory(tudatpy)
