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

#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

if (_ENABLE_EXTENDED_ALIGNED_STORAGE)
    add_definitions("-D_ENABLE_EXTENDED_ALIGNED_STORAGE")
endif ()

# `PyLong_Type' can not be used when making a shared object; recompile with -fPIC
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

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
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif ()
message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

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
        COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")

# TODO: Make Tudat export definitions to the config for inheritence to this project.
add_definitions(-DTUDAT_BUILD_WITH_SPICE_INTERFACE=1)

#if (DEFINED _GLIBCXX_USE_CXX11_ABI)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=${_GLIBCXX_USE_CXX11_ABI}")
#endif ()
#set(TUDAT_BUNDLE_BUILD ON)
#find_package(Python3 COMPONENTS Interpreter Development)

if (TUDAT_INSTALL)
    add_definitions(-DTUDAT_INSTALL=1)
endif ()

if (TUDAT_CONDA_INSTALL)
    add_definitions(-DTUDAT_CONDA_INSTALL=1)
endif ()

if (TUDAT_TEST_INSTALL)
    add_definitions(-DTUDAT_TEST_INSTALL=1)
endif ()

# CSpice dependency.
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    find_package(CSpice REQUIRED 1.0.0)
else ()
    set(CSpice_INCLUDE_DIRS ${CSpice_BUILD_INCLUDE_DIRS})
endif ()


# Sofa dependency.
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    find_package(Sofa 1.0.0 REQUIRED)
else ()
    set(Sofa_INCLUDE_DIRS ${Sofa_BUILD_INCLUDE_DIRS})
endif ()

# nrlmsise-00 dependency.
#if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
find_package(NRLMSISE00 0.1 REQUIRED)
#else ()
#    set(NRLMSISE00_INCLUDE_DIRS ${NRLMSISE00_BUILD_INCLUDE_DIRS})
#endif ()

# Sofa dependency if in build settings.
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    find_package(Tudat REQUIRED)
    # TODO: Fi   gure out when and how to increment minimum Tudat version
    add_definitions(${Tudat_DEFINITIONS})
else ()
    set(Tudat_INCLUDE_DIRS ${Tudat_BUILD_INCLUDE_DIRS})
endif ()

#endif ()
message(STATUS Tudat: [Tudat_PROPAGATION_LIBRARIES] ${Tudat_PROPAGATION_LIBRARIES})
message(STATUS Tudat: [Tudat_INCLUDE_DIRS] ${Tudat_INCLUDE_DIRS})

set(POSITION_INDEPENDENT_CODE True)

enable_testing()

# python.
include(YACMAPythonSetup)

# python version check.
if (${PYTHON_VERSION_MAJOR} LESS 3 OR (${PYTHON_VERSION_MAJOR} EQUAL 3 AND ${PYTHON_VERSION_MINOR} LESS 4))
    message(FATAL_ERROR "Minimum supported python version is 3.4.")
endif ()

# pybind11.
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    find_package(pybind11 REQUIRED)
else ()
    set(pybind11_INCLUDE_DIRS ${pybind11_BUILD_INCLUDE_DIRS})
    set(pybind11_INCLUDE_DIR ${pybind11_BUILD_INCLUDE_DIRS})
endif ()

# Configure the sphinx config file.
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py" @ONLY)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

# Add the module directory.
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/docs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_subdirectory(tudatpy)
