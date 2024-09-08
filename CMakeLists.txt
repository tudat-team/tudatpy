cmake_minimum_required(VERSION 3.15)

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

# Get version + label from file
file(STRINGS "version" tudatpy_VERSION_LABEL)
list(GET tudatpy_VERSION_LABEL 0 tudatpy_VERSION_LABEL)
string(
    REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)"
    tudatpy_VERSION "${tudatpy_VERSION_LABEL}"
)

# Project declaration.
project(tudatpy VERSION "${tudatpy_VERSION}" LANGUAGES CXX C)
set(TUDATPY_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# Add cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

# Initial log
message("")
message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${PROJECT_VERSION}) >>")
message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")
message("")


# Additional compiler configuration
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
endif ()

# Build configuration
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
add_definitions(-DTUDAT_BUILD_WITH_SPICE_INTERFACE=1)

# Possibly deprecated
if (TUDAT_INSTALL)
    add_definitions(-DTUDAT_INSTALL=1)
endif ()
if (TUDAT_CONDA_INSTALL)
    add_definitions(-DTUDAT_CONDA_INSTALL=1)
endif ()
if (TUDAT_TEST_INSTALL)
    add_definitions(-DTUDAT_TEST_INSTALL=1)
endif ()

# Dependencies
######################################

# Boost
find_package(
    Boost 1.78.0
    COMPONENTS thread date_time system unit_test_framework filesystem regex
    REQUIRED
)

# Eigen
# Second line allows to include <Eigen/...> instead of <eigen3/Eigen/...>
find_package(Eigen3 REQUIRED)
include_directories(AFTER SYSTEM ${EIGEN3_INCLUDE_DIR})

# SPICE
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

# NRLMSISE-00
find_package(NRLMSISE00 0.1 REQUIRED)

# Tudat
if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    find_package(Tudat REQUIRED)
    add_definitions(${Tudat_DEFINITIONS})
else ()
    set(Tudat_INCLUDE_DIRS ${Tudat_BUILD_INCLUDE_DIRS})
endif ()

# Python and Pybind11
# Find Python and pybind11.
if (NOT pybind11_FOUND)
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG REQUIRED)
else()
    set(pybind11_INCLUDE_DIRS ${pybind11_BUILD_INCLUDE_DIRS})
endif()

################################################

# Misc
set(POSITION_INDEPENDENT_CODE True)
enable_testing()

# Docs
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/docs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Tudatpy subdirectory
include_directories(include)
add_link_options("-L${CMAKE_INSTALL_PREFIX}/lib")
add_subdirectory(src/tudatpy)
