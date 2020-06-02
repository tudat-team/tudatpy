cmake_minimum_required(VERSION 2.8.3)

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

#CMP0048 is set to NEW
project(TudatPy VERSION 0.2.0 LANGUAGES CXX C)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/cmake")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}-Wl,--no-undefined")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")
message(STATUS "tudatpy version: ${tudatpy_VERSION}")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# `PyLong_Type' can not be used when making a shared object; recompile with -fPIC
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")

# Add local module path
#list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules")
message(STATUS "CMake Module path(s): ${CMAKE_MODULE_PATH}")

# Set compiler based on preferences (e.g. USE_CLANG) and system.
#include(compiler)

# Set testing options based on platform.
#enable_testing()

# Find Eigen3 library on local system.
#find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set CMake flag to suppress Eigen warnings (platform-dependent solution).
#if (NOT APPLE OR APPLE_INCLUDE_FORCE)
#    include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
#else ()
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
#endif ()

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

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


#message(STATUS Tudat: [Tudat_INCLUDE_DIRS] ${Tudat_LIBRARIES})

#set(BOOST_ROOT /home/ggarrett/anaconda3/envs/test)

# Find Boost libraries on local system.
find_package(Boost 1.64.0
        COMPONENTS thread date_time system unit_test_framework filesystem regex REQUIRED)

#add_subdirectory(eigen)
#if (NOT EIGEN3_INCLUDE_DIR)
find_package(Eigen3 REQUIRED)
#endif ()

#add_subdirectory(cspice)
find_package(CSpice REQUIRED 1.0.0)

#add_subdirectory(sofa)
find_package(Sofa REQUIRED 1.0.0)

# TODO: Make Tudat export definitions to the config for inheritence to this project.
#set(Tudat_DEFINITIONS "DTUDAT_BUILD_WITH_SPICE_INTERFACE=ON")
add_definitions(-DTUDAT_BUILD_WITH_SPICE_INTERFACE=1)
add_definitions(-DTUDAT_BUILD_WITH_SOFA_INTERFACE=1)
add_definitions(-DTUDAT_BUILD_WITH_JSON_INTERFACE=1)

# Find Tudat library on local system.
find_package(Tudat 4.0 REQUIRED)
if (NOT ${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    if (DEFINED ${Tudat_BUILD_INCLUDE_DIRS})
        set(Tudat_INCLUDE_DIRS ${Tudat_BUILD_INCLUDE_DIRS})
    endif ()
endif ()

#endif ()
message(STATUS Tudat: [Tudat_PROPAGATION_LIBRARIES] ${Tudat_PROPAGATION_LIBRARIES})
message(STATUS Tudat: [Tudat_INCLUDE_DIRS] ${Tudat_INCLUDE_DIRS})


#include(TudatFindBoost)


set(POSITION_INDEPENDENT_CODE True)

enable_testing()

# python.
include(YACMAPythonSetup)

# python version check.
if (${PYTHON_VERSION_MAJOR} LESS 3 OR (${PYTHON_VERSION_MAJOR} EQUAL 3 AND ${PYTHON_VERSION_MINOR} LESS 4))
    message(FATAL_ERROR "Minimum supported python version is 3.4.")
endif ()

# pybind11.
find_package(pybind11)

# Configure the sphinx config file.
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py" @ONLY)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

# Add the module directory.
add_subdirectory(tudatpy)
