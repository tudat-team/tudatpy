cmake_minimum_required(VERSION 3.20)

# Display basic cmake configuration
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")

# Remove caches if present [Trick to recompile changes for BUILD-TREE]
# From tudat cmake, might not be needed.
if (EXISTS "${CMAKE_BINARY_DIR}/CMakeCache.txt")
    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeCache.txt")
endif ()


# Policies
cmake_policy(SET CMP0069 NEW)  # Missing documentation
cmake_policy(SET CMP0048 NEW)  # Missing documentation

# Get version from file
file(STRINGS "version" VERSION_LABEL)
list(GET VERSION_LABEL 0 VERSION_LABEL)
string(
    REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)"
    TUDATPY_VERSION "${VERSION_LABEL}"
)

# Project declaration.
project(tudatpy VERSION "${TUDATPY_VERSION}" LANGUAGES CXX C)
set(BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

# Set position independent code [From tudat cmake]
set(CMAKE_POSITION_INDEPENDENT_CODE True)

# Allow access to cmake modules
list(APPEND CMAKE_MODULE_PATH "${BASE_DIRECTORY}/cmake/tudat")
list(APPEND CMAKE_MODULE_PATH "${BASE_DIRECTORY}/cmake/tudat/yolo")
list(APPEND CMAKE_MODULE_PATH "${BASE_DIRECTORY}/cmake/tudatpy")

# Additional statements from original cmake files (might not be needed)
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
enable_testing()
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

# +============================================================================
# BUILD OPTIONS
#  Offer the user the choice of defining the build variation.
# +============================================================================
option(
    TUDAT_BUILD_TESTS
    "Build C++ tests for tudat."
    ON
)
option(
    TUDAT_DOWNLOAD_AND_BUILD_BOOST
    "Downloads and builds boost"
    OFF
)
option(
    TUDAT_BUILD_WITH_PROPAGATION_TESTS
    "Build tudat with propagation tests. (>30 s propagations - Total test time > 10 minutes.)"
    OFF
)
option(
    TUDAT_BUILD_WITH_ESTIMATION_TOOLS
    "Build tudat with estimation tools."
    ON
)
option(
    TUDAT_BUILD_TUDAT_TUTORIALS
    "Build C++ tutorials."
    ON
)
option(
    TUDAT_BUILD_STATIC_LIBRARY
    "Build tudat as a static library, instead of dynamic."
    ON
)
option(
    TUDAT_BUILD_WITH_FILTERS
    "Build Tudat with preliminary Kalman Filter functionality."
    OFF
)
option(
    TUDAT_BUILD_WITH_SOFA_INTERFACE
    "Build Tudat with sofa interface."
    ON
)
option(
    TUDAT_BUILD_WITH_JSON_INTERFACE
    "Build Tudat with json interface."
    OFF
)
option(
    TUDAT_BUILD_WITH_FFTW3
    "Build Tudat with FFTW3."
    OFF
)
option(
    TUDAT_BUILD_WITH_PAGMO
    "Build Tudat with pagmo."
    OFF
)
if(CMAKE_CXX_SIMULATE_ID MATCHES "MSVC")
    # Build extended precision propagation tools.
    option(
        TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS
        "Build tudat with extended precision propagation tools."
        OFF
    )
else()
    option(
        TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS
        "Build tudat with extended precision propagation tools."
        OFF
    )
endif()

# Build as part of a GitHub Actions workflow
# Option enables the use of ccache by MSVC
# see https://github.com/ccache/ccache/wiki/MS-Visual-Studio
option(
    TUDAT_BUILD_GITHUB_ACTIONS
    "Build as part of GitHub Actions workflow."
    OFF
)
if (TUDAT_BUILD_GITHUB_ACTIONS AND MSVC)
    find_program(ccache_exe ccache)
    if(ccache_exe)
        message(STATUS "CONFIGURING MSVC FOR CCACHE")
        file(COPY_FILE ${ccache_exe} ${CMAKE_BINARY_DIR}/cl.exe ONLY_IF_DIFFERENT)

        # By default Visual Studio generators will use /Zi which is not compatible
        # with ccache, so tell Visual Studio to use /Z7 instead.
        message(STATUS "Setting MSVC debug information format to 'Embedded'")
        set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<$<CONFIG:Debug,RelWithDebInfo>:Embedded>")

        set(CMAKE_VS_GLOBALS
        "CLToolExe=cl.exe"
        "CLToolPath=${CMAKE_BINARY_DIR}"
        "TrackFileAccess=false"
        "UseMultiToolTask=true"
        "DebugInformationFormat=OldStyle"
        )
    endif()
endif()

include(show_build_configuration)

# +============================================================================
# TUDAT DEFINITIONS (Based on user configuration)
# +============================================================================
set(Tudat_DEFINITIONS "")
set(
    Tudat_DEFINITIONS
    "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_FILTERS=${TUDAT_BUILD_WITH_FILTERS}"
)
set(
    Tudat_DEFINITIONS
    "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_SOFA_INTERFACE=${TUDAT_BUILD_WITH_SOFA_INTERFACE}"
)
set(
    Tudat_DEFINITIONS
    "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_FFTW3=${TUDAT_BUILD_WITH_FFTW3}"
)
set(
    Tudat_DEFINITIONS
    "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_JSON_INTERFACE=${TUDAT_BUILD_WITH_JSON_INTERFACE}"
)
set(
    Tudat_DEFINITIONS
    "${Tudat_DEFINITIONS} -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=${TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS}"
)

# TODO: NOT SURE IF THIS BELONGS HERE (Comes from tudatpy cmake)
if (_ENABLE_EXTENDED_ALIGNED_STORAGE)
    add_definitions("-D_ENABLE_EXTENDED_ALIGNED_STORAGE")
endif ()
if (NOT TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    add_definitions(-DTUDAT_BUILD_WITH_ESTIMATION_TOOLS=0)
else ()
    add_definitions(-DTUDAT_BUILD_WITH_ESTIMATION_TOOLS=1)
endif ()
if (NOT TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS)
    add_definitions(-DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=0)
else ()
    add_definitions(-DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=1)
endif ()

# +============================================================================
# INSTALL TREE CONFIGURATION (Project name independent)
#  Offer the user the choice of overriding the installation directories.
# +============================================================================
# TODO: Consider if this can work with <tudat/paths.hpp> going forward.

include(install_tree_configuration)

# +============================================================================
# DEPENDENCIES SPECIFIC TO BUILD OPTIONS
# +============================================================================
if (NOT TUDAT_BUILD_STATIC_LIBRARY)
    # https://stackoverflow.com/questions/30398238/cmake-rpath-not-working-could-not-find-shared-object-file
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

# Set required boost components.
set(
    _TUDAT_REQUIRED_BOOST_LIBS
    filesystem system regex date_time thread chrono atomic unit_test_framework
)
# if (TUDAT_BUILD_TESTS)
#     # Add the unit test framework, if needed.
#     list(APPEND _TUDAT_REQUIRED_BOOST_LIBS unit_test_framework)
# endif ()

# Download boost if option is used.
if (
    TUDAT_DOWNLOAD_AND_BUILD_BOOST
    OR (TUDAT_BOOST_VERSION AND TUDAT_BOOST_SHA1)
)
    set(TUDAT_DOWNLOAD_AND_BUILD_BOOST ON)
    # Set the required boost components (for tudat)
    set(BoostComponents ${_TUDAT_REQUIRED_BOOST_LIBS})
    # Pick a boost version
    if (TUDAT_BOOST_VERSION AND TUDAT_BOOST_SHA1)
        set(BoostVersion ${TUDAT_BOOST_VERSION})
        set(BoostSHA1 ${TUDAT_BOOST_SHA1})
    else ()
        set(BoostVersion 1.67.0)
        set(BoostSHA1 694ae3f4f899d1a80eb7a3b31b33be73c423c1ae)
    endif ()
    message(STATUS "Downloading Boost Version ${BoostVersion}")
    include(cmake_modules/add_boost.cmake)
    message(STATUS "Downloading Boost Version ${BoostVersion} - done")
    #    set(BOOST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/boost)
    #    set(Boost_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/boost/stage/lib)
    #    set(Boost_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/boost)
    #    set(Boost_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/boost)
endif ()

# Set cmake build-type.
# If not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file), the default built type is "Release".
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_BUILD_TYPE Release)
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE Debug)
endif ()

# Set compiler based on user input
include(compiler)

# Find Boost and check GLBICXX_USE_CXX11_ABI (TODO: Better description)
include(TudatFindBoost)
message(STATUS "Checking for _GLIBCXX_USE_CXX11_ABI definition...")
if (DEFINED _GLIBCXX_USE_CXX11_ABI)
    if (_GLIBCXX_USE_CXX11_ABI)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=1")
        set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=1")
        message(STATUS "-D_GLIBCXX_USE_CXX11_ABI=1 added to CMAKE_CXX_FLAGS!")
    elseif (NOT _GLIBCXX_USE_CXX11_ABI)
        message(Boost_VERSION_MAJOR: ${Boost_VERSION_MAJOR})
        message(Boost_VERSION_MINOR: ${Boost_VERSION_MINOR})
        message(TUDAT_BUILD_GNU: ${TUDAT_BUILD_GNU})
        if (
            (${Boost_VERSION_MAJOR} LESS_EQUAL "1")
            AND (${Boost_VERSION_MINOR} LESS_EQUAL "67") AND ${TUDAT_BUILD_GNU}
        )
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
            set(
                Tudat_DEFINITIONS
                "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0"
            )
            message(
                STATUS
                "-D_GLIBCXX_USE_CXX11_ABI=0 added to CMAKE_CXX_FLAGS!"
            )
        elseif (
            (${Boost_VERSION_MAJOR} GREATER_EQUAL "1")
            AND (${Boost_VERSION_MINOR} GREATER_EQUAL "70")
            AND ${TUDAT_BUILD_GNU}
        )
            #            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
            #            set(Tudat_DEFINITIONS "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0")
            #            message(STATUS "-D_GLIBCXX_USE_CXX11_ABI=0 added to CMAKE_CXX_FLAGS!")
        elseif (${TUDAT_BUILD_MSVC})
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
            set(
                Tudat_DEFINITIONS
                "${Tudat_DEFINITIONS} -D_GLIBCXX_USE_CXX11_ABI=0"
            )
            message(
                STATUS
                "-D_GLIBCXX_USE_CXX11_ABI=0 added to CMAKE_CXX_FLAGS!"
            )
        else ()
            message(STATUS "_GLIBCXX_USE_CXX11_ABI was set as 0, but since the boost version is >=1.68|>1.70 [${Boost_VERSION}] and compiler_id is GNU, it's being ignored!")
            add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
        endif ()
    endif ()
else ()
    message(STATUS "_GLIBCXX_USE_CXX11_ABI was not found.")
endif ()

if (UNIX OR (WIN32 AND CLANG))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif ()

# Get tudat resource manager dependency.
find_package(TudatResources REQUIRED)
message(STATUS "TudatResources_INCLUDE_DIRS: ${TudatResources_INCLUDE_DIRS}")
include_directories(TudatResources_INCLUDE_DIRS)

# Find Eigen3 library on local system and include directories.
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
if (APPLE OR (WIN32 AND CLANG))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif ()

# PAGMO dependency
if (TUDAT_BUILD_WITH_PAGMO)
    add_definitions(-DTUDAT_BUILD_WITH_PAGMO=1)
    find_package(Pagmo REQUIRED)
else ()
    add_definitions(-DTUDAT_BUILD_WITH_PAGMO=0)
endif ()

# CSPICE dependency
add_definitions(-DTUDAT_BUILD_WITH_SPICE_INTERFACE=1)
find_package(CSpice REQUIRED 1.0.0)

# SOFA dependency
if (TUDAT_BUILD_WITH_SOFA_INTERFACE)
    find_package(Sofa REQUIRED 1.0.0)
    add_definitions(-DTUDAT_BUILD_WITH_SOFA_INTERFACE=1)
endif ()

# NRLMSISE-00 dependency
find_package(nrlmsise00 REQUIRED 0.1)
message(STATUS ${NRLMSISE00_LIBRARIES})
add_definitions(-DTUDAT_BUILD_WITH_NRLMSISE=1)

# FFTW3 dependency
if (TUDAT_BUILD_WITH_FFTW3)
    if(NOT APPLE)
        find_package(FFTW3)
        message(STATUS "FFTW " ${FFTW3_LIBRARIES})
        include_directories(SYSTEM AFTER "${FFTW3_INCLUDE_DIRS}")
        if(WIN32)
            set (FFTW3_LIBRARIES "${FFTW3_LIBRARY_DIRS}/fftw3.lib")
        else()
            set (FFTW3_LIBRARIES "${FFTW3_LIBRARY_DIRS}/libfftw3.so")
        endif()
        add_definitions(-DTUDAT_BUILD_WITH_FFTW3=1)
    else( )
        find_package(FFTW3)
        message(STATUS ${FFTW3_LIBRARIES})
        include_directories(SYSTEM AFTER "${FFTW3_INCLUDE_DIRS}")
        set(
            CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -isystem \"${FFTW3_INCLUDE_DIR}\""
        )
        set (FFTW3_LIBRARIES "${FFTW3_LIBRARY_DIRS}/libfftw3.dylib")
        add_definitions(-DTUDAT_BUILD_WITH_FFTW3=1)
    endif( )
else( )
    add_definitions(-DTUDAT_BUILD_WITH_FFTW3=0)
endif( )

# JSON dependency
if (TUDAT_BUILD_WITH_JSON_INTERFACE)
    find_package(nlohmann_json REQUIRED 3.7.3)
endif ()

# TODO: Not sure if this is needed. Comes from CMake of tudatpy
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_STATIC_RUNTIME OFF)

# Include YOLO functionality
# TODO: Check what this does and if we need it
include(YOLOProjectAddTestCase)
include(YOLOProjectAddLibrary)
include(YOLOProjectAddExecutable)
include(YOLOProjectAddExternalData)

# Define list of libraries in tudat and run setup for submodules
include(TudatLinkLibraries)
add_subdirectory(src/tudat/astro)
add_subdirectory(src/tudat/basics)
add_subdirectory(src/tudat/interface)
add_subdirectory(src/tudat/math)
add_subdirectory(src/tudat/simulation)
add_subdirectory(src/tudat/io)
add_subdirectory(src/tudat/utils)
if (TUDAT_BUILD_TESTS)
    add_subdirectory(tests/tudat)
endif ()

# Additional contents from original tudat cmake
# ========================================================================

file(
    RELATIVE_PATH TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX
    "${CMAKE_INSTALL_PREFIX}"
    "${INSTALL_DATA_DIR}/tudat"
)

# Set up config.hpp file of tudat for current version
configure_file(
        "${BASE_DIRECTORY}/configuration/tudat/config.hpp.in"
        "${PROJECT_BINARY_DIR}/include/tudat/config.hpp" @ONLY
)

# Installation of the header files.
install(
    DIRECTORY include/
    DESTINATION "${INSTALL_INCLUDE_DIR}"
)
install(
    FILES "${PROJECT_BINARY_DIR}/include/tudat/config.hpp"
    DESTINATION ${INSTALL_INCLUDE_DIR}/tudat
)

# Add all targets to the build-tree export set
export(
    TARGETS ${Tudat_LIBRARIES}
    NAMESPACE Tudat::
    FILE "${CMAKE_CURRENT_BINARY_DIR}/tudat_export.cmake"
)

# Installation of the export modules.
install(EXPORT tudat_export NAMESPACE Tudat:: DESTINATION "${RELATIVE_INSTALL_CMAKE_DIR}")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE Tudat)

# +============================================================================
# PROJECT CONFIG
#  For external use of project.
# +============================================================================
include(CMakePackageConfigHelpers)

configure_package_config_file(
        "${BASE_DIRECTORY}/configuration/tudat/tudat-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/tudat-config.cmake"
        INSTALL_DESTINATION "${INSTALL_CMAKE_DIR}"
)

configure_file(
        "${BASE_DIRECTORY}/configuration/tudat/tudat-config.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/tudat-config.cmake" @ONLY
)

# Configure file path for tudat data loading.
#set(Tudat_DATA_DIRS ${INSTALL_CMAKE_DIR}/${RELATIVE_INSTALL_DATA_DIR})


file(
        COPY "${BASE_DIRECTORY}/cmake/tudat/TudatFindBoost.cmake"
        DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
)

install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/tudat-config.cmake"
        "${BASE_DIRECTORY}/cmake/tudat/TudatFindBoost.cmake"
        DESTINATION "${INSTALL_CMAKE_DIR}"
)

# Take care of versioning.

# NOTE: SameMinorVersion available only
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/tudat-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
)
install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/tudat-config-version.cmake"
    DESTINATION "${INSTALL_CMAKE_DIR}"
)

# ========
message(
    STATUS
    "TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX: ${TUDAT_DATA_DIR_RELATIVE_TO_INSTALL_PREFIX}"
)

# +============================================================================
# Tudatpy specific configuration
# +============================================================================

message(STATUS Tudat: [Tudat_PROPAGATION_LIBRARIES] ${Tudat_PROPAGATION_LIBRARIES})
message(STATUS Tudat: [Tudat_INCLUDE_DIRS] ${Tudat_INCLUDE_DIRS})

# Pybind11 dependency
include(YACMAPythonSetup)
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
set(pybind11_INCLUDE_DIRS ${pybind11_BUILD_INCLUDE_DIRS}) # TODO: Check
set(pybind11_INCLUDE_DIR ${pybind11_BUILD_INCLUDE_DIRS})  # TODO: Check

# Source directory for tudatpy
set(TUDATPY_SOURCE_DIR ${BASE_DIRECTORY}/src/tudatpy)

# Add docs to build directory
file(
    COPY ${BASE_DIRECTORY}/docs/tudatpy
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/docs
)

# Include tudatpy CMake configuration
include_directories(include/tudatpy)
add_subdirectory(src/tudatpy)

# +============================================================================
# CLEAN UP (Project name independent)
#  Cleanup in case project is not top level.
# +============================================================================
unset(INSTALL_LIB_DIR)
unset(INSTALL_BIN_DIR)
unset(INSTALL_INCLUDE_DIR)
unset(INSTALL_DATA_DIR)
unset(INSTALL_MAN_DIR)
unset(INSTALL_TESTS_DIR)
unset(INSTALL_CMAKE_DIR)

unset(RELATIVE_INSTALL_LIB_DIR)
unset(RELATIVE_INSTALL_BIN_DIR)
unset(RELATIVE_INSTALL_INCLUDE_DIR)
unset(RELATIVE_INSTALL_DATA_DIR)
unset(RELATIVE_INSTALL_MAN_DIR)
unset(RELATIVE_INSTALL_TESTS_DIR)
unset(RELATIVE_INSTALL_CMAKE_DIR)
