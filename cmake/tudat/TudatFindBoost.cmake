
# Set CMake policies to suppress warnings
cmake_policy(SET CMP0144 NEW)  # find_package uses upper-case <PACKAGENAME>_ROOT variables
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)  # Use FindBoost module (not yet removed)
endif()

# Run a first pass for finding the headers only,
# and establishing the Boost version.
set(_TUDAT_BOOST_MINIMUM_VERSION 1.72.0)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

# For Emscripten/WASM builds, we only use Boost headers (no compiled libraries)
# The Emscripten port provides headers only, and most Boost functionality
# Tudat needs can work header-only or with Emscripten's threading model
if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten" OR EMSCRIPTEN OR TUDAT_BUILD_WASM)
    message(STATUS "WASM build: Using Boost headers only")
    # Lower minimum version for WASM since we only need headers
    find_package(Boost 1.67.0 QUIET REQUIRED)
    # Skip looking for compiled Boost libraries - set them as found
    set(Boost_FILESYSTEM_FOUND TRUE)
    set(Boost_SYSTEM_FOUND TRUE)
    set(Boost_REGEX_FOUND TRUE)
    set(Boost_DATE_TIME_FOUND TRUE)
    set(Boost_THREAD_FOUND TRUE)
    set(Boost_CHRONO_FOUND TRUE)
    set(Boost_ATOMIC_FOUND TRUE)
    # Create interface-only targets for the "libraries"
    foreach(_comp filesystem system regex date_time thread chrono atomic)
        if(NOT TARGET Boost::${_comp})
            add_library(Boost::${_comp} INTERFACE IMPORTED)
            set_target_properties(Boost::${_comp} PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")
        endif()
    endforeach()
    set(_TUDAT_REQUIRED_BOOST_LIBS "")
else()
    find_package(Boost ${_TUDAT_BOOST_MINIMUM_VERSION} QUIET REQUIRED)
endif()

if(_TUDAT_FIND_BOOST_PYTHON)
    # NOTE: since Boost 1.67, the naming of the Boost.Python library has changed to include the
    # major and minor python version as a suffix. See the release notes:
    # https://www.boost.org/users/history/version_1_67_0.html
    if(${Boost_MAJOR_VERSION} GREATER 1 OR (${Boost_MAJOR_VERSION} EQUAL 1 AND ${Boost_MINOR_VERSION} GREATER 66))
        list(APPEND _TUDAT_REQUIRED_BOOST_LIBS "python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
    else()
        list(APPEND _TUDAT_REQUIRED_BOOST_LIBS python3)
    endif()
endif()
message(STATUS "Required Boost libraries: ${_TUDAT_REQUIRED_BOOST_LIBS}")
# For WASM builds, we already found Boost headers and created interface targets
if(NOT (CMAKE_SYSTEM_NAME STREQUAL "Emscripten" OR EMSCRIPTEN OR TUDAT_BUILD_WASM))
    find_package(Boost ${_TUDAT_BOOST_MINIMUM_VERSION} REQUIRED COMPONENTS ${_TUDAT_REQUIRED_BOOST_LIBS})
    if(NOT Boost_FOUND)
        message(FATAL_ERROR "Not all requested Boost components were found, exiting.")
    endif()
endif()
message(STATUS "Detected Boost version: ${Boost_VERSION}")
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
# Might need to recreate targets if they are missing (e.g., older cmake versions).
if(NOT TARGET Boost::boost)
    message(STATUS "The 'Boost::boost' target is missing, creating it.")
    add_library(Boost::boost INTERFACE IMPORTED)
    set_target_properties(Boost::boost PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")
endif()
if(NOT TARGET Boost::disable_autolinking)
    message(STATUS "The 'Boost::disable_autolinking' target is missing, creating it.")
    add_library(Boost::disable_autolinking INTERFACE IMPORTED)
    if(WIN32)
        set_target_properties(Boost::disable_autolinking PROPERTIES INTERFACE_COMPILE_DEFINITIONS "BOOST_ALL_NO_LIB")
    endif()
endif()
foreach(_TUDAT_BOOST_COMPONENT ${_TUDAT_REQUIRED_BOOST_LIBS})
    if(NOT TARGET Boost::${_TUDAT_BOOST_COMPONENT})
        message(STATUS "The 'Boost::${_TUDAT_BOOST_COMPONENT}' imported target is missing, creating it.")
        string(TOUPPER ${_TUDAT_BOOST_COMPONENT} _TUDAT_BOOST_UPPER_COMPONENT)
        if(Boost_USE_STATIC_LIBS)
            add_library(Boost::${_TUDAT_BOOST_COMPONENT} STATIC IMPORTED)
        else()
            add_library(Boost::${_TUDAT_BOOST_COMPONENT} UNKNOWN IMPORTED)
        endif()
        set_target_properties(Boost::${_TUDAT_BOOST_COMPONENT} PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")
        set_target_properties(Boost::${_TUDAT_BOOST_COMPONENT} PROPERTIES
                IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
                IMPORTED_LOCATION "${Boost_${_TUDAT_BOOST_UPPER_COMPONENT}_LIBRARY}")
    endif()
endforeach()

unset(_TUDAT_BOOST_MINIMUM_VERSION)
unset(_TUDAT_REQUIRED_BOOST_LIBS)
