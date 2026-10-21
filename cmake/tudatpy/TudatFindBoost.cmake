# Tudat only uses header-only Boost libraries.
set(_TUDAT_BOOST_MINIMUM_VERSION 1.84.0)
set(Boost_NO_WARN_NEW_VERSIONS ON)
find_package(Boost ${_TUDAT_BOOST_MINIMUM_VERSION} QUIET REQUIRED)

message(STATUS "Detected Boost version: ${Boost_VERSION}")
message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")

# Recreate the header target if an older FindBoost module did not provide it.
if(NOT TARGET Boost::boost)
    message(STATUS "The 'Boost::boost' target is missing, creating it.")
    add_library(Boost::boost INTERFACE IMPORTED)
    set_target_properties(Boost::boost PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
            INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")
endif()

if(NOT TARGET Boost::disable_autolinking)
    message(STATUS "The 'Boost::disable_autolinking' target is missing, creating it.")
    add_library(Boost::disable_autolinking INTERFACE IMPORTED)
    if(WIN32)
        set_target_properties(Boost::disable_autolinking PROPERTIES
                INTERFACE_COMPILE_DEFINITIONS "BOOST_ALL_NO_LIB")
    endif()
endif()

unset(_TUDAT_BOOST_MINIMUM_VERSION)
