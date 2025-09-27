# FindCSpice.cmake - Attempt to find CSpice libraries and include paths
# This module defines:
#   CSpice_FOUND - True if headers and requested libraries were found
#   CSpice_INCLUDE_DIRS - Where to find the headers
#   CSpice_LIBRARIES - List of libraries when using CSpice
#   CSpice_VERSION - The version of the found CSpice
set(SEARCH_PATHS ${CMAKE_PREFIX_PATH} ${CMAKE_INSTALL_PREFIX} $ENV{CONDA_PREFIX})
# Create a new list for suffixed paths
set(SUFFIXED_PATHS)

# Duplicate paths with 'library' suffix and add both original and suffixed to SEARCH_PATHS
foreach(PATH IN LISTS SEARCH_PATHS)
    list(APPEND SUFFIXED_PATHS "${PATH}/Library")
endforeach()

if(NOT TARGET CSpice::cspice)
    # Define the version of CSpice you are looking for
    set(CSpice_VERSION "67.0.0")

    # Search for the include directory containing cspice/SpiceUsr.h
    find_path(CSpice_INCLUDE_DIR_PARENT
            NAMES cspice/SpiceUsr.h
            PATHS ${SEARCH_PATHS} ${SUFFIXED_PATHS}
            PATH_SUFFIXES include
            DOC "Parent directory where cspice/SpiceUsr.h can be found"
    )

    # Adjust CSpice_INCLUDE_DIRS to the parent of the found include directory
    if(CSpice_INCLUDE_DIR_PARENT)
        set(CSpice_INCLUDE_DIRS "${CSpice_INCLUDE_DIR_PARENT}")
    endif()

    # Search for the main CSpice library
    find_library(CSpice_LIBRARY
            NAMES cspice
            PATHS ${SEARCH_PATHS} ${SUFFIXED_PATHS}
            PATH_SUFFIXES lib
            DOC "Main CSpice library"
    )

    # Search for the CSpice support library
    find_library(CSpice_SUPPORT_LIBRARY
            NAMES csupport.66 csupport.67 csupport csupport.a
            PATHS ${SEARCH_PATHS} ${SUFFIXED_PATHS}
            PATH_SUFFIXES lib
            DOC "CSpice support library"
    )

    # Aggregate found components
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(CSpice
            REQUIRED_VARS CSpice_LIBRARY CSpice_INCLUDE_DIRS
            VERSION_VAR CSpice_VERSION
    )

    if(CSpice_FOUND)
        # Create imported targets for CSpice
        add_library(CSpice::cspice UNKNOWN IMPORTED)
        set_target_properties(CSpice::cspice PROPERTIES
                IMPORTED_LOCATION "${CSpice_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${CSpice_INCLUDE_DIRS}"
        )

        add_library(CSpice::csupport UNKNOWN IMPORTED)
        set_target_properties(CSpice::csupport PROPERTIES
                IMPORTED_LOCATION "${CSpice_SUPPORT_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${CSpice_INCLUDE_DIRS}"
        )

        # Aggregate the libraries to a variable for external usage
        set(CSpice_LIBRARIES CSpice::cspice CSpice::csupport)

        message(STATUS "Found CSpice: ${CSpice_INCLUDE_DIRS} (found version ${CSpice_VERSION})")
    else()
        message(STATUS "Could not find CSpice")
    endif()

    mark_as_advanced(CSpice_INCLUDE_DIRS CSpice_LIBRARY CSpice_SUPPORT_LIBRARY CSpice_LIBRARIES CSpice_VERSION)
endif()
