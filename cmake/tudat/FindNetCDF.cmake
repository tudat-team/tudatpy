# FindNetCDF.cmake - Find NetCDF library (C and Fortran components)
# This module defines:
#   NetCDF_FOUND - True if both C and Fortran components are found
#   NetCDF_INCLUDE_DIRS - Include directories
#   NetCDF_LIBRARIES - Libraries to link
#   NetCDF_C_LIBRARIES - C library
#   NetCDF_Fortran_LIBRARIES - Fortran library
#   NetCDF::NetCDF_C - Imported target for C library
#   NetCDF::NetCDF_Fortran - Imported target for Fortran library

# Prioritize conda environment if available
if(DEFINED ENV{CONDA_PREFIX})
    set(CONDA_PREFIX $ENV{CONDA_PREFIX})
    list(INSERT CMAKE_PREFIX_PATH 0 ${CONDA_PREFIX})
    message(STATUS "Using conda environment: ${CONDA_PREFIX}")
endif()

# Get search paths from environment and CMake
set(SEARCH_PATHS 
    ${CMAKE_PREFIX_PATH} 
    ${CMAKE_INSTALL_PREFIX} 
    $ENV{CONDA_PREFIX}
    $ENV{NetCDF_ROOT}
    $ENV{NETCDF_ROOT}
    $ENV{NetCDF_DIR}
    /usr
    /usr/local
)

# Add standard library suffixes
set(LIB_SUFFIXES lib lib64 lib/x86_64-linux-gnu)

# Add suffixed paths for Windows conda environments
if(WIN32)
    set(WIN_SUFFIXED_PATHS)
    foreach(PATH IN LISTS SEARCH_PATHS)
        list(APPEND WIN_SUFFIXED_PATHS "${PATH}/Library")
    endforeach()
    list(APPEND SEARCH_PATHS ${WIN_SUFFIXED_PATHS})
    list(APPEND LIB_SUFFIXES Library/lib)
endif()

# Find NetCDF C library
find_library(NetCDF_C_LIBRARY
    NAMES netcdf
    HINTS ${SEARCH_PATHS}
    PATH_SUFFIXES ${LIB_SUFFIXES}
    DOC "NetCDF C library"
    NO_DEFAULT_PATH
)

# If not found with NO_DEFAULT_PATH, try system paths
if(NOT NetCDF_C_LIBRARY)
    find_library(NetCDF_C_LIBRARY
        NAMES netcdf
        PATH_SUFFIXES ${LIB_SUFFIXES}
        DOC "NetCDF C library"
    )
endif()

# Find NetCDF Fortran library
find_library(NetCDF_Fortran_LIBRARY
    NAMES netcdff
    HINTS ${SEARCH_PATHS}
    PATH_SUFFIXES ${LIB_SUFFIXES}
    DOC "NetCDF Fortran library"
    NO_DEFAULT_PATH
)

# If not found with NO_DEFAULT_PATH, try system paths
if(NOT NetCDF_Fortran_LIBRARY)
    find_library(NetCDF_Fortran_LIBRARY
        NAMES netcdff
        PATH_SUFFIXES ${LIB_SUFFIXES}
        DOC "NetCDF Fortran library"
    )
endif()

# Find NetCDF include directory (contains netcdf.h)
find_path(NetCDF_INCLUDE_DIR
    NAMES netcdf.h
    HINTS ${SEARCH_PATHS}
    PATH_SUFFIXES include
    DOC "NetCDF C include directory"
    NO_DEFAULT_PATH
)

# If not found with NO_DEFAULT_PATH, try system paths
if(NOT NetCDF_INCLUDE_DIR)
    find_path(NetCDF_INCLUDE_DIR
        NAMES netcdf.h
        PATH_SUFFIXES include
        DOC "NetCDF C include directory"
    )
endif()

# Find NetCDF Fortran module directory (contains netcdf.mod)
# In conda, this is usually the same as the C include directory
find_path(NetCDF_Fortran_INCLUDE_DIR
    NAMES netcdf.mod
    HINTS ${SEARCH_PATHS}
    PATH_SUFFIXES include
    DOC "NetCDF Fortran module directory"
    NO_DEFAULT_PATH
)

# If not found with NO_DEFAULT_PATH, try system paths
if(NOT NetCDF_Fortran_INCLUDE_DIR)
    find_path(NetCDF_Fortran_INCLUDE_DIR
        NAMES netcdf.mod
        PATH_SUFFIXES include
        DOC "NetCDF Fortran module directory"
    )
endif()

# In conda environments, Fortran modules are often in the same directory as C headers
if(NOT NetCDF_Fortran_INCLUDE_DIR AND NetCDF_INCLUDE_DIR)
    if(EXISTS "${NetCDF_INCLUDE_DIR}/netcdf.mod")
        set(NetCDF_Fortran_INCLUDE_DIR ${NetCDF_INCLUDE_DIR})
        message(STATUS "Found NetCDF Fortran modules in C include directory")
    endif()
endif()

# Set variables if all required components are found
if(NetCDF_C_LIBRARY AND NetCDF_Fortran_LIBRARY AND NetCDF_INCLUDE_DIR)
    set(NetCDF_FOUND TRUE)
    set(NetCDF_LIBRARIES ${NetCDF_C_LIBRARY} ${NetCDF_Fortran_LIBRARY})
    set(NetCDF_C_LIBRARIES ${NetCDF_C_LIBRARY})
    set(NetCDF_Fortran_LIBRARIES ${NetCDF_Fortran_LIBRARY})
    
    # Set include dirs
    set(NetCDF_INCLUDE_DIRS ${NetCDF_INCLUDE_DIR})
    if(NetCDF_Fortran_INCLUDE_DIR AND NOT "${NetCDF_Fortran_INCLUDE_DIR}" STREQUAL "${NetCDF_INCLUDE_DIR}")
        list(APPEND NetCDF_INCLUDE_DIRS ${NetCDF_Fortran_INCLUDE_DIR})
    endif()
    
    # Remove duplicates
    list(REMOVE_DUPLICATES NetCDF_INCLUDE_DIRS)
    
    # Create imported target for NetCDF C
    if(NOT TARGET NetCDF::NetCDF_C)
        add_library(NetCDF::NetCDF_C UNKNOWN IMPORTED)
        set_target_properties(NetCDF::NetCDF_C PROPERTIES
            IMPORTED_LOCATION "${NetCDF_C_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIR}"
        )
        
        # Add HDF5 and other dependencies if needed (common in conda)
        if(DEFINED ENV{CONDA_PREFIX})
            # NetCDF typically depends on HDF5, curl, and zlib in conda
            find_library(HDF5_LIBRARY NAMES hdf5 HINTS ${CONDA_PREFIX} PATH_SUFFIXES ${LIB_SUFFIXES})
            find_library(HDF5_HL_LIBRARY NAMES hdf5_hl HINTS ${CONDA_PREFIX} PATH_SUFFIXES ${LIB_SUFFIXES})
            
            if(HDF5_LIBRARY)
                set_property(TARGET NetCDF::NetCDF_C APPEND PROPERTY 
                    INTERFACE_LINK_LIBRARIES ${HDF5_LIBRARY})
            endif()
            if(HDF5_HL_LIBRARY)
                set_property(TARGET NetCDF::NetCDF_C APPEND PROPERTY 
                    INTERFACE_LINK_LIBRARIES ${HDF5_HL_LIBRARY})
            endif()
        endif()
    endif()
    
    # Create imported target for NetCDF Fortran
    if(NOT TARGET NetCDF::NetCDF_Fortran)
        add_library(NetCDF::NetCDF_Fortran UNKNOWN IMPORTED)
        set_target_properties(NetCDF::NetCDF_Fortran PROPERTIES
            IMPORTED_LOCATION "${NetCDF_Fortran_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "NetCDF::NetCDF_C"
        )
    endif()
else()
    set(NetCDF_FOUND FALSE)
endif()

# Handle standard find_package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
    REQUIRED_VARS 
        NetCDF_C_LIBRARY 
        NetCDF_Fortran_LIBRARY 
        NetCDF_INCLUDE_DIR
    FOUND_VAR NetCDF_FOUND
    FAIL_MESSAGE "Could not find NetCDF. Please ensure netcdf-fortran is installed: conda install -c conda-forge netcdf-fortran"
)

# Mark variables as advanced
mark_as_advanced(
    NetCDF_C_LIBRARY
    NetCDF_Fortran_LIBRARY
    NetCDF_INCLUDE_DIR
    NetCDF_Fortran_INCLUDE_DIR
)

# Print detailed info when NetCDF is found
if(NetCDF_FOUND)
    message(STATUS "NetCDF found:")
    message(STATUS "  NetCDF C library: ${NetCDF_C_LIBRARY}")
    message(STATUS "  NetCDF Fortran library: ${NetCDF_Fortran_LIBRARY}")
    message(STATUS "  NetCDF include dirs: ${NetCDF_INCLUDE_DIRS}")
    if(DEFINED ENV{CONDA_PREFIX})
        message(STATUS "  Using conda environment")
    endif()
endif()