# CMakeLists.txt for MCD (Mars Climate Database)
# All MCD-specific configuration is contained in this file

# Fortran should already be enabled by parent CMakeLists.txt

# Set Fortran compiler flags - based on working Makefile
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # Use only optimization, no default real type changes
    # The MCD code has explicit type declarations and mixing them causes issues
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -ffree-form -ffree-line-length-none")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -free")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
endif()

# Find NetCDF (required for MCD)
find_package(NetCDF REQUIRED COMPONENTS Fortran C)

if(NOT NetCDF_FOUND)
    message(FATAL_ERROR "NetCDF not found. MCD requires NetCDF Fortran and C libraries. "
                       "Please install NetCDF using: conda install -c conda-forge netcdf-fortran")
endif()

message(STATUS "MCD: NetCDF C library: ${NetCDF_C_LIBRARIES}")
message(STATUS "MCD: NetCDF Fortran library: ${NetCDF_Fortran_LIBRARIES}")
message(STATUS "MCD: NetCDF include dirs: ${NetCDF_INCLUDE_DIRS}")

# Check if MCD.F90 exists
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/MCD.F90)
    message(FATAL_ERROR "MCD.F90 not found at ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# Check if mcd.h exists
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mcd.h)
    message(WARNING "mcd.h not found at ${CMAKE_CURRENT_SOURCE_DIR}/mcd.h")
endif()

# Create MCD library with OBJECT library type first to avoid export issues
add_library(mcd_fortran_obj OBJECT
    ${CMAKE_CURRENT_SOURCE_DIR}/MCD.F90
)

# Set Fortran module directory for object library
set_target_properties(mcd_fortran_obj PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mcd_modules
    POSITION_INDEPENDENT_CODE ON
)

# Make modules accessible - use PRIVATE to avoid export issues
target_include_directories(mcd_fortran_obj PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}/mcd_modules
)

# Link NetCDF to object library
if(TARGET NetCDF::NetCDF_Fortran)
    target_link_libraries(mcd_fortran_obj PUBLIC
        NetCDF::NetCDF_Fortran
        NetCDF::NetCDF_C
    )
    
    # Also propagate include directories from NetCDF
    target_include_directories(mcd_fortran_obj PRIVATE
        ${NetCDF_INCLUDE_DIRS}
    )
else()
    message(FATAL_ERROR "NetCDF targets not found. MCD requires NetCDF to be found first.")
endif()

# Now create the actual static library from the object library
add_library(mcd_fortran STATIC $<TARGET_OBJECTS:mcd_fortran_obj>)

# Set properties for the final library
set_target_properties(mcd_fortran PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mcd_modules
    PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/mcd.h"
)

# NO include directories on mcd_fortran itself - handled at link time by consumers
# Link NetCDF libraries to the final library
target_link_libraries(mcd_fortran PUBLIC
    NetCDF::NetCDF_Fortran
    NetCDF::NetCDF_C
)

# Set MCD data directory path (absolute path to the data directory)
set(MCD_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")

# Validate that the data directory exists
if(NOT EXISTS "${MCD_DATA_DIR}")
    message(WARNING "MCD data directory does not exist: ${MCD_DATA_DIR}")
    message(WARNING "MCD atmosphere model will not work correctly without data files.")
else()
    message(STATUS "MCD data directory found: ${MCD_DATA_DIR}")
endif()

# Export MCD_DATA_DIR as a compile definition for use by C++ code
# This passes the absolute path to the data directory at compile time
target_compile_definitions(mcd_fortran PUBLIC
    MCD_DATA_PATH="${MCD_DATA_DIR}/"
)

# Also propagate this definition to any target that links against mcd_fortran
# This ensures tudat_aerodynamics gets the definition
set_target_properties(mcd_fortran PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "MCD_DATA_PATH=\"${MCD_DATA_DIR}/\""
)

message(STATUS "MCD_DATA_PATH set to: ${MCD_DATA_DIR}/")

# Install library - ADD to tudat_export to resolve dependency
install(TARGETS mcd_fortran
    EXPORT tudat_export
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include/tudat/external/mcd
)

# Install Fortran modules
install(DIRECTORY ${CMAKE_BINARY_DIR}/mcd_modules/
    DESTINATION include/tudat/external/mcd/modules
    FILES_MATCHING PATTERN "*.mod"
)

message(STATUS "MCD Fortran library configured:")
message(STATUS "  Library target: mcd_fortran")
message(STATUS "  Source: ${CMAKE_CURRENT_SOURCE_DIR}/MCD.F90")
message(STATUS "  Module dir: ${CMAKE_BINARY_DIR}/mcd_modules")
message(STATUS "  Header: ${CMAKE_CURRENT_SOURCE_DIR}/mcd.h")
message(STATUS "  Compiler flags: ${CMAKE_Fortran_FLAGS}")
