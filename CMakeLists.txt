#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.


# Specify minimum cmake version required to compile Tudat.
cmake_minimum_required(VERSION 3.3)

if (EXISTS "${CMAKE_BINARY_DIR}/CMakeCache.txt")
    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeCache.txt")
endif ()

# Set cmake policy.
cmake_policy(SET CMP0069 NEW)

# Specific Tudat project name.
project(Tudat VERSION 4.0.0 LANGUAGES C CXX)


set(CMAKE_BUILD_TYPE DEBUG)

# TODO: Fix temporary TUDAT_CUSTOM_ROOT_PATH fix.
set(TUDAT_CUSTOM_ROOT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src")

#set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" /home/ggarrett/anaconda3/envs/py37/lib/cmake)
set(CMAKE_LIBRARY_PATH "/home/ggarrett/anaconda3/envs/test/lib")
set(CMAKE_INCLUDE_PATH "/home/ggarrett/anaconda3/envs/test/include")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

# Add local module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")
enable_testing()

# Build option: enable the test suite.
option(TUDAT_BUILD_TESTS "Build the test suite." ON)

# Build option: include default data suite.
option(TUDAT_FETCH_DATA_SUITE "Downloads default data suite during build." OFF)

# Build extended precision propagation tools.
option(TUDAT_BUILD_WITH_PROPAGATION_TESTS "Build tudat with propagation tests. (>30 s propagations - Total test time > 10 minutes.)" OFF)


option(TUDAT_BUILD_WITH_ESTIMATION_TOOLS "Build tudat with estimation tools." OFF)

# Build option: enable tutorials.
option(TUDAT_BUILD_TUDAT_TUTORIALS "Build tutorials." ON)

# Build option: enable tutorials.
option(TUDAT_BUILD_PAGMO_TUTORIALS "Build pagmo2 tutorials." OFF)

# Build option: enable tutorials.
option(TUDAT_BUILD_TBB_SUPPORT "Build tbb support for Pagmo." OFF)

# Build static library instead of dynamic.
option(TUDAT_BUILD_STATIC_LIBRARY "Build tudat as a static library, instead of dynamic." OFF)

# Build static library instead of dynamic.
option(TUDAT_BUILD_WITH_FILTERS "Build Tudat with preliminary Kalman Filter functionality." OFF)

# Build sofa interface.
option(TUDAT_BUILD_WITH_SOFA_INTERFACE "Build Tudat with sofa interface." ON)

# Build json interface.
option(TUDAT_BUILD_WITH_JSON_INTERFACE "Build Tudat with json interface." OFF)

# Build spice interface.
option(TUDAT_BUILD_WITH_SPICE_INTERFACE "Build Tudat with spice interface." ON)

# Build extended precision propagation tools.
option(TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS "Build tudat with extended precision propagation tools." OFF)

#set(CMAKE_INSTALL_PREFIX "/home/ggarrett/Repositories/new/tudatBundle/tudat/install_test")

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include CACHE PATH "Installation directory for header files")
set(INSTALL_DATA_DIR ${CMAKE_INSTALL_PREFIX}/data/ CACHE PATH "Installation directory for data files")
set(INSTALL_MAN_DIR ${CMAKE_INSTALL_PREFIX} CACHE PATH "Installation directory for man documentation")
set(INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR}/cmake/tudat CACHE PATH "Installation directory for cmake config files")

# Make relative paths absolute (needed later on)
foreach (p LIB BIN INCLUDE DATA CMAKE)
    set(var INSTALL_${p}_DIR)
    set(RELATIVE_INSTALL_${p}_DIR ${INSTALL_${p}_DIR})
    if (NOT IS_ABSOLUTE "${${var}}")
        set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif ()
endforeach ()

# Find relative paths for cmake-config file. (for find_package)
set(Tudat_CMAKE_DIR "${INSTALL_CMAKE_DIR}")
file(RELATIVE_PATH rel_include_dir "${Tudat_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
file(RELATIVE_PATH rel_lib_dir "${Tudat_CMAKE_DIR}" "${INSTALL_LIB_DIR}")
file(RELATIVE_PATH rel_data_dir "${Tudat_CMAKE_DIR}" "${INSTALL_DATA_DIR}")
file(RELATIVE_PATH rel_bin_dir "${Tudat_CMAKE_DIR}" "${INSTALL_BIN_DIR}")
file(RELATIVE_PATH rel_man_dir "${Tudat_CMAKE_DIR}" "${INSTALL_MAN_DIR}")

# Set relative paths for config.cmake.
foreach (p include lib data bin man)
    string(TOUPPER ${p} P)
    set(RELATIVE_INSTALL_${P}_DIR ${rel_${p}_dir})
endforeach ()

# Message (should be for debug mode).
message(STATUS RELATIVE_INSTALL_LIB_DIR: ${RELATIVE_INSTALL_LIB_DIR})
message(STATUS RELATIVE_INSTALL_BIN_DIR: ${RELATIVE_INSTALL_BIN_DIR})
message(STATUS RELATIVE_INSTALL_INCLUDE_DIR: ${RELATIVE_INSTALL_INCLUDE_DIR})
message(STATUS RELATIVE_INSTALL_DATA_DIR: ${RELATIVE_INSTALL_DATA_DIR})
message(STATUS RELATIVE_INSTALL_MAN_DIR: ${RELATIVE_INSTALL_MAN_DIR})

# Enforce TBB if pagmo tutorials are to be built.
if (TUDAT_BUILD_PAGMO_TUTORIALS)
    if (TUDAT_BUILD_TBB_SUPPORT)

    else ()
        message("Forcing TBB support for Pagmo tutorial build.")
        set(TUDAT_BUILD_TBB_SUPPORT ON)
    endif ()
endif ()

# Set cmake build-type. If it not supplied by the user (either directly as an argument of through
# the "UserSettings.txt" file, the default built type is "Release".
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_BUILD_TYPE Release)
elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE Debug)
endif ()

message(STATUS "<< ${PROJECT_NAME} (${CMAKE_BUILD_TYPE} - ${BUILD_STYLE}) >>")


# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(compiler)

# NOTE: make sure the include directories from the current build
# are included first, so that if there is already a pagmo installation
# in the prefix path we don't risk including the headers from that
# one instead.
#target_include_directories(tudat PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
#        $<INSTALL_INTERFACE:include>)

# Find Eigen3 library on local system.
find_package(Eigen3 REQUIRED)

# Include Eigen3 directories.
# Set cmake flag to suppress Eigen warnings (platform-dependent solution).
if (NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${EIGEN3_INCLUDE_DIR}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${EIGEN3_INCLUDE_DIR}\"")
endif ()

# Adds Boost_UNIT_TEST_FRAMEWORK_LIBRARY if required.
if (TUDAT_BUILD_TESTS)
    set(_TUDAT_FIND_BOOST_UNIT_TEST_FRAMEWORK ON)
endif ()
#
# Boost
#
# Configure Boost libraries.
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
include(TudatFindBoost)

# Include Boost directories.
# Set CMake flag to suppress Boost warnings (platform-dependent solution).
if (NOT APPLE OR APPLE_INCLUDE_FORCE)
    include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif ()

#
# DOCUMENTATION.
#
# Add an option to toggle the generation of the API documentation.
# If documentation should be built, find Doxygen package and setup config file.
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if (BUILD_DOCUMENTATION)
    find_package(Doxygen)

    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR
                "Doxygen is needed to build the documentation. Please install it correctly.")
    endif ()

    # Configure the Template Doxyfile for our specific project.
    configure_file(${CODEROOT}/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)

    # Add a custom target to run Doxygen when ever the project is built.
    # IF you do NOT want the documentation to be generated EVERY time you build the project
    # then leave out the 'ALL' keyword from the above command.
    add_custom_target(docs ALL COMMAND ${DOXYGEN_EXECUTABLE}
            ${PROJECT_BINARY_DIR}/Doxyfile
            SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)

    # Install documentation.
    install(DIRECTORY "${CODEROOT}/Documentation"
            DESTINATION "${TUDAT_DISTRIBUTION_PATH}")
endif ()

#
# CSPICE.
#
if (TUDAT_BUILD_WITH_SPICE_INTERFACE)
    find_package(CSpice REQUIRED 1.0.0)
    add_definitions(-DTUDAT_BUILD_WITH_SPICE_INTERFACE=1)
    # Include CSpice directories.
    # Set CMake flag to suppress CSpice warnings (platform-dependent solution).
    if (NOT APPLE OR APPLE_INCLUDE_FORCE)
        include_directories(SYSTEM AFTER "${CSpice_INCLUDE_DIRS}")
    else ()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${CSpice_INCLUDE_DIRS}\"")
    endif ()
else ()
    add_definitions(-DTUDAT_BUILD_WITH_SPICE_INTERFACE=1)
endif ()

# JSON
if (TUDAT_BUILD_WITH_JSON_INTERFACE)
    find_package(nlohmann_json REQUIRED 3.7.3)
endif ()

# SOFA
if (TUDAT_BUILD_WITH_SOFA_INTERFACE)
    find_package(Sofa REQUIRED 1.0.0)
    add_definitions(-DTUDAT_BUILD_WITH_SOFA_INTERFACE=1)
    # Include Sofa directories.
    # Set CMake flag to suppress CSpice warnings (platform-dependent on).
    if (NOT APPLE OR APPLE_INCLUDE_FORCE)
        include_directories(SYSTEM AFTER "${Sofa_INCLUDE_DIRS}")
    else ()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Sofa_INCLUDE_DIRS}\"")
    endif ()
else ()
    add_definitions(-DTUDAT_BUILD_WITH_SPICE_INTERFACE=1)
endif ()


#
# NRLMSISE-00
#
if (NOT BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS)
    message(STATUS "Extended precision propagation disabled!")
    add_definitions(-DBUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=0)
else ()
    message(STATUS "Extended precision propagation enabled!")
    add_definitions(-DBUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=1)
endif ()


# Set compiler based on preferences (e.g. USE_CLANG) and system.
include(tudatLinkLibraries)


if (USE_TUDAT_EXAMPLE_APPLICATIONS)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/examples")
endif ()

function(TUDAT_ADD_TEST_DATA arg1)
    # We move the test files into the testing directory
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${arg1}/" DESTINATION "${CMAKE_BINARY_DIR}/${arg1}/")
endfunction()

function(TUDAT_ADD_TESTCASE arg1)
    if (TUDAT_BUILD_TESTS)
        # arg1: test name
        # arg2: test path
        #        get_filename_component(name_top "${CMAKE_CURRENT_SOURCE_DIR}/.." NAME)
        get_filename_component(name_bot ${CMAKE_CURRENT_SOURCE_DIR} NAME)
        set(executable_name "test_${name_bot}_${arg1}")
        add_executable(${executable_name} ${CMAKE_CURRENT_SOURCE_DIR}/tests/unitTest${arg1}.cpp)
        if (CMAKE_BUILD_TYPE MATCHES DEBUG)
            if (${arg1} MATCHES "tudat_sofa_interface")
                message("Tudat TESTCASE: ${arg1}")
                message("ARGN: ${ARGN}")
            endif ()
        endif ()
        set_target_properties(${executable_name}
                PROPERTIES
                LINKER_LANGUAGE CXX
                ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
                )
        # All tests require Boost_LIBRARIES, so they are kept constant here.
        target_link_libraries(${executable_name} PUBLIC
                "${ARGN}"
                "${Boost_LIBRARIES}"
                )
        # Has to do with a fundamental change in the way that vectors, maps
        # etc work compared to the unit-test-framework. Not the best link for
        # resource but:
        # https://stackoverflow.com/questions/33644088/linker-error-while-building-unit-tests-with-boost
        target_compile_definitions(${executable_name} PUBLIC
                D_GLIBCXX_USE_CXX11_ABI=0
                )
        #        target_compile_options(${executable_name} PRIVATE
        #                "$<$<CONFIG:Debug>:${TUDAT_CXX_FLAGS_DEBUG}>"
        #                "$<$<CONFIG:Release>:${TUDAT_CXX_FLAGS_RELEASE}>"
        #                "$<$<CONFIG:RelWithDebInfo>:${TUDAT_CXX_FLAGS_RELEASE}>"
        #                "$<$<CONFIG:MinSizeRel>:${TUDAT_CXX_FLAGS_RELEASE}>"
        #                )
        # Let's setup the target C++ standard, but only if the user did not provide it manually.
        if (NOT CMAKE_CXX_STANDARD)
            set_property(TARGET ${executable_name} PROPERTY CXX_STANDARD 17)
        endif ()
        set_property(TARGET ${executable_name} PROPERTY CXX_STANDARD_REQUIRED YES)
        set_property(TARGET ${executable_name} PROPERTY CXX_EXTENSIONS NO)
        add_test(${executable_name} "${CMAKE_BINARY_DIR}/tests/${executable_name}")
    endif ()
endfunction()

function(TUDAT_ADD_LIBRARY arg1 arg2 arg3)
    # arg1: library name
    # arg3: sources
    # arg4: headers
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        if ("${arg1}" STREQUAL "sofa_interface")
            message("Added library: " tudat_${arg1})
            message("Tudat LIBRARY: ${arg1}")
            message("Tudat LIBRARY sources: ${arg2}")
            message("Tudat LIBRARY headers: ${arg3}")
            message("ARGN: ${ARGN}")
        endif ()
    endif ()
    if (TUDAT_BUILD_STATIC_LIBRARY)
        add_library("tudat_${arg1}" STATIC ${arg2})
    else ()
        add_library("tudat_${arg1}" SHARED ${arg2})
    endif ()

    # Link libraries if extra arguments provided.
    target_link_libraries("tudat_${arg1}" PUBLIC "${ARGN}")

    # BUILD-TREE OUTPUT.
    set_target_properties(tudat_${arg1}
            PROPERTIES
            LINKER_LANGUAGE CXX
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            )

    # INSTALL TREE OUTPUT.
    install(TARGETS tudat_${arg1}
            EXPORT tudat_targets
            LIBRARY DESTINATION ${RELATIVE_INSTALL_LIB_DIR}
            ARCHIVE DESTINATION ${RELATIVE_INSTALL_LIB_DIR}
            RUNTIME DESTINATION ${RELATIVE_INSTALL_BIN_DIR}
            #            PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
            INCLUDES DESTINATION "${RELATIVE_INSTALL_INCLUDE_DIR}"
            COMPONENT Development)
endfunction()


function(TUDAT_ADD_EXECUTABLE arg1 arg2)
    # arg1: executable name
    # arg2: sources
    if (CMAKE_BUILD_TYPE MATCHES DEBUG)
        message("Tudat EXECUTABLE: ${arg1}")
        message("ARGN: ${ARGN}")
    endif ()
    add_executable(${arg1} ${arg2})
    target_link_libraries(${arg1}
            ${ARGN}
            ${Boost_LIBRARIES}
            )
    set("tudat_${arg1}_INCLUDES"
            ${CMAKE_CURRENT_SOURCE_DIR}
            )
    target_include_directories(${arg1}
            PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}"
            )

    # BUILD-TREE OUTPUT.
    set_target_properties(${arg1}
            PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_INCLUDE_DIR}"
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            )

    # INSTALL TREE OUTPUT.
    install(TARGETS ${arg1}
            EXPORT tudat_targets
            LIBRARY DESTINATION ${RELATIVE_INSTALL_LIB_DIR}
            ARCHIVE DESTINATION ${RELATIVE_INSTALL_LIB_DIR}
            RUNTIME DESTINATION ${RELATIVE_INSTALL_BIN_DIR}
            #            PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}"
            INCLUDES DESTINATION "${RELATIVE_INSTALL_INCLUDE_DIR}"
            COMPONENT Development)

endfunction()

function(TUDAT_ADD_DATA_DIR arg1)
    get_filename_component(dir_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${arg1} DESTINATION ${INSTALL_DATA_DIR}/${dir_name})
endfunction()

# Include libraries of Tudat.
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

include_directories(include/)
add_subdirectory(src/astro)
add_subdirectory(src/basics)
add_subdirectory(src/interface)
add_subdirectory(src/math)
add_subdirectory(src/simulation)
add_subdirectory(src/io)
add_subdirectory(src/external)
add_subdirectory(examples)

# The interesting stuff goes here
# ===============================

get_property(tudat_library_list GLOBAL PROPERTY tudat_library_list_PROPERTY)
#set(Tudat_LIBRARIES ${tudat_library_alias_list})

# Add all targets to the build-tree export set
export(TARGETS ${Tudat_LIBRARIES}
        FILE "${PROJECT_BINARY_DIR}/tudat_targets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE Tudat)

## Create the tudat-config.cmake and tudat-config-version.cmake files
#file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
#        "${INSTALL_INCLUDE_DIR}")
## ... for the build tree
#set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
#configure_file(tudat-config.cmake.in
#        "${PROJECT_BINARY_DIR}/tudat-config.cmake" @ONLY)
## ... for the install tree
#set(CONF_INCLUDE_DIRS "\${Tudat_CMAKE_DIR}/${rel_include_dir}")
#configure_file(tudat-config.cmake.in
#        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/tudat-config.cmake" @ONLY)

# Configure cmake files.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tudat-config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/tudat-config.cmake" @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tudat-config.cmake" DESTINATION ${INSTALL_CMAKE_DIR})

# Install the export set for use with the install-tree
install(EXPORT tudat_targets DESTINATION "${INSTALL_CMAKE_DIR}")

# Take care of versioning.
include(CMakePackageConfigHelpers)
# NOTE: SameMinorVersion available only
# since CMake 3.11.
if (${CMAKE_VERSION} VERSION_LESS "3.11.0")
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/tudat-config-version.cmake" VERSION ${Tudat_VERSION}
            COMPATIBILITY SameMajorVersion)
else ()
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/tudat-config-version.cmake" VERSION ${Tudat_VERSION}
            COMPATIBILITY SameMinorVersion)
endif ()
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tudat-config-version.cmake" DESTINATION "${RELATIVE_INSTALL_CMAKE_DIR}")

install(DIRECTORY include/ DESTINATION "${RELATIVE_INSTALL_INCLUDE_DIR}")


#message(chacha ${tudat_library_alias_list})
#add_executable(testing_test test.cpp
#        )
#target_link_libraries(testing_test
#        ${Tudat_PROPAGATION_LIBRARIES}
#        tudat_electromagnetism
#        tudat_spice_interface
#        ${CSpice_LIBRARIES}
#        ${Boost_LIBRARIES}
#        ${Sofa_LIBRARIES})

#get_property(Tudat_LIBRARY_LIST GLOBAL PROPERTY tudat_library_list_PROPERTY)

#target_link_libraries(tudat INTERFACE "${Tudat_LIBRARY_LIST}")

#add_library(tudat_libraries INTERFACE)
#target_compile_definitions(tudat_libraries INTERFACE ${Tudat_DEFINITIONS})
#target_include_directories(tudat_libraries INTERFACE
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/tudat>
#        $<INSTALL_INTERFACE:include/tudat>  # <prefix>/include/mylib
#        )
#target_link_libraries(tudat_libraries INTERFACE ${tudat_library_alias_list})
#export(TARGETS tudat_libraries
#        APPEND FILE ${INSTALL_CMAKE_DIR}/tudat_targets_export.cmake)
#
#if (TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
#    add_library(tudat_estimation INTERFACE)
#    target_compile_definitions(tudat_estimation INTERFACE ${Tudat_DEFINITIONS})
#    target_include_directories(tudat_estimation INTERFACE
#            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/tudat>
#            $<INSTALL_INTERFACE:include/tudat>  # <prefix>/include/mylib
#            )
#    target_link_libraries(tudat_estimation INTERFACE ${Tudat_ESTIMATION_LIBRARIES})
#    export(TARGETS tudat_estimation
#            APPEND FILE ${INSTALL_CMAKE_DIR}/tudat_estimation_export.cmake)
#endif ()

#install(EXPORT tudat_targets_export
#        NAMESPACE Tudat::
#        DESTINATION ${RELATIVE_INSTALL_CMAKE_DIR}
#        COMPONENT Development)

# Installation of the header files.


# Installation of the data files.
#install(DIRECTORY data/ DESTINATION "${RELATIVE_INSTALL_DATA_DIR}")

