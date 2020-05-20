
cmake_minimum_required(VERSION 3.16)


# Set library and include path to conda test environment.
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/home/ggarrett/anaconda3/envs/test/lib")
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/home/ggarrett/anaconda3/envs/test/include")

project(TestTudatProject LANGUAGES CXX)

# Find package test.
find_package(Tudat)
include_directories(${Tudat_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})

# MESSAGES.
message(Tudat_PROPAGATION_LIBRARIES: "${Tudat_PROPAGATION_LIBRARIES}")
message(Tudat_INCLUDE_DIRS: "${Tudat_INCLUDE_DIRS}")
message(Tudat_LIBRARY_DIRS: "${Tudat_LIBRARY_DIRS}")

# BOOST
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/home/ggarrett/Repositories/new/tudat_bundle/tudat/cmake_modules")
set(_TUDAT_FIND_BOOST_UNIT_TEST_FRAMEWORK ON)
include(TudatFindBoost)
message(Boost_LIBRARIES: "${Boost_LIBRARIES}")


# Tudat function for adding executables.
function(TUDAT_ADD_EXECUTABLE arg1 arg2)
    # arg1: executable name
    # arg2: sources
    if (CMAKE_BUILD_TYPE MATCHES DEBUG)
        message("Tudat EXECUTABLE: ${arg1}")
        message("ARGN: ${ARGN}")
    endif ()
    add_executable(${arg1} ${arg2})
    target_link_libraries(${arg1} ${ARGN} ${Boost_LIBRARIES})
    set_target_properties(${arg1}
            PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_INCLUDE_DIR}"
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            )
endfunction()

TUDAT_ADD_EXECUTABLE(RAAAAAAAAAAAAAAAAAAAAAAAAAAAW
        test.cpp
        ${Tudat_PROPAGATION_LIBRARIES}
        )

# THIS IS REQUIRED. WAS USED
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
