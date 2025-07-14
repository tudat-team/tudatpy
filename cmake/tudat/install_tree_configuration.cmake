set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set(INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(INSTALL_DATA_DIR "${CMAKE_INSTALL_PREFIX}/data")
set(INSTALL_MAN_DIR "${CMAKE_INSTALL_PREFIX}")
#set(INSTALL_TESTS_DIR "${CMAKE_INSTALL_PREFIX}/tests")
set(INSTALL_CMAKE_DIR "${INSTALL_LIB_DIR}/cmake/${PROJECT_NAME_LOWER}")

message("")
message(STATUS "INSTALLATION PREFIX: ${CMAKE_INSTALL_PREFIX}")
message("")


# Make relative paths absolute (needed later on)
foreach (p LIB BIN INCLUDE DATA CMAKE)
    set(var INSTALL_${p}_DIR)
    set(RELATIVE_INSTALL_${p}_DIR ${INSTALL_${p}_DIR})
    if (NOT IS_ABSOLUTE "${${var}}")
        set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif ()
endforeach ()

# Find relative paths for cmake-config file. (for find_package)
file(RELATIVE_PATH rel_include_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
file(RELATIVE_PATH rel_lib_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_LIB_DIR}")
file(RELATIVE_PATH rel_data_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_DATA_DIR}")
file(RELATIVE_PATH rel_bin_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_BIN_DIR}")
file(RELATIVE_PATH rel_man_dir "${INSTALL_CMAKE_DIR}" "${INSTALL_MAN_DIR}")

# Set relative paths for config.cmake.
foreach (p include lib data bin man)
    string(TOUPPER ${p} P)
    set(RELATIVE_INSTALL_${P}_DIR ${rel_${p}_dir})
endforeach ()
