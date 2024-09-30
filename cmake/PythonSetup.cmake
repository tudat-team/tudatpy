if(YACMAPythonSetupIncluded)
    return()
endif()

if (PREFIX)
    message(STATUS "Using PREFIX: ${PREFIX}")
elseif(CONDA_PREFIX)
    message(STATUS "Using CONDA_PREFIX: ${CONDA_PREFIX}")
    set(PREFIX ${CONDA_PREFIX})
elseif(CMAKE_PREFIX_PATH)
    message(STATUS "Using CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
    set(PREFIX ${CMAKE_PREFIX_PATH})
else()
    message(FATAL_ERROR "We need a prefix")
endif()

set(Python_ROOT ${PREFIX})
find_package(Python REQUIRED COMPONENTS Interpreter Development)

# Python include directory
if (NOT Python_INCLUDE_DIRS)
    message(FATAL_ERROR "Could not determine the Python include dir.")
else()
    set(
        YACMA_PYTHON_INCLUDE_DIR "${Python_INCLUDE_DIRS}"
        CACHE PATH "Path to the Python include dir."
    )
endif()
add_library(YACMA::PythonIncludeDir INTERFACE IMPORTED)
set_target_properties(
    YACMA::PythonIncludeDir PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${YACMA_PYTHON_INCLUDE_DIR}
)

# Python executable
if (NOT Python_EXECUTABLE)
    message(FATAL_ERROR "Could not determine the Python executable.")
else()
    set(PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
endif()

# Python version
if (NOT Python_VERSION)
    message(FATAL_ERROR "Could not determine the Python version.")
else()
    set(PYTHON_VERSION_STRING "${Python_VERSION}")
endif()

# Python libraries
if (WIN32)
    if (NOT Python_LIBRARIES)
        message(FATAL_ERROR "Could not determine the Python libraries.")
    else()
        set(PYTHON_LIBRARIES "${Python_LIBRARIES}")
    endif()
endif()


# Platform-specific installation directory
if (NOT Python_SITEARCH)
    message(FATAL_ERROR "Could not determine the Python installation directory.")
else()
    string(REPLACE "/" ";" _install_dir ${Python_SITEARCH})
    list(FIND _install_dir "lib" _lib_index)
    list(SUBLIST _install_dir ${_lib_index} -1 _install_dir)
    list(JOIN _install_dir "/" YACMA_PYTHON_MODULES_INSTALL_PATH)
endif()

# Display info
message(STATUS "Python interpreter: ${PYTHON_EXECUTABLE}")
message(STATUS "Python interpreter version: ${PYTHON_VERSION_STRING}")
if (WIN32)
    message(STATUS "Python libraries: ${PYTHON_LIBRARIES}")
endif()
message(STATUS "Python include dir: ${YACMA_PYTHON_INCLUDE_DIR}")
message(STATUS "Python modules install path: ${YACMA_PYTHON_MODULES_INSTALL_PATH}")


# Module extension
if (UNIX)
    if (APPLE)
        message(STATUS "OS X platform detected: Extension is .so")
        set(_YACMA_PY_MODULE_EXTENSION ".so")
    else()
        message(STATUS "Generic UNIX platform detected")
    endif()

    if(NOT YACMA_PYTHON_MODULES_INSTALL_PATH)
        message(WARNING "THIS SHOULD NOT BE CALLED")
        execute_process(
            COMMAND ${PYTHON_EXECUTABLE} -c
            "from __future__ import print_function\nimport distutils.sysconfig\nimport os\nprint(os.path.split(distutils.sysconfig.get_python_lib())[-1])"
            OUTPUT_VARIABLE _YACMA_PY_PACKAGES_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        message(STATUS "Python packages dir is: ${_YACMA_PY_PACKAGES_DIR}")
        set(
            YACMA_PYTHON_MODULES_INSTALL_PATH
            "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/${_YACMA_PY_PACKAGES_DIR}"
            CACHE PATH "Install path for Python modules."
        )
        mark_as_advanced(YACMA_PYTHON_MODULES_INSTALL_PATH)
    endif()
    elseif(WIN32)
        message(STATUS "Windows platform detected.")
        message(STATUS "Output extension for compiled modules will be '.pyd'.")
        set(_YACMA_PY_MODULE_EXTENSION "pyd")
        if(NOT YACMA_PYTHON_MODULES_INSTALL_PATH)
            message(WARNING "THIS SHOULD NOT BE CALLED!")
            # On Windows, we will install directly into the install path of the Python interpreter.
            execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
                    OUTPUT_VARIABLE _YACMA_PYTHON_MODULES_INSTALL_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
            set(YACMA_PYTHON_MODULES_INSTALL_PATH "${_YACMA_PYTHON_MODULES_INSTALL_PATH}" CACHE PATH "Install path for Python modules.")
            mark_as_advanced(YACMA_PYTHON_MODULES_INSTALL_PATH)
        endif()
    else()
        message(FATAL_ERROR "Platform not supported.")
endif()

# Check the install path was actually detected.
if("${YACMA_PYTHON_MODULES_INSTALL_PATH}" STREQUAL "")
    message(FATAL_ERROR "Python module install path not detected correctly.")
endif()

function(YACMA_PYTHON_MODULE name)
    message(STATUS "Setting up the compilation of the Python module '${name}'.")
    # If we need an explicit link to the Python library, we compile it as a normal shared library.
    # Otherwise, we compile it as a module.
    if(_YACMA_PYTHON_MODULE_NEED_LINK)
        add_library("${name}" SHARED ${ARGN})
    else()
        add_library("${name}" MODULE ${ARGN})
    endif()
    # Any "lib" prefix normally added by CMake must be removed.
    set_target_properties("${name}" PROPERTIES PREFIX "")
    if(NOT ${_YACMA_PY_MODULE_EXTENSION} STREQUAL "")
        # If needed, set a custom extension for the module.
        message(STATUS "Setting up custom extension '${_YACMA_PY_MODULE_EXTENSION}' for the Python module '${name}'.")
        set_target_properties("${name}" PROPERTIES SUFFIX ".${_YACMA_PY_MODULE_EXTENSION}")
    endif()
    # We need extra flags to be set when compiling Python modules, at least
    # with clang and gcc. See:
    # https://bugs.python.org/issue11149
    # http://www.python.org/dev/peps/pep-3123/
    # NOTE: do not use the yacma compiler linker settings bits, so this module
    # can be used stand-alone.
    if(CMAKE_COMPILER_IS_GNUCXX OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" AND NOT MSVC))
        message(STATUS "Setting up extra compiler flag '-fwrapv' for the Python module '${name}'.")
        target_compile_options(${name} PRIVATE "-fwrapv")
    endif()
    if(APPLE AND ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        # On OSX + Clang this link flag is apparently necessary in order to avoid
        # undefined references to symbols defined in the Python library. See also:
        # https://github.com/potassco/clingo/issues/79
        # https://stackoverflow.com/questions/25421479/clang-and-undefined-symbols-when-building-a-library
        # https://cmake.org/pipermail/cmake/2017-March/065115.html
        set_target_properties(${name} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif()

    # Add the Python include dirs.
    target_include_directories("${name}" SYSTEM PRIVATE ${YACMA_PYTHON_INCLUDE_DIR})

    # Link to the Python libs, if necessary.
    if(_YACMA_PYTHON_MODULE_NEED_LINK)
        target_link_libraries("${name}" PRIVATE ${PYTHON_LIBRARIES})
    endif()
endfunction()

# Mark as included.
set(YACMAPythonSetupIncluded YES)
