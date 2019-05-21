cmake_minimum_required(VERSION 2.8.3)


#CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(TudatPy)

#message("CHAAAA ${CMAKE_SOURCE_DIR}")
include_directories ("${CMAKE_SOURCE_DIR}/tudat")
include_directories ("${CMAKE_SOURCE_DIR}/eigen")
include_directories ("${CMAKE_SOURCE_DIR}/cspice")
find_package(Spice)
find_package(Tudat 2.0 REQUIRED)
include(tudatLinkLibraries)

include_directories(SYSTEM AFTER "${SPICE_INCLUDE_DIR}")
#
# TODO: Allow for choice between python 3.6 and python 2.7 for compilation.
FIND_PACKAGE(PythonInterp 3.6 REQUIRED)

# TODO: Boost component for python must only contain major. Modify and push
# suggestion on github.

if (PYTHONINTERP_FOUND)
    if (UNIX AND NOT APPLE)
        if (PYTHON_VERSION_MAJOR EQUAL 3)
            # NOTE: As far as is known, building boost.python3 provides the lib files
            # in the form of python3 and not python3.6 as is suggested by original
            # code.
            FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR})
            FIND_PACKAGE(PythonInterp 3)
            FIND_PACKAGE(PythonLibs 3 REQUIRED)
        else()
            FIND_PACKAGE(Boost COMPONENTS python)
            FIND_PACKAGE(PythonInterp)
            FIND_PACKAGE(PythonLibs REQUIRED)
        endif()
    else()
        if (PYTHON_VERSION_MAJOR EQUAL 3)
            FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
            FIND_PACKAGE(PythonInterp 3)
            FIND_PACKAGE(PythonLibs 3 REQUIRED)
        else()
            FIND_PACKAGE(Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
            FIND_PACKAGE(PythonInterp)
            FIND_PACKAGE(PythonLibs REQUIRED)
        endif()
    endif()
else()
    message("Python not found")
endif()

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")

ENABLE_TESTING()
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
LINK_LIBRARIES(${Boost_LIBRARIES} ${PYTHON_LIBRARIES}) # Deprecated but so convenient!

ADD_SUBDIRECTORY(src/constants)
ADD_SUBDIRECTORY(src/simulation)
