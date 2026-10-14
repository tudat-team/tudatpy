#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add source files.
set(io_SOURCES
        "basicInputOutput.cpp"
        "dictionaryComparer.cpp"
        "dictionaryTools.cpp"
        "fieldValue.cpp"
        "fixedWidthParser.cpp"
        "linearFieldTransform.cpp"
        "missileDatcomData.cpp"
        "missileDatcomReader.cpp"
        "parsedDataVectorUtilities.cpp"
        "separatedParser.cpp"
        "textParser.cpp"
        "twoLineElementData.cpp"
        "twoLineElementsTextFileReader.cpp"
        "streamFilters.cpp"
        "parseSolarActivityData.cpp"
        "extractSolarActivityData.cpp"
        "solarActivityData.cpp"
        "preProcessFdetsFile.cpp"
        "preProcessIfmsFile.cpp"
        "multiDimensionalArrayReader.cpp"
        "aerodynamicCoefficientReader.cpp"
        "tabulatedAtmosphereReader.cpp"
        "util.cpp"
        "readOdfFile.cpp"
        "preProcessOdfFile.cpp"
        "readPsfFile.cpp"
        "readCrdFile.cpp"
        "readSinexFile.cpp"
        "mapTextFileReader.cpp"
        "readTabulatedMediaCorrections.cpp"
        "readTabulatedWeatherData.cpp"
        "readTrackingTxtFile.cpp"
        "readVariousPdsFiles.cpp"
        "readViennaMappingFunctionData.cpp"
        "readIonexFile.cpp"
        "readInpopEphemerisFile.cpp"
        "readNeQuick2Data.cpp"
        "comaModelInputOutput.cpp"
)

# Add header files.
set(io_HEADERS
        "basicInputOutput.h"
        "dictionaryComparer.h"
        "dictionaryEntry.h"
        "dictionaryTools.h"
        "extractor.h"
        "fieldTransform.h"
        "fieldType.h"
        "fieldValue.h"
        "fixedWidthParser.h"
        "linearFieldTransform.h"
        "missileDatcomData.h"
        "missileDatcomReader.h"
        "parsedDataVectorUtilities.h"
        "parser.h"
        "separatedParser.h"
        "textParser.h"
        "twoLineElementData.h"
        "twoLineElementsTextFileReader.h"
        "basicInputOutput.h"
        "mapTextFileReader.h"
        "matrixTextFileReader.h"
        "streamFilters.h"
        "parseSolarActivityData.h"
        "extractSolarActivityData.h"
        "solarActivityData.h"
        "multiDimensionalArrayReader.h"
        "multiDimensionalArrayWriter.h"
        "aerodynamicCoefficientReader.h"
        "readHistoryFromFile.h"
        "tabulatedAtmosphereReader.h"
        "util.h"
        "preProcessFdetsFile.h"
        "preProcessIfmsFile.h"
        "readOdfFile.h"
        "preProcessOdfFile.h"
        "readPsfFile.h"
        "readCrdFile.h"
        "readSinexFile.h"
        "readBinaryFile.h"
        "readTabulatedMediaCorrections.h"
        "readTabulatedWeatherData.h"
        "readTrackingTxtFile.h"
        "readVariousPdsFiles.h"
        "readInpopEphemerisFile.h"
        "comaModelInputOutput.h"
        "trackingData.h"
        )

# Raw readers compile once in a target that has no environment or propagation models.
# The legacy aggregate retains only the readers that construct or consume model objects.
set(io_model_SOURCES
        "readInpopEphemerisFile.cpp"
        "readTabulatedMediaCorrections.cpp"
        "readTabulatedWeatherData.cpp")
set(io_model_HEADERS
        "readInpopEphemerisFile.h"
        "readTabulatedMediaCorrections.h"
        "readTabulatedWeatherData.h")
list(REMOVE_ITEM io_SOURCES ${io_model_SOURCES})
list(REMOVE_ITEM io_HEADERS ${io_model_HEADERS})

TUDAT_ADD_LIBRARY("input_output_raw"
        "${io_SOURCES}"
        "${io_HEADERS}"
        PUBLIC_LINKS
        tudat_interpolators
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_basics
        PRIVATE_LINKS
        tudat_sofa_interface)

TUDAT_ADD_LIBRARY("input_output"
        "${io_model_SOURCES}"
        "${io_model_HEADERS}"
        PUBLIC_LINKS
        tudat_input_output_raw
        tudat_ephemerides)

# Check actual target links, including transitive and private static-library links.
include("${PROJECT_SOURCE_DIR}/cmake/tudat/CheckRawInputOutputDependencies.cmake")
cmake_language(DEFER DIRECTORY "${PROJECT_SOURCE_DIR}" CALL tudat_check_raw_input_output_dependencies)

if (TUDAT_BUILD_WITH_SERIALIZATION)
    add_subdirectory(serialization)
endif ()

# Pure EOP file loading is separately linkable, without dynamics/estimation,
# interpolation, ephemerides, SOFA, or serialization dependencies.
if(TUDAT_BUILD_STATIC_LIBRARY)
    add_library(tudat_earth_orientation_data STATIC readEarthOrientationData.cpp)
else()
    add_library(tudat_earth_orientation_data SHARED readEarthOrientationData.cpp)
endif()
target_compile_features(tudat_earth_orientation_data PUBLIC cxx_std_17)
target_include_directories(tudat_earth_orientation_data PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
set_target_properties(tudat_earth_orientation_data PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
add_library(Tudat::tudat_earth_orientation_data ALIAS tudat_earth_orientation_data)
install(TARGETS tudat_earth_orientation_data EXPORT tudat_export
    LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
    ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
    RUNTIME DESTINATION "${INSTALL_BIN_DIR}")

# BLQ input records/parsing are independent of dynamics, estimation and the legacy I/O library.
if(TUDAT_BUILD_STATIC_LIBRARY)
    add_library(tudat_ocean_loading_data STATIC readOceanLoadingData.cpp)
else()
    add_library(tudat_ocean_loading_data SHARED readOceanLoadingData.cpp)
endif()
target_include_directories(tudat_ocean_loading_data PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
if(TARGET Eigen3::Eigen)
    target_link_libraries(tudat_ocean_loading_data PUBLIC Eigen3::Eigen)
else()
    target_include_directories(tudat_ocean_loading_data SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIRS})
endif()
set_target_properties(tudat_ocean_loading_data PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
add_library(Tudat::tudat_ocean_loading_data ALIAS tudat_ocean_loading_data)
install(TARGETS tudat_ocean_loading_data EXPORT tudat_export
    LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
    ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
    RUNTIME DESTINATION "${INSTALL_BIN_DIR}")
