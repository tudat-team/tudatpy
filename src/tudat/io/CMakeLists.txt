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

if (TUDAT_BUILD_WITH_SERIALIZATION)
    add_subdirectory(serialization)
endif ()

# Directory-local registration for the six independently linkable neutral readers.
# Keep this separate from TUDAT_ADD_LIBRARY and its dynamics-oriented dependencies.
function(_tudat_add_neutral_reader target)
    if(TUDAT_BUILD_STATIC_LIBRARY)
        add_library(${target} STATIC ${ARGN})
    else()
        add_library(${target} SHARED ${ARGN})
    endif()
    target_include_directories(${target} PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
    set_target_properties(${target} PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
    add_library(Tudat::${target} ALIAS ${target})
    install(TARGETS ${target} EXPORT tudat_export
        LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
        ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
        RUNTIME DESTINATION "${INSTALL_BIN_DIR}")
endfunction()

_tudat_add_neutral_reader(tudat_earth_orientation_data readEarthOrientationData.cpp readHighFrequencyEopData.cpp)
_tudat_add_neutral_reader(tudat_ocean_loading_data readOceanLoadingData.cpp readOceanPoleTideData.cpp)
_tudat_add_neutral_reader(tudat_atmospheric_loading_data readAtmosphericTideData.cpp)
_tudat_add_neutral_reader(tudat_terrestrial_frame_data readTerrestrialFrameTransformationData.cpp)
_tudat_add_neutral_reader(tudat_station_history_data readSinexStationHistory.cpp)
_tudat_add_neutral_reader(tudat_station_loading_data readStationLoadingData.cpp)

# Preserve the existing feature, runtime-location and Eigen interfaces exactly.
foreach(target IN ITEMS tudat_earth_orientation_data tudat_atmospheric_loading_data
        tudat_terrestrial_frame_data tudat_station_history_data tudat_station_loading_data)
    target_compile_features(${target} PUBLIC cxx_std_17)
endforeach()
foreach(target IN ITEMS tudat_ocean_loading_data tudat_terrestrial_frame_data)
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endforeach()
foreach(target IN ITEMS tudat_ocean_loading_data tudat_atmospheric_loading_data tudat_terrestrial_frame_data)
    if(TARGET Eigen3::Eigen)
        target_link_libraries(${target} PUBLIC Eigen3::Eigen)
    else()
        target_include_directories(${target} SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIRS})
    endif()
endforeach()

# Only the compatibility adapter needs the old aggregate's other dependencies.
target_link_libraries(tudat_input_output_raw PUBLIC tudat_station_history_data)
