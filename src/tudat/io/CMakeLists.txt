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

if (TUDAT_BUILD_WITH_SERIALIZATION)
    add_subdirectory(serialization)
endif ()

# Neutral environment readers compile once without links to dynamics models.
set(environment_data_SOURCES
        "readEarthOrientationData.cpp"
        "readOceanLoadingData.cpp"
        "readOceanPoleTideData.cpp"
        "readAtmosphericTideData.cpp"
        "readTerrestrialFrameTransformationData.cpp"
        "readSinexStationHistory.cpp"
        "readStationLoadingData.cpp"
        "readIcgemGravityField.cpp"
        "readOceanPoleGravityData.cpp"
        "readHarmonicPotentialData.cpp"
        "readAod1bGravityField.cpp")

# Installed interfaces for the observations and their file readers.
set(environment_data_HEADERS
        "textParsing.h"
        "readEarthOrientationData.h"
        "readOceanLoadingData.h"
        "readOceanPoleTideData.h"
        "readAtmosphericTideData.h"
        "readTerrestrialFrameTransformationData.h"
        "readSinexStationHistory.h"
        "readStationLoadingData.h"
        "readIcgemGravityField.h"
        "readOceanPoleGravityData.h"
        "readHarmonicPotentialData.h"
        "readAod1bGravityField.h")

TUDAT_ADD_LIBRARY("environment_data"
        "${environment_data_SOURCES}"
        "${environment_data_HEADERS}")

# Preserve language, shared-library and public Eigen requirements for all consumers.
target_compile_features(tudat_environment_data PUBLIC cxx_std_17)
set_target_properties(tudat_environment_data PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(TARGET Eigen3::Eigen)
    target_link_libraries(tudat_environment_data PUBLIC Eigen3::Eigen)
else()
    target_include_directories(tudat_environment_data SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIRS})
endif()

# Only the compatibility adapter needs the old aggregate's other dependencies.
target_link_libraries(tudat_input_output_raw PUBLIC tudat_environment_data)
