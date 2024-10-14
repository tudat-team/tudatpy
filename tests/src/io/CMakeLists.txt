#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(MapTextFileReader
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(MatrixTextFileReader
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(StreamFilters
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(TwoLineElementsTextFileReader
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(BasicInputOutput
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(ParsedDataVectorUtilities
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(FieldValue
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(TextParser
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(MissileDatcomReader
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(MissileDatcomData
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(DictionaryInputSystem
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(Extractor
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(FixedWidthParser
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(SeparatedParser
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(LinearFieldTransform
        PRIVATE_LINKS
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(SolarActivityData
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(MultiArrayReader
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(MultiArrayWriter
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(AerodynamicCoefficientReader
        PRIVATE_LINKS
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basics
        )

TUDAT_ADD_TEST_CASE(OdfFileReader
        PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES} )

TUDAT_ADD_TEST_CASE(ReadTabulatedMediaCorrections
        PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES} )

TUDAT_ADD_TEST_CASE(ReadTabulatedWeatherData
        PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES} )

TUDAT_ADD_TEST_CASE(TrackingTxtFileReader
        PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES} )

TUDAT_ADD_TEST_CASE(IfmsFileReader
        PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES} )