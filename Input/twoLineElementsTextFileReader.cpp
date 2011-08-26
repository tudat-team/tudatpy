/*! \file twoLineElementsTextFileReader.cpp
 *    This header file defines a class which can:
 *     - read a 3-line Two-Line Element (TLE) catalog file and store its data,
 *     - check if the data file is valid.
 *    Note that 2-line TLE handling still has to be added.
 *
 *    Path              : /Input/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : J. Leloux
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : j.leloux@tudelft.nl, j.leloux@student.tudelft.nl,
 *                        j.leloux@gmail.com
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 21 February, 2011
 *    Last modified     : 10 August, 2011
 *
 *    References
 *      Leloux, J. Filtering Techniques for Orbital Debris Conjunction Analysis
 *          - applied to SSN TLE catalog data and including astrodynamics and
 *          collision probability theory, MSc Literature Research, Delft
 *          University of Technology, 2010.
 *      Celestrak (a). Space Track TLE Retriever Help,
 *          http://celestrak.com/SpaceTrack/TLERetrieverHelp.asp, 2011. Last
 *          accessed: 5 August, 2011.
 *      Space Track. TLE Format, http://www.space-track.org/tle_format.html,
 *          2004. Last accessed: 5 August, 2011.
 *      Celestrak (b). FAQs: Two-Line Element Set Format,
 *          http://celestrak.com/columns/v04n03/, 2006. Last accessed:
 *          5 August, 2011.
 *      Celestrak (c). NORAD Two-Line Element Set Format,
 *          http://celestrak.com/NORAD/documentation/tle-fmt.asp, 2004. Last
 *          accessed: 5 August, 2011.
 *
 *    Notes
 *      Raw TLE data can be obtained from (Celestrak (a), 2011). Explanations
 *      of the TLE data format can be viewed in (Space Track, 2004),
 *      (Celestrak (b), 2006), and (Celestrak (c), 2004).
 *
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110221    J. Leloux         Startup of TLE header file and class.
 *      110301    J. Leloux         Adjusting header file to parent classes
 *                                  and tudat rules.
 *      110803    J. Leloux         First setup for codecheck.
 *      110805    K. Kumar          Layout and comment corrections; added
 *                                  get-function for vector container of TLE
 *                                  data.
 *      110810    J. Leloux         Tested new setup and changed descriptions.
 */

// Include statements.
#include "twoLineElementsTextFileReader.h"

// Using declarations.
using std::string;
using std::endl;
using std::stringstream;
using std::cout;
using std::cerr;
using std::vector;
using basic_functions::convertStringToTemplate;

//! Default constructor.
TwoLineElementsTextFileReader::TwoLineElementsTextFileReader( )
{
}

//! Default destructor.
TwoLineElementsTextFileReader::~TwoLineElementsTextFileReader( )
{
}

//! Set current year.
void TwoLineElementsTextFileReader::setCurrentYear( const unsigned int& currentYear )
{
    currentYear_ = currentYear;
}

//! Get TLE data.
vector< TwoLineElementData >& TwoLineElementsTextFileReader::getTwoLineElementData( )
{
    return twoLineElementData_;
}

//! Get number of objects.
unsigned int& TwoLineElementsTextFileReader::getNumberOfObjects( )
{
    return numberOfObjects_;
}

//! Convert and store TLE data.
void TwoLineElementsTextFileReader::storeTwoLineElementData( )
{
    // Calculate number of objects in catalog file
    numberOfObjects_ = static_cast< int >( lineCounter_ / 3 );

    // Set TLE data vector size.
    twoLineElementData_.resize( numberOfObjects_ );

    // Create vector of the three lines of a single object's TLE data as strings.
    vector< string > twoLineElementString_( 3 );

    // Create the object counter.
    unsigned int objectNumberCounter_ = 0;

    // Declare Keplerian elements variables.
    double inclination_;
    double rightAscensionOfAscendingNode_;
    double eccentricity_;
    double argumentOfPerigee_;
    double meanMotion_;

    // Declare approximate number of revolutions, remainder, and lost part.
    int approximateNumberOfRevolutions_;
    int approximateNumberOfRevolutionsRemainder_;
    int lostNumberOfRevolutions_;

    // Create Earth World Geodetic System (WGS) 72 gravity field.
    SphericalHarmonicsGravityField earthWorldGeodeticSystem72GravityField_;

    // Set predefined settings for World Geodetic System (WGS) 72 gravity field.
    earthWorldGeodeticSystem72GravityField_
            .setPredefinedSphericalHarmonicsGravityFieldSettings(
                SphericalHarmonicsGravityField::earthWorldGeodeticSystem72 );

    // Create Earth object and add WGS-72 gravity field.
    Planet earthWithWorldGeodeticSystem72GravityField_;
    earthWithWorldGeodeticSystem72GravityField_.setGravityFieldModel(
                &earthWorldGeodeticSystem72GravityField_ );

    // For every 3 lines, read the data from 3 consecutive strings of TLE data
    // and convert them to the TLE data variables
    for ( unsigned int i = 1; i < lineCounter_; i+=3 )
    {
        // General set-up for variable storing.
        //---------------------------------------------------------------------

        // Fill the vector with the line strings from the data container.
        twoLineElementString_.at( 0 ) = containerOfDataFromFile_.at( i );
        twoLineElementString_.at( 1 ) = containerOfDataFromFile_.at( i + 1 );
        twoLineElementString_.at( 2 ) = containerOfDataFromFile_.at( i + 2 );

        twoLineElementData_[ objectNumberCounter_ ].lineNumbers.push_back( i );
        twoLineElementData_[ objectNumberCounter_ ].lineNumbers.push_back( i + 1 );
        twoLineElementData_[ objectNumberCounter_ ].lineNumbers.push_back( i + 2 );

        // Initiate a stringstream for each line.
        stringstream line0StringStream_( stringstream::in | stringstream::out );
        stringstream line1StringStream_( stringstream::in | stringstream::out );
        stringstream line2StringStream_( stringstream::in | stringstream::out );

        // Insert the strings into the stringstreams
        line0StringStream_ << twoLineElementString_.at( 0 );
        line1StringStream_ << twoLineElementString_.at( 1 );
        line2StringStream_ << twoLineElementString_.at( 2 );

        // Line-0 variable storing.
        //---------------------------------------------------------------------

        // Declare string containing part of name of object.
        string namePart_;

        // Loop through the stringstream and read words that constitute name of
        // object. Store name parts in objectName storage container.
        while ( line0StringStream_ >> namePart_ )
        {
            twoLineElementData_[ objectNumberCounter_ ].objectName.push_back( namePart_ );
        }

        // Insert the entire line-0 string in the object name string.
        twoLineElementData_[ objectNumberCounter_ ].objectNameString
                = twoLineElementString_.at( 0 );

        // Line-1 variable storing.
        //---------------------------------------------------------------------
        // Fill all line-1 variables of the object structure with substrings
        // of the line-1 string using the template function to transform the
        // string to another type.
        // See reference for which columns are assigned to which variable.

        // Get line number integer of line-1 from string.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 0, 1 ),
            twoLineElementData_[ objectNumberCounter_ ].lineNumberLine1 );

        // Get object indentification number integer of line-1 from string
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 2, 5 ),
            twoLineElementData_[ objectNumberCounter_ ].objectIdentificationNumber );

        // Get classification character from string.
        twoLineElementData_[ objectNumberCounter_ ].tleClassification =
            twoLineElementString_.at( 1 ) [ 7 ];

        // Get launch year integer from string.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 9, 2 ),
            twoLineElementData_[ objectNumberCounter_ ].launchYear );

        // Calculate four-digit launch year from the above.
        if ( twoLineElementData_[ objectNumberCounter_ ].launchYear > 56 )
        {
            twoLineElementData_[ objectNumberCounter_ ].fourDigitlaunchYear =
                    twoLineElementData_[ objectNumberCounter_ ].launchYear + 1900;
        }

        else
        {
            twoLineElementData_[ objectNumberCounter_ ].fourDigitlaunchYear =
                    twoLineElementData_[ objectNumberCounter_ ].launchYear + 2000;
        }

        // Get launch number integer from string.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 11, 3 ),
                                 twoLineElementData_[ objectNumberCounter_ ].launchNumber );

        // Get launch part string from string.
        twoLineElementData_[ objectNumberCounter_ ].launchPart =
                twoLineElementString_.at( 1 ).substr( 14, 3 );

        // Get epoch year integer from string.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 18, 2 ),
            twoLineElementData_[ objectNumberCounter_ ].epochYear );

        // Calculate four-digit epoch year from the above.
        if ( twoLineElementData_[ objectNumberCounter_ ].epochYear > 56 )
        {
            twoLineElementData_[ objectNumberCounter_ ].fourDigitEpochYear =
                    twoLineElementData_[ objectNumberCounter_ ].epochYear + 1900;
        }

        else
        {
            twoLineElementData_[ objectNumberCounter_ ].fourDigitEpochYear =
                    twoLineElementData_[ objectNumberCounter_ ].epochYear + 2000;
        }

        // Get epoch day double from string.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 20, 12 ),
            twoLineElementData_[ objectNumberCounter_ ].epochDay );

        // Get "first-derivative of mean motion divided by two" double from string.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 33, 10 ),
            twoLineElementData_[ objectNumberCounter_].firstDerivativeOfMeanMotionDividedByTwo );

        // gGt coefficient of scientific notation of "second-derivative of mean motion divided
        // by six" double from string,
        // Apply implied leading decimal point.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 44, 6 ),
            twoLineElementData_[ objectNumberCounter_ ]
                                 .coefficientOfSecondDerivativeOfMeanMotionDividedBySix );
        twoLineElementData_[ objectNumberCounter_ ]
                .coefficientOfSecondDerivativeOfMeanMotionDividedBySix /= 100000;

        // Get exponent of scientific notation of "second-derivative of mean motion divided
        // by six" integer from string.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 50, 2 ),
            twoLineElementData_[ objectNumberCounter_]
                                 .exponentOfSecondDerivativeOfMeanMotionDividedBySix );

        // Calculate "second-derivative of mean motion divided by six" double from the above two.
        twoLineElementData_[ objectNumberCounter_ ].secondDerivativeOfMeanMotionDividedBySix =
                twoLineElementData_[ objectNumberCounter_]
                .coefficientOfSecondDerivativeOfMeanMotionDividedBySix
                * pow( 10, twoLineElementData_[ objectNumberCounter_ ]
                       .exponentOfSecondDerivativeOfMeanMotionDividedBySix );

        // Get coefficient of scientific notation of "B* divided by six" double
        // from string; apply implied leading decimal point.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 53, 6 ),
            twoLineElementData_[ objectNumberCounter_ ].coefficientOfBStar );
        twoLineElementData_[ objectNumberCounter_ ].coefficientOfBStar /= 100000;

        // Get exponent of scientific notation of B* integer from string
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 59, 2 ),
            twoLineElementData_[ objectNumberCounter_ ].exponentOfBStar );

        // Calculate B* double from the above two.
        twoLineElementData_[ objectNumberCounter_ ].bStar =
            twoLineElementData_[ objectNumberCounter_ ].coefficientOfBStar
                * pow( 10, twoLineElementData_[ objectNumberCounter_ ].exponentOfBStar);

        // Get orbital model integer from string.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 62, 1 ),
            twoLineElementData_[ objectNumberCounter_ ].orbitalModel );

        // Get TLE number integer from string.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 64, 4 ),
            twoLineElementData_[ objectNumberCounter_ ].tleNumber );

        // Get modulo-10 checksum integer from string.
        convertStringToTemplate( twoLineElementString_.at( 1 ).substr( 68, 1 ),
            twoLineElementData_[ objectNumberCounter_ ].modulo10CheckSumLine1 );

        // Line-2 variable storing.
        //---------------------------------------------------------------------
        // Fill all line-2 variables of the object structure partly with the
        // stringstream and partly with with substrings of the line-2 string
        // using the template function to transform the string to another type.
        // See reference for which columns are assigned to which variable.

        // Get line number integer of line-2 from stringstream.
        line2StringStream_ >> twoLineElementData_[ objectNumberCounter_ ].lineNumberLine2;

        // Get object identification number integer of line-2 from stringstream.
        line2StringStream_ >> twoLineElementData_[ objectNumberCounter_ ]
                              .objectIdentificationNumberLine2;

        // Get inclination double from stringstream.
        line2StringStream_ >> inclination_;
        twoLineElementData_[ objectNumberCounter_ ].TLEKeplerianElements
                .setInclination( inclination_ );

        // Get right ascension of ascending node double from stringstream.
        line2StringStream_ >> rightAscensionOfAscendingNode_;
        twoLineElementData_[ objectNumberCounter_ ].TLEKeplerianElements
                .setLongitudeOfAscendingNode( rightAscensionOfAscendingNode_ );

        // Get eccetricity double from stringstream.
        line2StringStream_ >> eccentricity_;
        eccentricity_ /= 10000000;
        twoLineElementData_[ objectNumberCounter_ ].TLEKeplerianElements
                .setEccentricity( eccentricity_ );

        // Get argument of perigee double from stringstream
        line2StringStream_ >> argumentOfPerigee_;
        twoLineElementData_[ objectNumberCounter_ ].TLEKeplerianElements
                .setArgumentOfPeriapsis( argumentOfPerigee_ );

        // Get mean anomaly double from stringstream.
        line2StringStream_ >> twoLineElementData_[ objectNumberCounter_ ].meanAnomaly;

        // Get mean motion double from line-2 string.
        convertStringToTemplate( twoLineElementString_[ 2 ].substr( 52, 11 ),
            twoLineElementData_[ objectNumberCounter_ ].meanMotionInRevolutionsPerDay );

        // Get revolution number integer from line-2 string.
        convertStringToTemplate( twoLineElementString_[ 2 ].substr( 63, 5 ),
            twoLineElementData_[ objectNumberCounter_ ].revolutionNumber );

        // Get modulo-10 checksum integer of line-2 from line-2 string.
        convertStringToTemplate( twoLineElementString_[ 2 ].substr( 68, 1 ),
            twoLineElementData_[ objectNumberCounter_ ].modulo10CheckSumLine2 );

        // Calculate the approximate total number of revolutions, as the counter resets to 0 after
        // passing by 99999, insert the current year in the following equation.
        approximateNumberOfRevolutions_ = twoLineElementData_[ objectNumberCounter_ ]
                .meanMotionInRevolutionsPerDay
                * ( currentYear_ - twoLineElementData_[ objectNumberCounter_ ].fourDigitlaunchYear )
                * PhysicalConstants::JULIAN_YEAR_IN_DAYS;
        approximateNumberOfRevolutionsRemainder_ = approximateNumberOfRevolutions_ % 100000;
        lostNumberOfRevolutions_ = approximateNumberOfRevolutions_
                - approximateNumberOfRevolutionsRemainder_;

        // Check if the counter has been reset after passing 99999.
        if ( ( twoLineElementData_[ objectNumberCounter_ ].revolutionNumber -
               approximateNumberOfRevolutionsRemainder_ ) <=
             ( approximateNumberOfRevolutionsRemainder_ + 100000 -
               twoLineElementData_[ objectNumberCounter_ ].revolutionNumber ) )
        {
            twoLineElementData_[ objectNumberCounter_ ].totalRevolutionNumber =
                    lostNumberOfRevolutions_ + twoLineElementData_[ objectNumberCounter_ ]
                    .revolutionNumber;
        }

        else
        {
            twoLineElementData_[ objectNumberCounter_ ].totalRevolutionNumber =
                    lostNumberOfRevolutions_ - 100000 + twoLineElementData_[ objectNumberCounter_ ].
                    revolutionNumber;
        }

        // Check if total number of revolutions is negative.
        if ( twoLineElementData_[ objectNumberCounter_ ].totalRevolutionNumber < 0 )
        {
            twoLineElementData_[ objectNumberCounter_ ].totalRevolutionNumber =
                    twoLineElementData_[ objectNumberCounter_ ].revolutionNumber;
        }

        // Semi-major axis of the object is calculated from the other TLE variables.
        meanMotion_ = twoLineElementData_[ objectNumberCounter_ ].meanMotionInRevolutionsPerDay
                 * 2.0 * M_PI / PhysicalConstants::JULIAN_DAY;
        twoLineElementData_[ objectNumberCounter_ ].TLEKeplerianElements.setSemiMajorAxis(
                    orbital_element_conversions::convertMeanMotionToSemiMajorAxis(
                        meanMotion_, &earthWithWorldGeodeticSystem72GravityField_ ) );

        // Perigee of the object is calculated from the other TLE variables.
        twoLineElementData_[ objectNumberCounter_ ].perigee =
                twoLineElementData_[ objectNumberCounter_ ].TLEKeplerianElements.getSemiMajorAxis( )
                * ( 1.0 - twoLineElementData_[ objectNumberCounter_ ].TLEKeplerianElements
                    .getEccentricity( ) );

        // Apogee of the object is calculated from the other TLE variables.
        twoLineElementData_[ objectNumberCounter_ ].apogee =
                twoLineElementData_[ objectNumberCounter_ ].TLEKeplerianElements.getSemiMajorAxis( )
                * ( 1.0 + twoLineElementData_[ objectNumberCounter_ ].TLEKeplerianElements
                    .getEccentricity( ) );

        // Increment object number counter by one.
        objectNumberCounter_++;
    }
}

//! Checks the integrity of the TLE input file.
multimap< int, string > TwoLineElementsTextFileReader::checkTwoLineElementsFileIntegrity( )
{
    // Create vector of the three lines of a TLE as strings
    vector< string > twoLineElementString_( 3 );

    // Boolean which turns to true if one of the tests is not passed.
    bool isObjectErroneous = false;

    // Vector of object numbers of objects with corrupted TLE data.
    vector< int > corruptedTwoLineElementDataPositions_;

    // Multimap of corrupted TLE errors.
    multimap< int, string > corruptedTwoLineElementDataErrors_;

    // Loop is started over the entire structure of objects,
    // multiple checks are performed per TLE, if a check is not passed boolean
    // is set to true for that object, and amount of corrupted TLEs integer is
    // incremented.
    for ( unsigned int i = 0; i < numberOfObjects_; i++ )
    {
        // Fill the vector with the line strings from the data container.
        twoLineElementString_.at( 0 ) = containerOfDataFromFile_[ i + ( i * 2 ) + 1 ];
        twoLineElementString_.at( 1 ) = containerOfDataFromFile_[ i + ( i * 2 ) + 2 ];
        twoLineElementString_.at( 2 ) = containerOfDataFromFile_[ i + ( i * 2 ) + 3 ];

        // Set boolean to false.
        isObjectErroneous = false;

        // If a line-1 line number is not 1, error message is given and boolean is set to true.
        if ( twoLineElementData_.at( i ).lineNumberLine1 != 1 )
        {
            corruptedTwoLineElementDataErrors_.insert(
                        pair< int, string >( i, "Incorrect line-1 leading integer." ) );

            isObjectErroneous = true;
        }

        // If a line-2 line number is not 2, error message is given and boolean is set to true.
        if ( twoLineElementData_.at( i ).lineNumberLine2 != 2 )
        {
            corruptedTwoLineElementDataErrors_.insert(
                        pair< int, string >( i, "Incorrect line-2 leading integer." ) );

            isObjectErroneous = true;
        }

        if ( twoLineElementData_.at( i ).tleClassification != 'U'
                  && twoLineElementData_.at( i ).tleClassification != 'C' )
        {
            corruptedTwoLineElementDataErrors_.insert(
                        pair< int, string >( i, "Invalid TLE classification." ) );

            isObjectErroneous = true;
        }

        // Check if orbital model is 0, else error message is given and boolean is set to true.
        if ( twoLineElementData_.at( i ).orbitalModel != 0 )
        {
            corruptedTwoLineElementDataErrors_.insert(
                        pair< int, string >( i, "Incorrect orbital model." ) );

            isObjectErroneous = true;
        }

        // Line-1 modulo-10 checker, add all integers of line 1 together,
        // and minus signs count as 1, rest is not counted.

        // Create calculated modulo-10 checksum of line 1 and temporary integer for addition.
        unsigned int line1Modulo10Sum_ = 0;
        int temporaryInteger_ = 0;

        // Run loop over the 69-character long, line-2 string, but not the last character,
        // since this is the modulo-10 checksum itself and is not added.
        for ( unsigned int j = 0; j < 68; j++ )
        {
            // Check if the char is not a space or a dot, or a plus or a minus sign or any of the
            // alphabetical positions, because these can not be added.
            if ( twoLineElementString_.at( 1 )[ j ] != ' ' &&
                 twoLineElementString_.at( 1 )[ j ] != '.' &&
                 twoLineElementString_.at( 1 )[ j ] != '+' &&
                 twoLineElementString_.at( 1 )[ j ] != '-' &&
                 j != 7 &&
                 j != 14 &&
                 j != 15 &&
                 j != 16 )
            {
                // Convert char to int.
                convertStringToTemplate( twoLineElementString_.at( 1 ).substr( j , 1 ),
                                         temporaryInteger_ );

                // Add int to checksum.
                line1Modulo10Sum_ += temporaryInteger_;
            }

            // Check is the char is a minus sign
            else if ( twoLineElementString_.at( 1 )[ j ] == '-' )
            {
                // Add 1 to checksum.
                line1Modulo10Sum_++;
            }
        }

        // Take modulo 10 of checksum.
        line1Modulo10Sum_ %= 10;

        // Comparison of calculated modulo-10 checksum with modulo-10 checksum from TLE for line 1,
        // if they do not compare, error message is given and boolean is set to true.
        if ( line1Modulo10Sum_ != twoLineElementData_.at( i ).modulo10CheckSumLine1 )
        {
            corruptedTwoLineElementDataErrors_.insert(
                        pair< int, string >( i, "Incorrect line-1 modulo-10 checksum." ) );

            isObjectErroneous = true;
        }

        // Line-2 modulo-10 checker, add all integers of line 2 together.
        unsigned int line2Modulo10Sum_ = 0;

        // Run loop over the 69-character long, line-2 string, but not the last character,
        // since this is the modulo-10 checksum itself and is not added.
        for ( unsigned int k = 0; k < 68; k++ )
        {
            // Check if the char is not a space or a dot, because these can not be added.
            if ( twoLineElementString_.at( 2 )[ k ] != ' ' &&
                 twoLineElementString_.at( 2 )[ k ] != '.' )
            {
                // Convert char to int
                convertStringToTemplate( twoLineElementString_[ 2 ].substr( k , 1 ),
                                         temporaryInteger_ );

                // Add int to checksum.
                line2Modulo10Sum_ += temporaryInteger_;
            }
        }

        // Take modulo 10 of checksum.
        line2Modulo10Sum_ %= 10;

        // Comparison of calculated modulo-10 checksum with modulo-10 checksum from TLE for line 2,
        // if they do not compare error message is given and boolean is set to true.
        if ( line2Modulo10Sum_ != twoLineElementData_.at( i ).modulo10CheckSumLine2 )
        {
            corruptedTwoLineElementDataErrors_.insert(
                        pair< int, string >( i, "Incorrect line-2 modulo-10 checksum." ) );

            isObjectErroneous = true;
        }

        // When the object's line-1 and line-2 identification numbers do not match an
        // error message is given and boolean is set to true.
        if ( twoLineElementData_.at( i ).objectIdentificationNumber !=
             twoLineElementData_.at( i ).objectIdentificationNumberLine2 )
        {
            corruptedTwoLineElementDataErrors_.insert(
                        pair< int, string >(
                            i, "Line-1 and line-2 object idenfitication number mismatch." ) );

            isObjectErroneous = true;
        }

        // If object is erroneous, add to list.
        if ( isObjectErroneous )
        {
            corruptedTwoLineElementDataPositions_.push_back( i );
        }
    }

    // Erase all the corrupted TLEs from the TLE data vector.

    for ( unsigned int q = 0; q < corruptedTwoLineElementDataPositions_.size( ); q++ )
    {
        twoLineElementData_.erase( twoLineElementData_.begin( ) +
                                   corruptedTwoLineElementDataPositions_.at( q ) - q );
        numberOfObjects_--;
    }

    // Return the amount of corrupted TLEs
    return corruptedTwoLineElementDataErrors_;
}

// End of file.
