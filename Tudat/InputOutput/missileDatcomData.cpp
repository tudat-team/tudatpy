/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Blake, W.B. Missile Datcom User's Manual - 1997 Fortran 90 Version, AFRL-VA-WP-TR-1998-3009
 *          Air Force Research Laboratory, 1998.
 *
 */

#include <iostream>
#include <iterator>
#include <sstream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/InputOutput/missileDatcomData.h"
#include "Tudat/InputOutput/missileDatcomReader.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace input_output
{

using unit_conversions::convertDegreesToRadians;
using std::iterator;

//! Constructor that reads and processes Missile Datcom output.
MissileDatcomData::MissileDatcomData( const std::string& fileNameAndPath )
{
    MissileDatcomReader myMissileDatcomReader( fileNameAndPath );
    convertDatcomData( myMissileDatcomReader.getMissileDatcomData( ) );
}

//! Function to convert the MissileDatcomData.
void MissileDatcomData::convertDatcomData( const std::vector< double >& datcomData )
{
    std::vector< double > datcomDataToBeConverted = datcomData;

    // Add an extra item at the begin, such that the first real value is at entry 1.
    // The same indices as described in the MissileDatcom user manual can now be used.
    std::vector< double >::iterator it;
    it = datcomDataToBeConverted.begin( );
    datcomDataToBeConverted.insert( it,-0.0 );

    // The first "card" is the Flight condition data array, which is the same for every "case".
    numberOfAnglesOfAttack_ = datcomDataToBeConverted[ 1 ];

    for ( int iteratorDataVector_= 2; iteratorDataVector_ < ( 2 + numberOfAnglesOfAttack_ );
          iteratorDataVector_++ )
    {
        angleOfAttack_.push_back( datcomDataToBeConverted[ iteratorDataVector_ ] );
    }

    sideslipAngle_ = datcomDataToBeConverted[ 22 ];

    rollAngle_ = datcomDataToBeConverted[ 23 ];

    numberOfMachNumbers_ = datcomDataToBeConverted[ 24 ];

    for ( int iteratorDataVector_= 25; iteratorDataVector_ < ( 25 + numberOfMachNumbers_ );
          iteratorDataVector_++ )
    {
        machNumber_.push_back( datcomDataToBeConverted[ iteratorDataVector_ ] );
    }

    for ( int iteratorDataVector_= 45; iteratorDataVector_ < ( 25 + numberOfMachNumbers_ );
          iteratorDataVector_++ )
    {
        altitudes_.push_back( datcomDataToBeConverted[ iteratorDataVector_ ] );
    }

    for ( int iteratorDataVector_= 66; iteratorDataVector_ < ( 66 + numberOfMachNumbers_ );
          iteratorDataVector_++ )
    {
        reynoldsNumbers_.push_back( datcomDataToBeConverted[ iteratorDataVector_ ] );
    }

    for ( int iteratorDataVector_= 86; iteratorDataVector_ < (86 +numberOfMachNumbers_);
          iteratorDataVector_++ )
    {
        freeStreamVelocities_.push_back( datcomDataToBeConverted[ iteratorDataVector_ ] );
    }

    for ( int iteratorDataVector_= 106; iteratorDataVector_ < (106 + numberOfMachNumbers_);
          iteratorDataVector_++ )
    {
        freeStreamStaticTemperatures_.push_back( datcomDataToBeConverted[ iteratorDataVector_ ] );
    }

    // NOTE apperently the last entry of flight card is not written to de for004.dat file.
    // According to the user manual the last entry should be 145, this entry is not found in the
    // for004.dat file.
    for ( int iteratorDataVector_= 126; iteratorDataVector_ <= 144; iteratorDataVector_++ )
    {
        freeStreamStaticPressure_.push_back( datcomDataToBeConverted[ iteratorDataVector_ ] );
    }

    // Loop over Mach, angle of attack and coefficient index to fill the static coefficient array.
    for (int machIndex = 0; machIndex < numberOfMachNumbers_; machIndex++ )
    {
        for (int angleOfAttackIndex = 0; angleOfAttackIndex < numberOfAnglesOfAttack_;
            angleOfAttackIndex++ )
        {
            for (int coefficientIndex = 0; coefficientIndex < 6; coefficientIndex++ )
            {
                // The second index begins at entry 145
                // The first section has 144 entries, the second 220 and the last 400. Therefore
                // there is a repetition after 144+220+400=764 entries.
                // every coefficients has 20 entries.
                int dataVectorIndex_ = 145 + 764 * machIndex + angleOfAttackIndex
                                   + 20 * coefficientIndex;
                staticCoefficients_[ machIndex ][ angleOfAttackIndex ][ coefficientIndex ] =
                        datcomDataToBeConverted[ dataVectorIndex_ ];
            }

            for ( int coefficientIndex = 6; coefficientIndex < 11; coefficientIndex++ )
            {
                // The second index begins at entry 145
                // The first section has 144 entries, the second 220 and the last 400. Therefore
                // there is a repetition after 144+220+400=764 entries.
                // every coefficients has 20 entries.
                // The coefficients are in per degree and needs to be converted to per radian
                int dataVectorIndex_ = 145 + 764 * machIndex + angleOfAttackIndex
                                   + 20 * coefficientIndex;
                staticCoefficients_[ machIndex ][ angleOfAttackIndex ][ coefficientIndex ] =
                        datcomDataToBeConverted[ dataVectorIndex_ ]
                        / convertDegreesToRadians( 1.0 );
            }
        }
    }

    // Loop over Mach, angle of attack and coefficient index to fill the dynamic coefficient array.
    for ( int machIndex = 0; machIndex < numberOfMachNumbers_; machIndex++ )
    {
        for ( int angleOfAttackIndex = 0; angleOfAttackIndex < numberOfAnglesOfAttack_;
            angleOfAttackIndex++ )
        {
            for ( int coefficientIndex = 0; coefficientIndex < 20; coefficientIndex++ )
            {
                // The third index begins at entry 365.
                // The first section has 144 entries, the second 220 and the last 400. Therefore
                // there is a repetition after 144+220+400=764 entries.
                // every coefficients has 20 entries.
                // WARNINNG: possible error in datcom user manual. It looks like to values are
                // always in per radian (instead per degree, which is mentioned in the manual).
                int dataVectorIndex_ = 365 + 764 * machIndex + angleOfAttackIndex
                                   + 20 * coefficientIndex;
                dynamicCoefficients_[ machIndex ][ angleOfAttackIndex ][ coefficientIndex ] =
                     datcomDataToBeConverted[ dataVectorIndex_ ];
            }
        }
    }
}

//! Function to access the static coefficient database.
double MissileDatcomData::getStaticCoefficient( int machIndex, int angleOfAttackIndex,
                                                StaticCoefficientNames coefficientIndex )
{
    return staticCoefficients_[ machIndex ][ angleOfAttackIndex ][ coefficientIndex ];
}

//! Function to access the dynamic coefficient database.
double MissileDatcomData::getDynamicCoefficient( int machIndex, int angleOfAttackIndex,
                                                 DynamicCoefficientNames coefficientIndex )
{
    return dynamicCoefficients_[ machIndex ][ angleOfAttackIndex ][ coefficientIndex ];
}

//! Write the database to CSV files.
void MissileDatcomData::writeCoefficientsToFile( const std::string& fileNameBase,
                                                 const int basePrecision,
                                                 const int exponentWidth )
{
    // Write coefficients to a file spearately for each angle of attack.
    for ( unsigned int i = 0; i < angleOfAttack_.size( ) ; i++ )
    {
        // Make filename for specific angle of attack.
        std::ostringstream stringstream;
        stringstream << fileNameBase << "_" << i;
        std::string fileName = stringstream.str( );

        // Open output file.
        std::ofstream outputFile;
        outputFile.open( fileName.c_str( ) );
        outputFile.precision( 10 );

        // Iterate over all Mach numbers in database.
        for ( unsigned int j = 0; j < machNumber_.size( ) ; j++ )
        {
            // Write angle of attack and mach number.
            outputFile << angleOfAttack_[ i ] << " " << machNumber_[ j ] << " ";

            // Write static coefficients.
            for ( int k = 0; k < 11; k++ )
            {
                outputFile << printInFormattedScientificNotation(
                                  staticCoefficients_[ j ][ i ][ k ], basePrecision,
                                  exponentWidth ) << " ";
            }

            // Write dynamic coefficients.
            for ( int k = 0; k < 19; k++ )
            {

                outputFile << printInFormattedScientificNotation(
                                  dynamicCoefficients_[ j ][ i ][ k ], basePrecision,
                                  exponentWidth ) << " ";
            }

            outputFile << printInFormattedScientificNotation(
                              dynamicCoefficients_[ j ][ i ][ 19 ], basePrecision,
                              exponentWidth ) << std::endl;
        }

        outputFile.close( );
    }
}

} // namespace input_output
} // namespace tudat
