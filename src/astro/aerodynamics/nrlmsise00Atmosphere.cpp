/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/aerodynamics/nrlmsise00Atmosphere.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include <iostream>

//! Tudat library namespace.
namespace tudat
{
namespace aerodynamics
{

void NRLMSISE00Atmosphere::computeProperties(
        const double altitude, const double longitude,
        const double latitude, const double time )
{
    // Compute the hash key
    size_t hashKey = hashFunc( altitude, longitude, latitude, time );

    // If hash key is same do nothing
    if (hashKey == hashKey_)
    {
        return;
    }
    hashKey_ = hashKey;

    // Retrieve input data.
    inputData_ = nrlmsise00InputFunction_(
                altitude, longitude, latitude, time );
    std::copy( inputData_.apVector.begin( ), inputData_.apVector.end( ), aph_.a );
    std::copy( inputData_.switches.begin( ), inputData_.switches.end( ), flags_.switches);

    input_.g_lat  = latitude * 180.0 / mathematical_constants::PI; // rad to deg
    input_.g_long = longitude * 180.0 / mathematical_constants::PI; // rad to deg
    input_.alt    = altitude * 1.0E-3; // m to km
    input_.year   = inputData_.year;
    input_.doy    = inputData_.dayOfTheYear;
    input_.sec    = inputData_.secondOfTheDay;
    input_.lst    = inputData_.localSolarTime;
    input_.f107   = inputData_.f107;
    input_.f107A  = inputData_.f107a;
    input_.ap     = inputData_.apDaily;
    input_.ap_a   = &aph_;

    // Call NRLMSISE00
    gtd7(&input_, &flags_, &output_);

    // Retrieve density and temperature
    density_ = output_.d[ 5 ] * 1000.0; // GM/CM3 to kg/M3
    temperature_ = output_.t[1];

    // Get number densities
    numberDensities_.resize(8);
    numberDensities_[0] = output_.d[0] * 1.0E6 ; // HE NUMBER DENSITY    (M-3)
    numberDensities_[1] = output_.d[1] * 1.0E6 ; // O NUMBER DENSITY     (M-3)
    numberDensities_[2] = output_.d[2] * 1.0E6 ; // N2 NUMBER DENSITY    (M-3)
    numberDensities_[3] = output_.d[3] * 1.0E6 ; // O2 NUMBER DENSITY    (M-3)
    numberDensities_[4] = output_.d[4] * 1.0E6 ; // AR NUMBER DENSITY    (M-3)
    numberDensities_[5] = output_.d[6] * 1.0E6 ; // H NUMBER DENSITY     (M-3)
    numberDensities_[6] = output_.d[7] * 1.0E6 ; // N NUMBER DENSITY     (M-3)
    numberDensities_[7] = output_.d[8] * 1.0E6 ; // Anomalous oxygen NUMBER DENSITY  (M-3)

    // Get average number density
    double sumOfNumberDensity = 0.0 ;
    for( unsigned int i = 0 ; i < numberDensities_.size( ) ; i++)
    {
        sumOfNumberDensity += numberDensities_[ i ];
    }
    averageNumberDensity_ = sumOfNumberDensity / double( numberDensities_.size( ) );

    // Mean molar mass (Thermodynamics an Engineering Approach, Michael A. Boles)
    meanMolarMass_ = numberDensities_[0] * gasComponentProperties_.molarMassHelium;
    meanMolarMass_ += numberDensities_[1] * gasComponentProperties_.molarMassAtomicOxygen;
    meanMolarMass_ += numberDensities_[2] * gasComponentProperties_.molarMassNitrogen;
    meanMolarMass_ += numberDensities_[3] * gasComponentProperties_.molarMassOxygen;
    meanMolarMass_ += numberDensities_[4] * gasComponentProperties_.molarMassArgon;
    meanMolarMass_ += numberDensities_[5] * gasComponentProperties_.molarMassAtomicHydrogen;
    meanMolarMass_ += numberDensities_[6] * gasComponentProperties_.molarMassAtomicNitrogen;
    meanMolarMass_ += numberDensities_[7] * gasComponentProperties_.molarMassOxygen;
    meanMolarMass_ = meanMolarMass_ / sumOfNumberDensity ;

    // Speed of sound
    speedOfSound_ = aerodynamics::computeSpeedOfSound(
                temperature_, specificHeatRatio_, molarGasConstant_ / meanMolarMass_ );

    // Collision diameter
    weightedAverageCollisionDiameter_ = numberDensities_[0]* gasComponentProperties_.diameterHelium ;
    weightedAverageCollisionDiameter_ += numberDensities_[1]* gasComponentProperties_.diameterAtomicOxygen ;
    weightedAverageCollisionDiameter_ += numberDensities_[2]* gasComponentProperties_.diameterNitrogen ;
    weightedAverageCollisionDiameter_ += numberDensities_[3]* gasComponentProperties_.diameterOxygen ;
    weightedAverageCollisionDiameter_ += numberDensities_[4]* gasComponentProperties_.diameterArgon ;
    weightedAverageCollisionDiameter_ += numberDensities_[5]* gasComponentProperties_.diameterAtomicHydrogen ;
    weightedAverageCollisionDiameter_ += numberDensities_[6]* gasComponentProperties_.diameterAtomicNitrogen ;
    weightedAverageCollisionDiameter_ += numberDensities_[7]* gasComponentProperties_.diameterAtomicOxygen ;
    weightedAverageCollisionDiameter_ = weightedAverageCollisionDiameter_ / sumOfNumberDensity;

    // Mean free path.
    meanFreePath_ = aerodynamics::computeMeanFreePath( weightedAverageCollisionDiameter_, averageNumberDensity_ );

    // Calculate pressure using ideal gas law (Thermodynamics an Engineering Approach, Michael A. Boles)
    if( useIdealGasLaw_ )
    {
        pressure_ = density_ * molarGasConstant_ * temperature_ / meanMolarMass_ ;
    }
    else
    {
        pressure_ = TUDAT_NAN;
    }
}

//! Overloaded ostream to print class information.
std::ostream& operator << ( std::ostream& stream,
                            NRLMSISE00Input& nrlmsiseInput ){
    stream << "This is a NRLMSISE Input data object." << std::endl;
    stream << "The input data is stored as: " << std::endl;

    stream << "Year              = " << nrlmsiseInput.year << std::endl;
    stream << "Day of the year   = " << nrlmsiseInput.dayOfTheYear << std::endl;
    stream << "Second of the day = " << nrlmsiseInput.secondOfTheDay << std::endl;
    stream << "Local solar time  = " << nrlmsiseInput.localSolarTime << std::endl;
    stream << "f107              = " << nrlmsiseInput.f107 << std::endl;
    stream << "f107a             = " << nrlmsiseInput.f107a << std::endl;
    stream << "apDaily           = " << nrlmsiseInput.apDaily << std::endl;

    for( unsigned int i = 0 ; i < nrlmsiseInput.apVector.size( ) ; i++ )
    {
        stream << "apVector[ " << i << " ]     = " << nrlmsiseInput.apVector[i] << std::endl;
    }

    for( unsigned int i = 0 ; i < nrlmsiseInput.switches.size( ) ; i++ )
    {
        stream << "switches[ " << i << " ]     = " << nrlmsiseInput.switches[i] << std::endl;
    }

    return stream;
}

//! Get the full model output
std::pair< std::vector< double >, std::vector< double > >
NRLMSISE00Atmosphere::getFullOutput( const double altitude, const double longitude,
                                     const double latitude, const double time )
{
    // Compute the properties
    computeProperties( altitude, longitude, latitude, time );
    std::pair< std::vector< double >, std::vector< double >> output;

    // Copy array members of struct to vectors on the pair.
    output.first = std::vector< double >(
                output_.d, output_.d + sizeof output_.d / sizeof output_.d[ 0 ] );
    output.second = std::vector< double >(
                output_.t, output_.t + sizeof output_.t / sizeof output_.t[ 0 ] );
    return output;
}

}  // namespace aerodynamics
}  // namespace tudat
