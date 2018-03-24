/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/External/SofaInterface/fundamentalArguments.h"

namespace tudat
{

namespace sofa_interface
{

//! Function to calculate mean elongation of the Moon from the Sun
double calculateMeanElongationOfMoonFromSun( const double julianCenturiesSinceJ2000 )
{
    return iauFad03( julianCenturiesSinceJ2000 );
}

//! Function to calculate mean argument of latitude of the Moon
double calculateMeanArgumentOfLatitudeOfMoon( const double julianCenturiesSinceJ2000 )
{
    return iauFaf03( julianCenturiesSinceJ2000 );
}

//! Function to calculate mean anomaly of the Moon
double calculateMeanAnomalyOfMoon( const double julianCenturiesSinceJ2000 )
{
    return iauFal03( julianCenturiesSinceJ2000 );
}

//! Function to calculate mean anomaly of the Sun
double calculateMeanAnomalyOfSun( const double julianCenturiesSinceJ2000 )
{
    return iauFalp03( julianCenturiesSinceJ2000 );
}

//! Function to calculate longitude of the Moon's ascending node
double calculateLongitudeOfMoonsAscendingNode( const double julianCenturiesSinceJ2000 )
{
    return iauFaom03( julianCenturiesSinceJ2000 );
}

//! Function to calculate the Delaunay fundamental arguments at the requested time.
Eigen::Matrix< double, 5, 1 > calculateDelaunayFundamentalArguments(
        const double tdbTime )
{
    // Calculate current centuries since J2000.
    double julianCenturiesSinceJ2000 = basic_astrodynamics::convertSecondsSinceEpochToJulianCenturiesSinceEpoch( tdbTime );

    return ( Eigen::Matrix< double, 5, 1 >( ) <<
             calculateMeanAnomalyOfMoon( julianCenturiesSinceJ2000 ),
             calculateMeanAnomalyOfSun( julianCenturiesSinceJ2000 ),
             calculateMeanArgumentOfLatitudeOfMoon( julianCenturiesSinceJ2000 ),
             calculateMeanElongationOfMoonFromSun( julianCenturiesSinceJ2000 ),
             calculateLongitudeOfMoonsAscendingNode( julianCenturiesSinceJ2000 ) ).finished( );
}

//! Function to calculate the Delaunay fundamental arguments and (GMST + pi) at the requested time.
Eigen::Vector6d  calculateDelaunayFundamentalArgumentsWithGmst(
        const double tdbTime, const double terrestrialTime, const double universalTime1 )
{
    // Calculate GMST
    double greenwichMeanSiderealTime = iauGmst06(
                basic_astrodynamics::JULIAN_DAY_ON_J2000, universalTime1 / physical_constants::JULIAN_DAY,
                basic_astrodynamics::JULIAN_DAY_ON_J2000, terrestrialTime / physical_constants::JULIAN_DAY );

    // Calculate fundamental arguments.
    Eigen::Vector6d fundamentalArguments;
    fundamentalArguments << greenwichMeanSiderealTime + mathematical_constants::PI,
            calculateDelaunayFundamentalArguments( tdbTime );

    return fundamentalArguments;
}

//! Function to calculate the Delaunay fundamental arguments and (GMST + pi) at the requested time.
Eigen::Vector6d  calculateDelaunayFundamentalArgumentsWithGmst(
        const double tdbTime )
{
    return calculateDelaunayFundamentalArgumentsWithGmst( tdbTime, tdbTime, convertTTtoUTC( tdbTime ) );
}

//! Function to calculate the Doodson arguments at the requested time.
Eigen::Vector6d calculateDoodsonFundamentalArguments(
        const double tdbTime, const double terrestrialTime, const double universalTime1 )
{
    using mathematical_constants::PI;

    // Calculate current centuries since J2000.
    double julianCenturiesSinceJ2000 = basic_astrodynamics::convertSecondsSinceEpochToJulianCenturiesSinceEpoch( tdbTime );

    // Calculate GMST
    double greenwichMeanSiderealTime = iauGmst06(
                basic_astrodynamics::JULIAN_DAY_ON_J2000, universalTime1 / physical_constants::JULIAN_DAY,
                basic_astrodynamics::JULIAN_DAY_ON_J2000, terrestrialTime / physical_constants::JULIAN_DAY );

    // Calculate relevant angles.
    double meanElongationOfMoonFromSun = calculateMeanElongationOfMoonFromSun( julianCenturiesSinceJ2000 );
    double meanArgumentOfLatitudeOfMoon = calculateMeanArgumentOfLatitudeOfMoon( julianCenturiesSinceJ2000 );
    double meanAnomalyOfMoon = calculateMeanAnomalyOfMoon( julianCenturiesSinceJ2000 );
    double meanAnomalyOfSun = calculateMeanAnomalyOfSun( julianCenturiesSinceJ2000 );
    double longitudeOfMoonsAscendingNode = calculateLongitudeOfMoonsAscendingNode( julianCenturiesSinceJ2000 );

    // Calculate Doodson arguments.
    Eigen::Vector6d doodsonArguments;
    doodsonArguments( 1 ) = meanArgumentOfLatitudeOfMoon + longitudeOfMoonsAscendingNode;
    doodsonArguments( 0 ) = greenwichMeanSiderealTime + PI - doodsonArguments( 1 );
    doodsonArguments( 2 ) = doodsonArguments( 1 ) - meanElongationOfMoonFromSun;
    doodsonArguments( 3 ) = doodsonArguments( 1 ) - meanAnomalyOfMoon;
    doodsonArguments( 5 ) = doodsonArguments( 2 ) - meanAnomalyOfSun;
    doodsonArguments( 4 ) = -longitudeOfMoonsAscendingNode;

    return doodsonArguments;
}

//! Function to calculate the Doodson arguments at the requested time.
Eigen::Vector6d calculateDoodsonFundamentalArguments( const double tdbTime )
{
    return calculateDoodsonFundamentalArguments( tdbTime, tdbTime,  convertTTtoUTC( tdbTime ) );
}

} // namespace sofa_interfaces

} // namespace tudat
