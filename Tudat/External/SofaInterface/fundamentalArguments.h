/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SOFAFUNDAMENTALARGUMENTS_H
#define TUDAT_SOFAFUNDAMENTALARGUMENTS_H

extern "C"
{
    #include <sofa/src/sofa.h>
    #include <sofa/src/sofam.h>
}

#include <vector>

#include <Eigen/Core>

#include "Tudat/External/SofaInterface/sofaTimeConversions.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace sofa_interface
{

//! Matrix to convert 5 Doodson arguments (without GMST) to Delaunay arguments.
/*!
 * Matrix to convert 5 Doodson arguments (without GMST) to Delaunay arguments. For order:
 * \sa calculateDoodsonFundamentalArguments
 * \sa calculateDelaunayFundamentalArguments
 */
static const Eigen::Matrix< double, 5, 5 > delaunayToDoodsonArguments =
        ( Eigen::Matrix< double, 5, 5 >( ) <<
           0.0,  0.0, 1.0,  0.0,  1.0,
           0.0,  0.0, 1.0, -1.0,  1.0,
          -1.0,  0.0, 1.0,  0.0,  1.0,
           0.0,  0.0, 0.0,  0.0, -1.0,
           0.0, -1.0, 1.0, -1.0,  1.0 ).finished( );

//! Matrix to convert delaunay arguments to 5 Doodson arguments (without GMST).
/*!
 * Matrix to convert delaunay arguments to 5 Doodson arguments (without GMST). For order:
 * \sa calculateDoodsonFundamentalArguments
 * \sa calculateDelaunayFundamentalArguments
 */
static const Eigen::Matrix< double, 5, 5 > doodsonToDelaunayArguments =
        ( Eigen::Matrix< double, 5, 5 >( ) <<
          1.0,  0.0, -1.0,  0.0,  0.0,
          0.0,  1.0,  0.0,  0.0, -1.0,
          1.0,  0.0,  0.0,  1.0,  0.0,
          1.0, -1.0,  0.0,  0.0,  0.0,
          0.0,  0.0,  0.0, -1.0,  0.0 ).finished( );

//! Function to calculate mean elongation of the Moon from the Sun
/*!
 *  Function to calculate mean elongation of the Moon from the Sun fundamental argument using Sofa implementation
 *  \param julianCenturiesSinceJ2000 Julian centuries (period of 36525 days of 86400 seconds) since J2000
 *  \return Mean elongation of the Moon from the Sun
 */
double calculateMeanElongationOfMoonFromSun( const double julianCenturiesSinceJ2000 );

//! Function to calculate mean argument of latitude of the Moon
/*!
 *  Function to calculate mean argument of latitude of the Moon fundamental argument using Sofa implementation
 *  \param julianCenturiesSinceJ2000 Julian centuries (period of 36525 days of 86400 seconds) since J2000
 *  \return Mean argument of latitude of the Moon
 */
double calculateMeanArgumentOfLatitudeOfMoon( const double julianCenturiesSinceJ2000 );

//! Function to calculate mean anomaly of the Moon
/*!
 *  Function to calculate mean anomaly of the Moon fundamental argument using Sofa implementation
 *  \param julianCenturiesSinceJ2000 Julian centuries (period of 36525 days of 86400 seconds) since J2000
 *  \return Mean anomaly of the Moon
 */
double calculateMeanAnomalyOfMoon( const double julianCenturiesSinceJ2000 );

//! Function to calculate mean anomaly of the Sun
/*!
 *  Function to calculate mean anomaly of the Sun fundamental argument using Sofa implementation
 *  \param julianCenturiesSinceJ2000 Julian centuries (period of 36525 days of 86400 seconds) since J2000
 *  \return Mean anomaly of the Sun
 */
double calculateMeanAnomalyOfSun( const double julianCenturiesSinceJ2000 );

//! Function to calculate longitude of the Moon's ascending node
/*!
 *  Function to calculate longitude of the Moon's ascending node fundamental argument using Sofa implementation
 *  \param julianCenturiesSinceJ2000 Julian centuries (period of 36525 days of 86400 seconds) since J2000
 *  \return Longitude of the Moon's ascending node
 */
double calculateLongitudeOfMoonsAscendingNode( const double julianCenturiesSinceJ2000 );

//! Function to calculate the Delaunay fundamental arguments at the requested time.
/*!
 * Function to calculate the Delaunay fundamental arguments at the requested time, using Sofa implementation of
 * angle calculations. Using the typical notation of teh angles, the output order is: l, l', F, D, Omega.
 * \param tdbTime Time in TDB at which the arguments are to be calculated.
 * \return Delaunay fundamental arguments at the requested time.
 */
Eigen::Matrix< double, 5, 1 > calculateDelaunayFundamentalArguments(
        const double tdbTime );

//! Function to calculate the Delaunay fundamental arguments and (GMST + pi) at the requested time.
/*!
 * Function to calculate the Delaunay fundamental arguments and (GMST + pi) at the requested time, using Sofa implementation
 * of angle calculations.
 * Using the typical notation of teh angles, the output order is: (GMST + pi), l, l', F, D, Omega.
 * We stress that the first argument is GMST with a 180 degree phase shift
 * \param tdbTime Time in TDB at which the arguments are to be calculated, in seconds since J2000.
 * \param terrestrialTime Time in TT at which the arguments are to be calculated, in seconds since J2000.
 * \param universalTime1 Time in UT1 at which the arguments are to be calculated, in seconds since J2000.
 * \return Delaunay fundamental arguments and (GMST + pi) at the requested time.
 */
Eigen::Vector6d  calculateDelaunayFundamentalArgumentsWithGmst(
        const double tdbTime, const double terrestrialTime, const double universalTime1 );

//! Function to calculate the Delaunay fundamental arguments and (GMST + pi) at the requested time.
/*!
 * Function to calculate the Delaunay fundamental arguments and (GMST + pi) at the requested time, using Sofa implementation
 * of angle calculations. When using this function, TT and TDB are assumed equal, as are UTC and UT1, for the purposes of
 * calculating GMST (which requires TT and UT1).
 * Using the typical notation of teh angles, the output order is: (GMST + pi), l, l', F, D, Omega.
 * We stress that the first argument is GMST with a 180 degree phase shift
 * \param tdbTime Time in TDB at which the arguments are to be calculated, in seconds since J2000.
 * \return Delaunay fundamental arguments and (GMST + pi) at the requested time.
 */
Eigen::Vector6d  calculateDelaunayFundamentalArgumentsWithGmst(
        const double tdbTime );

//! Function to calculate the Doodson arguments at the requested time.
/*!
 * Function to calculate the Doodson arguments at the requested time, using Sofa implementation
 * of angle calculations.
 * Using the typical notation of teh angles, the output order is: (GMST + pi - s), s, h, p, N', p_{s}.
 * \param tdbTime Time in TDB at which the arguments are to be calculated, in seconds since J2000.
 * \param terrestrialTime Time in TT at which the arguments are to be calculated, in seconds since J2000.
 * \param universalTime1 Time in UT1 at which the arguments are to be calculated, in seconds since J2000.
 * \return Doodson arguments and at the requested time.
 */
Eigen::Vector6d calculateDoodsonFundamentalArguments(
        const double tdbTime, const double terrestrialTime, const double universalTime1 );

//! Function to calculate the Doodson arguments at the requested time.
/*!
 * Function to calculate the Doodson arguments at the requested time, using Sofa implementation
 * of angle calculations. he Delaunay fundamental arguments and (GMST + pi) at the requested time, using Sofa implementation
 * of angle calculations. When using this function, TT and TDB are assumed equal, as are UTC and UT1, for the purposes of
 * calculating GMST (which requires TT and UT1).
 * Using the typical notation of teh angles, the output order is: (GMST + pi - s), s, h, p, N', p_{s}.
 * \param tdbTime Time in TDB at which the arguments are to be calculated, in seconds since J2000.
 * \return Doodson arguments and at the requested time.
 */
Eigen::Vector6d calculateDoodsonFundamentalArguments( const double tdbTime );

} // namespace sofa_interfaces

} // namespace tudat

#endif // TUDAT_SOFAFUNDAMENTALARGUMENTS_H
