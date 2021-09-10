/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/ground_stations/pointingAnglesCalculator.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{

namespace ground_stations
{

//! Function to calculate the elevation angle from body-fixed point to given point.
double PointingAnglesCalculator::calculateElevationAngle(
        const Eigen::Vector3d inertialVectorAwayFromStation,
        const double time )
{
    // Transform vector to local topocentric frame.
    Eigen::Vector3d vectorInTopoCentricFrame = convertVectorFromInertialToTopocentricFrame(
                inertialVectorAwayFromStation, time );

    // Calculate and return elevation angle.
    return mathematical_constants::PI / 2.0 - linear_algebra::computeAngleBetweenVectors(
                 vectorInTopoCentricFrame, Eigen::Vector3d::UnitZ( ) );
}

//! Function to calculate the azimuth angle from body-fixed point to given point.
double PointingAnglesCalculator::calculateAzimuthAngle( const Eigen::Vector3d inertialVectorAwayFromStation,
                                const double time )
{
    // Transform vector to local topocentric frame.
    Eigen::Vector3d vectorInTopoCentricFrame = convertVectorFromInertialToTopocentricFrame(
                inertialVectorAwayFromStation, time );

    // Calculate and return azimuth angle.
    return std::atan2( vectorInTopoCentricFrame.y( ), vectorInTopoCentricFrame.x( ) );
}

//! Function to calculate the elevation and azimuth angles from body-fixed point to given point.
std::pair< double, double > PointingAnglesCalculator::calculatePointingAngles(
        const Eigen::Vector3d inertialVectorAwayFromStation, const double time )
{
    // Transform vector to local topocentric frame.
    Eigen::Vector3d vectorInTopoCentricFrame = convertVectorFromInertialToTopocentricFrame(
                inertialVectorAwayFromStation, time );

    // Calculate elevation angle.
    double elevationAngle = mathematical_constants::PI / 2.0 - linear_algebra::computeAngleBetweenVectors(
                vectorInTopoCentricFrame, Eigen::Vector3d::UnitZ( ) );

    // Calculate azimuth angle.
    double azimuthAngle = std::atan2( vectorInTopoCentricFrame.y( ), vectorInTopoCentricFrame.x( ) );

    // Return angles.
    return std::make_pair( elevationAngle, azimuthAngle );
}

//! Function to convert vector in inertial frame to topocentric frame.
Eigen::Vector3d PointingAnglesCalculator::convertVectorFromInertialToTopocentricFrame(
        const Eigen::Vector3d& inertialVector, const double time )
{
    // Calculate anf combine constituent rotations.
    return rotationFromBodyFixedToTopoCentricFrame_( time ) * rotationFromInertialToBodyFixedFrame_( time ) * inertialVector;
}

double calculateGroundStationElevationAngle(
        const std::shared_ptr< PointingAnglesCalculator > angleCalculator,
        const std::vector< Eigen::Vector6d > linkEndStates,
        const std::vector< double > linkEndTimes,
        const std::pair< int, int >& linkEndIndices )
{
    double stationTime = linkEndTimes.at( linkEndIndices.first );
    Eigen::Vector3d targetRelativeState = ( linkEndStates.at( linkEndIndices.second ) -
            linkEndStates.at( linkEndIndices.first ) ).segment( 0, 3 );
    return angleCalculator->calculateElevationAngle( targetRelativeState, stationTime );
}

double calculateGroundStationAzimuthAngle(
        const std::shared_ptr< PointingAnglesCalculator > angleCalculator,
        const std::vector< Eigen::Vector6d > linkEndStates,
        const std::vector< double > linkEndTimes,
        const std::pair< int, int >& linkEndIndices )
{
    double stationTime = linkEndTimes.at( linkEndIndices.first );
    Eigen::Vector3d targetRelativeState = ( linkEndStates.at( linkEndIndices.second ) -
            linkEndStates.at( linkEndIndices.first ) ).segment( 0, 3 );
    return angleCalculator->calculateAzimuthAngle( targetRelativeState, stationTime );
}


}

}
