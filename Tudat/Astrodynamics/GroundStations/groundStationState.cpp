#include <boost/assign/list_of.hpp>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/GroundStations/groundStationState.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

namespace ground_stations
{

//! Function to generate unit vectors of topocentric frame.
std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors(
            const Eigen::Matrix3d& toPlanetFixedFrameMatrix )
{
    std::vector< Eigen::Vector3d > geocentricUnitVectors;
    geocentricUnitVectors.reserve( 3 );
    geocentricUnitVectors[ 0 ] = toPlanetFixedFrameMatrix.block( 0, 0, 3, 1 );
    geocentricUnitVectors[ 1 ] = toPlanetFixedFrameMatrix.block( 0, 1, 3, 1 );
    geocentricUnitVectors[ 2 ] = toPlanetFixedFrameMatrix.block( 0, 2, 3, 1 );
    return geocentricUnitVectors;
}


//! Function to generate unit vectors of topocentric frame.
std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors(
        const double latitude, const double longitude )
{
    return getGeocentricLocalUnitVectors(
                Eigen::Matrix3d( reference_frames::getEnuLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
                                 longitude, latitude ) ) );
}

//! Constructor
GroundStationState::GroundStationState(
        const Eigen::Vector3d stationPosition,
        const coordinate_conversions::PositionElementTypes inputElementType,
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface ):
    bodySurface_( bodySurface )
{
    resetGroundStationPositionAtEpoch( stationPosition, inputElementType );
}

//! Function to obtain the Cartesian state of the ground station in the local frame at a given time.
Eigen::Vector6d GroundStationState::getCartesianStateInTime(
        const double secondsSinceEpoch,
        const double inputReferenceEpoch )
{
    return ( Eigen::Vector6d( ) << cartesianPosition_, Eigen::Vector3d::Zero( ) ).finished( );
}

//! Function to (re)set the nominal state of the station
void GroundStationState::resetGroundStationPositionAtEpoch(
                const Eigen::Vector3d stationPosition,
                const coordinate_conversions::PositionElementTypes inputElementType )
{
    using namespace coordinate_conversions;
    using mathematical_constants::PI;

    // Set Cartesian and spherical position
    cartesianPosition_ = coordinate_conversions::convertPositionElements(
                stationPosition, inputElementType, coordinate_conversions::cartesian_position, bodySurface_ );
    sphericalPosition_ = coordinate_conversions::convertPositionElements(
                stationPosition, inputElementType, coordinate_conversions::spherical_position, bodySurface_ );

    // If possible, set geodetic position, otherwise, set to NaN.
    try
    {
        geodeticPosition = coordinate_conversions::convertPositionElements(
                    stationPosition, inputElementType, coordinate_conversions::geodetic_position, bodySurface_ );
    }
    catch( std::runtime_error )
    {
        geodeticPosition = Eigen::Vector3d::Constant( TUDAT_NAN );
    }

    setTransformationAndUnitVectors( );

}

//! Function to reset the rotation from the body-fixed to local topocentric frame, and associated unit vectors
void GroundStationState::setTransformationAndUnitVectors( )
{
    geocentricUnitVectors_ = getGeocentricLocalUnitVectors( getLongitude( ), getLatitude( ) );
    bodyFixedToTopocentricFrameRotation_ = getRotationQuaternionFromBodyFixedToTopocentricFrame(
                bodySurface_, getLatitude( ), getLongitude( ), cartesianPosition_  );
}

//! Function to calculate the rotation from a body-fixed to a topocentric frame.
Eigen::Quaterniond getRotationQuaternionFromBodyFixedToTopocentricFrame(
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodyShapeModel,
        const double geocentricLatitude,
        const double geocentricLongitude,
        const Eigen::Vector3d localPoint )
{
    // Declare unit vectors of topocentric frame, to be calculated.
    std::vector< Eigen::Vector3d > topocentricUnitVectors;

    bool isSurfaceModelRecognized = 1;

    // Identify type of body shape model
    if( boost::dynamic_pointer_cast< basic_astrodynamics::SphericalBodyShapeModel >( bodyShapeModel ) != NULL )
    {
        // For a sphere the topocentric and geocentric frames are equal.
        topocentricUnitVectors = getGeocentricLocalUnitVectors(
                    geocentricLatitude, geocentricLongitude );
    }
    else if( boost::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( bodyShapeModel ) != NULL )
    {
        boost::shared_ptr< basic_astrodynamics::OblateSpheroidBodyShapeModel > oblateSphericalShapeModel =
                boost::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( bodyShapeModel );

        // Calculate geodetic latitude.
        double flattening = oblateSphericalShapeModel->getFlattening( );
        double equatorialRadius = oblateSphericalShapeModel->getEquatorialRadius( );
        double geodeticLatitude = coordinate_conversions::calculateGeodeticLatitude(
                    localPoint, equatorialRadius, flattening, 1.0E-4 );

        // Calculte unit vectors of topocentric frame.
        topocentricUnitVectors = getGeocentricLocalUnitVectors( geodeticLatitude, geocentricLongitude );
    }
    else
    {
        isSurfaceModelRecognized = 0;
        throw std::runtime_error( "Error when making transformation to topocentric frame, shape model not recognized" );
    }

    // Create rotation matrix

    Eigen::Matrix3d bodyFixedToTopocentricFrame;

    if( isSurfaceModelRecognized == 1 )
    {
        bodyFixedToTopocentricFrame.block( 0, 0, 1, 3 ) = topocentricUnitVectors[ 0 ].transpose( );
        bodyFixedToTopocentricFrame.block( 1, 0, 1, 3 ) = topocentricUnitVectors[ 1 ].transpose( );
        bodyFixedToTopocentricFrame.block( 2, 0, 1, 3 ) = topocentricUnitVectors[ 2 ].transpose( );
    }
    else
    {
        bodyFixedToTopocentricFrame = Eigen::Matrix3d::Identity( );
    }



    // Convert to quaternion and return.
    return Eigen::Quaterniond( bodyFixedToTopocentricFrame );
}

}


}
