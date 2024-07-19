#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/mathematicalConstants.h"

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

namespace tudat
{

namespace ground_stations
{

//! Function to generate unit vectors of topocentric frame.
std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors(
        const Eigen::Matrix3d& toPlanetFixedFrameMatrix )
{
    std::vector< Eigen::Vector3d > geocentricUnitVectors;
    geocentricUnitVectors.resize( 3 );
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
        const Eigen::Vector3d& stationPosition,
        const coordinate_conversions::PositionElementTypes inputElementType,
        const std::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface,
        const std::shared_ptr< StationMotionModel > stationMotionModel ):
    bodyShapeModel_( bodySurface ),
    stationMotionModel_( stationMotionModel )
{
    resetGroundStationPositionAtEpoch( stationPosition, inputElementType );
}

//! Function to obtain the Cartesian state of the ground station in the local frame at a given time.
Eigen::Vector6d GroundStationState::getCartesianStateInTime(
        const double secondsSinceEpoch, const std::string& targetFrameOrigin )
{
    Eigen::Vector6d currentStationState =
            ( stationMotionModel_ != nullptr ) ? ( nominalCartesianState_ + stationMotionModel_->getBodyFixedStationMotion(
                                                       secondsSinceEpoch, shared_from_this( ), targetFrameOrigin ) ) : nominalCartesianState_;

    return currentStationState;
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
                stationPosition, inputElementType, coordinate_conversions::cartesian_position, bodyShapeModel_ );
    nominalCartesianState_.segment( 0, 3 ) = cartesianPosition_;
    nominalCartesianState_.segment( 3, 3 ).setZero( );

    sphericalPosition_ = coordinate_conversions::convertPositionElements(
                stationPosition, inputElementType, coordinate_conversions::spherical_position, bodyShapeModel_ );

    // If possible, set geodetic position, otherwise, set to NaN.
    try
    {
        geodeticPosition = coordinate_conversions::convertPositionElements(
                    stationPosition, inputElementType, coordinate_conversions::geodetic_position, bodyShapeModel_ );
    }
    catch( std::runtime_error const& )
    {
        geodeticPosition = Eigen::Vector3d::Constant( TUDAT_NAN );
    }

    setTransformationAndUnitVectors( );

}

//! Function to reset the rotation from the body-fixed to local topocentric frame, and associated unit vectors
void GroundStationState::setTransformationAndUnitVectors( )
{
    geocentricUnitVectors_ = getGeocentricLocalUnitVectors( getNominalLatitude( ), getNominalLongitude( ) );
    bodyFixedToTopocentricFrameRotation_ = getRotationQuaternionFromBodyFixedToTopocentricFrame(
                bodyShapeModel_, getNominalLatitude( ), getNominalLongitude( ), cartesianPosition_  );
}

//! Function to calculate the rotation from a body-fixed to a topocentric frame.
Eigen::Quaterniond getRotationQuaternionFromBodyFixedToTopocentricFrame(
        const std::shared_ptr< basic_astrodynamics::BodyShapeModel > bodyShapeModel,
        const double geocentricLatitude,
        const double geocentricLongitude,
        const Eigen::Vector3d localPoint )
{
    // Declare unit vectors of topocentric frame, to be calculated.
    std::vector< Eigen::Vector3d > topocentricUnitVectors;

    bool isSurfaceModelRecognized = 1;

    // Identify type of body shape model
    if( std::dynamic_pointer_cast< basic_astrodynamics::SphericalBodyShapeModel >( bodyShapeModel ) != nullptr )
    {
        // For a sphere the topocentric and geocentric frames are equal.
        topocentricUnitVectors = getGeocentricLocalUnitVectors(
                    geocentricLatitude, geocentricLongitude );
    }
    else if( std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( bodyShapeModel ) != nullptr )
    {
        std::shared_ptr< basic_astrodynamics::OblateSpheroidBodyShapeModel > oblateSphericalShapeModel =
                std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( bodyShapeModel );

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
        // Assume spherical shape
        topocentricUnitVectors = getGeocentricLocalUnitVectors(
            geocentricLatitude, geocentricLongitude );
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

Eigen::Vector6d PiecewiseConstantStationMotionModel::getBodyFixedStationMotion(
        const double time,
        const std::shared_ptr< ground_stations::GroundStationState > groundStationState,
        const std::string& targetFrameOrigin )
{
    timeLookupScheme_ = std::make_shared< interpolators::BinarySearchLookupScheme< double > >(
                displacementTimes_ );
    Eigen::Vector6d stationMotion = Eigen::Vector6d::Zero( );
    if( !( time < firstDisplacementTime_ ) )
    {
        if( time >= finalDisplacementTime_ )
        {
            stationMotion.segment( 0, 3 ) += displacementVectors_.at(
                        displacementVectors_.size( ) - 1 );
        }
        else
        {
            stationMotion.segment( 0, 3 ) += displacementVectors_.at(
                        timeLookupScheme_->findNearestLowerNeighbour( time ) );
        }
    }
    return stationMotion;
}

Eigen::Vector6d BodyCentricToBarycentricRelativisticStationMotion::getBodyFixedStationMotion(
    const double time,
    const std::shared_ptr< ground_stations::GroundStationState > groundStationState,
    const std::string& targetFrameOrigin )
{
    Eigen::Vector6d stationMotion = Eigen::Vector6d::Zero( );

    if( targetFrameOrigin == "SSB" )
    {
        currentRotationToBodyFixedFrame_ = inertialToBodyFixedRotationFunction_( time );
        inertialNominalStationPosition_ = currentRotationToBodyFixedFrame_.inverse( ) * groundStationState->getNominalCartesianPosition( );

        centralBodyBarycentricState_ = bodyBarycentricStateFunction_( time );
        stationMotion.segment( 0, 3 ) += ( centralBodyBarycentricState_.segment( 3, 3 ).dot( inertialNominalStationPosition_ ) ) * centralBodyBarycentricState_.segment( 3, 3 );

        if( useGeneralRelativisticCorrection_ )
        {
            stationMotion.segment( 0, 3 ) += ( centralBodyGravitationalParameterFunction_( ) /
                                               ( centralBodyBarycentricState_.segment( 0, 3 ) -
                                                 centralBodyBarycentricPositionFunction_( time )).norm( )) *
                                             inertialNominalStationPosition_;
        }

        stationMotion.segment( 0, 3 ) *= physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
        stationMotion.segment( 0, 3 ) = currentRotationToBodyFixedFrame_ * stationMotion.segment( 0, 3 );
    }
    return stationMotion;
}
}

}

