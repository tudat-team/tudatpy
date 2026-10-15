#include <array>
#include <exception>
#include <set>
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/numericalDerivative.h"
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
std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors( const Eigen::Matrix3d& toPlanetFixedFrameMatrix )
{
    std::vector< Eigen::Vector3d > geocentricUnitVectors;
    geocentricUnitVectors.resize( 3 );
    geocentricUnitVectors[ 0 ] = toPlanetFixedFrameMatrix.block( 0, 0, 3, 1 );
    geocentricUnitVectors[ 1 ] = toPlanetFixedFrameMatrix.block( 0, 1, 3, 1 );
    geocentricUnitVectors[ 2 ] = toPlanetFixedFrameMatrix.block( 0, 2, 3, 1 );
    return geocentricUnitVectors;
}

//! Function to generate unit vectors of topocentric frame.
std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors( const double latitude, const double longitude )
{
    return getGeocentricLocalUnitVectors( Eigen::Matrix3d(
            reference_frames::getEnuLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion( longitude, latitude ) ) );
}

//! Constructor
GroundStationState::GroundStationState( const Eigen::Vector3d& stationPosition,
                                        const coordinate_conversions::PositionElementTypes inputElementType,
                                        const std::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface,
                                        const std::shared_ptr< StationMotionModel > stationMotionModel ):
    bodyShapeModel_( bodySurface ), stationMotionModel_( stationMotionModel )
{
    resetGroundStationPositionAtEpoch( stationPosition, inputElementType );
}

//! Function to obtain the Cartesian state of the ground station in the local frame at a given time.
Eigen::Vector6d GroundStationState::getCartesianStateInTime( const Time& secondsSinceEpoch, const std::string& targetFrameOrigin )
{
    Eigen::Vector6d currentStationState = ( stationMotionModel_ != nullptr )
            ? ( nominalCartesianState_ +
                stationMotionModel_->getBodyFixedStationMotion( secondsSinceEpoch, shared_from_this( ), targetFrameOrigin ) )
            : nominalCartesianState_;

    return currentStationState;
}

Eigen::Vector6d GroundStationState::getCartesianStateInTime( const double secondsSinceEpoch, const std::string& targetFrameOrigin )
{
    return getCartesianStateInTime( Time( secondsSinceEpoch ), targetFrameOrigin );
}

//! Function to (re)set the nominal state of the station
void GroundStationState::resetGroundStationPositionAtEpoch( const Eigen::Vector3d stationPosition,
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
            bodyShapeModel_, getNominalLatitude( ), getNominalLongitude( ), cartesianPosition_ );
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
        topocentricUnitVectors = getGeocentricLocalUnitVectors( geocentricLatitude, geocentricLongitude );
    }
    else if( std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( bodyShapeModel ) != nullptr )
    {
        std::shared_ptr< basic_astrodynamics::OblateSpheroidBodyShapeModel > oblateSphericalShapeModel =
                std::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( bodyShapeModel );

        // Calculate geodetic latitude.
        double flattening = oblateSphericalShapeModel->getFlattening( );
        double equatorialRadius = oblateSphericalShapeModel->getEquatorialRadius( );
        double geodeticLatitude = coordinate_conversions::calculateGeodeticLatitude( localPoint, equatorialRadius, flattening, 1.0E-4 );

        // Calculte unit vectors of topocentric frame.
        topocentricUnitVectors = getGeocentricLocalUnitVectors( geodeticLatitude, geocentricLongitude );
    }
    else
    {
        // Assume spherical shape
        topocentricUnitVectors = getGeocentricLocalUnitVectors( geocentricLatitude, geocentricLongitude );
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
        const Time& time,
        const std::shared_ptr< ground_stations::GroundStationState > groundStationState,
        const std::string& targetFrameOrigin )
{
    Eigen::Vector6d stationMotion = Eigen::Vector6d::Zero( );
    if( !( time < Time( firstDisplacementTime_ ) ) )
    {
        if( time >= Time( finalDisplacementTime_ ) )
        {
            stationMotion.segment( 0, 3 ) += displacementVectors_.at( displacementVectors_.size( ) - 1 );
        }
        else
        {
            stationMotion.segment( 0, 3 ) +=
                    displacementVectors_.at( std::upper_bound( displacementTimes_.begin( ),
                                                               displacementTimes_.end( ),
                                                               time,
                                                               []( const Time& query, double node ) { return query < Time( node ); } ) -
                                             displacementTimes_.begin( ) - 1 );
        }
    }
    return stationMotion;
}

namespace
{
constexpr double iersTtToTdbSpatialScale = 1.48082686741e-8;
}  // namespace

BodyCentricToBarycentricRelativisticStationMotion::BodyCentricToBarycentricRelativisticStationMotion(
        const std::function< Eigen::Vector6d( const Time& ) >& bodyBarycentricStateFunction,
        const std::function< Eigen::Quaterniond( const Time& ) >& inertialToBodyFixedRotationFunction,
        const std::function< Eigen::Matrix3d( const Time& ) >& rotationDerivativeFunction,
        const std::vector< RelativisticStationExternalBody >& externalBodies,
        const RelativisticStationCoordinateCompatibility coordinateCompatibility,
        const double derivativeStep ):
    bodyBarycentricStateFunction_( bodyBarycentricStateFunction ),
    inertialToBodyFixedRotationFunction_( inertialToBodyFixedRotationFunction ), rotationDerivativeFunction_( rotationDerivativeFunction ),
    externalBodies_( externalBodies ), coordinateCompatibility_( coordinateCompatibility ), derivativeStep_( derivativeStep )
{
    if( !bodyBarycentricStateFunction_ || !inertialToBodyFixedRotationFunction_ )
        throw std::invalid_argument( "Relativistic station transformation requires host state and rotation callbacks." );
    if( coordinateCompatibility != unscaled_tcg_to_tcb_station_coordinates && coordinateCompatibility != tt_to_tdb_station_coordinates )
        throw std::invalid_argument( "Select explicit coordinate compatibility for the IERS spatial station approximation." );
    if( !std::isfinite( derivativeStep ) || derivativeStep <= 0 )
        throw std::invalid_argument( "Relativistic station derivative step must be finite and positive." );
    std::set< std::string > names;
    for( const auto& body : externalBodies_ )
        if( body.name.empty( ) || !names.insert( body.name ).second || !body.barycentricState || !body.gravitationalParameter )
            throw std::invalid_argument( "Relativistic external bodies require unique names, state and GM callbacks." );
}

double BodyCentricToBarycentricRelativisticStationMotion::potential( const Time& time )
{
    if( externalBodies_.empty( ) ) return 0;
    const Eigen::Vector6d host = bodyBarycentricStateFunction_( time );
    if( !host.head< 3 >( ).allFinite( ) ) throw std::runtime_error( "Non-finite host position in relativistic station potential." );
    double value = 0;
    const auto contribution = [ & ]( const Eigen::Vector3d& position, const double gm ) {
        const double distance = ( host.head< 3 >( ) - position ).norm( );
        if( !position.allFinite( ) || !std::isfinite( gm ) || gm < 0 || !std::isfinite( distance ) || distance <= 0 )
            throw std::domain_error( "Relativistic station external potential requires positive separation and finite nonnegative GM." );
        return gm / distance;
    };
    for( const auto& body : externalBodies_ )
        value += contribution( body.barycentricState( time ).head< 3 >( ), body.gravitationalParameter( ) );
    if( !std::isfinite( value ) ) throw std::runtime_error( "Non-finite relativistic station potential." );
    return value;
}

double BodyCentricToBarycentricRelativisticStationMotion::potentialDerivative( const Time& time )
{
    if( externalBodies_.empty( ) ) return 0;
    const Eigen::Vector6d host = bodyBarycentricStateFunction_( time );
    if( !host.allFinite( ) ) throw std::runtime_error( "Non-finite host state in relativistic station potential rate." );
    double rate = 0;
    for( const auto& body : externalBodies_ )
    {
        const Eigen::Vector6d external = body.barycentricState( time );
        if( !external.allFinite( ) ) throw std::runtime_error( "External potential velocity is unavailable or non-finite." );
        const Eigen::Vector3d separation = host.head< 3 >( ) - external.head< 3 >( );
        const double distance = separation.norm( ), gm = body.gravitationalParameter( );
        if( !std::isfinite( distance ) || distance <= 0 || !std::isfinite( gm ) || gm < 0 )
            throw std::domain_error( "Invalid external potential separation or GM." );
        rate -= gm * separation.dot( host.tail< 3 >( ) - external.tail< 3 >( ) ) / ( distance * distance * distance );
    }
    if( !std::isfinite( rate ) ) throw std::runtime_error( "Non-finite external potential rate." );
    return rate;
}

Eigen::Matrix3d BodyCentricToBarycentricRelativisticStationMotion::correctionMatrix( const Time& time )
{
    const Eigen::Vector6d host = bodyBarycentricStateFunction_( time );
    const Eigen::Quaterniond rotation = inertialToBodyFixedRotationFunction_( time );
    if( !host.tail< 3 >( ).allFinite( ) || !rotation.coeffs( ).allFinite( ) || std::fabs( rotation.norm( ) - 1 ) > 1e-12 )
        throw std::runtime_error( "Invalid host velocity or rotation in relativistic station position." );
    const Eigen::Vector3d localVelocity = rotation * host.tail< 3 >( );
    const double scale = coordinateCompatibility_ == tt_to_tdb_station_coordinates ? iersTtToTdbSpatialScale : 0;
    return -( potential( time ) * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT + scale ) * Eigen::Matrix3d::Identity( ) -
            0.5 * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * localVelocity * localVelocity.transpose( );
}

Eigen::Matrix3d BodyCentricToBarycentricRelativisticStationMotion::correctionMatrixDerivative( const Time& time )
{
    const Eigen::Matrix3d rotation = inertialToBodyFixedRotationFunction_( time ).toRotationMatrix( );

    // Rotation rate from its callback if available, otherwise numerically; the host acceleration is always numerical.
    const std::function< Eigen::Matrix3d( const Time& ) > rotationMatrix = [ this ]( const Time& t ) -> Eigen::Matrix3d {
        return inertialToBodyFixedRotationFunction_( t ).toRotationMatrix( );
    };
    const std::function< Eigen::Vector3d( const Time& ) > hostVelocity = [ this ]( const Time& t ) -> Eigen::Vector3d {
        return bodyBarycentricStateFunction_( t ).tail< 3 >( );
    };
    const Eigen::Matrix3d rotationRate = rotationDerivativeFunction_
            ? rotationDerivativeFunction_( time )
            : numerical_derivatives::computeFourthOrderCentralDifference( rotationMatrix, time, derivativeStep_ );
    const Eigen::Vector3d acceleration = numerical_derivatives::computeFourthOrderCentralDifference( hostVelocity, time, derivativeStep_ );
    const Eigen::Vector3d velocity = bodyBarycentricStateFunction_( time ).tail< 3 >( );
    const Eigen::Vector3d localVelocity = rotation * velocity, localRate = rotationRate * velocity + rotation * acceleration;
    const Eigen::Matrix3d rate =
            -potentialDerivative( time ) * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * Eigen::Matrix3d::Identity( ) -
            0.5 * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                    ( localRate * localVelocity.transpose( ) + localVelocity * localRate.transpose( ) );
    if( !rate.allFinite( ) ) throw std::runtime_error( "Relativistic station spatial-map rate is unavailable or non-finite." );
    return rate;
}

Eigen::Vector3d BodyCentricToBarycentricRelativisticStationMotion::correctCompleteLocalPosition( const Time& time,
                                                                                                 const Eigen::Vector3d& position,
                                                                                                 const std::string& origin )
{
    if( origin != "SSB" ) return Eigen::Vector3d::Zero( );
    if( !position.allFinite( ) ) throw std::runtime_error( "Non-finite intrinsic station position." );
    return correctionMatrix( time ) * position;
}

Eigen::Vector6d BodyCentricToBarycentricRelativisticStationMotion::correctCompleteLocalState( const Time& time,
                                                                                              const Eigen::Vector6d& state,
                                                                                              const std::string& origin )
{
    if( origin != "SSB" ) return Eigen::Vector6d::Zero( );
    if( !state.allFinite( ) ) throw std::runtime_error( "Intrinsic station state is unavailable or non-finite." );
    const Eigen::Matrix3d correction = correctionMatrix( time ), rate = correctionMatrixDerivative( time );
    return ( Eigen::Vector6d( ) << correction * state.head< 3 >( ), rate * state.head< 3 >( ) + correction * state.tail< 3 >( ) )
            .finished( );
}

Eigen::Vector6d BodyCentricToBarycentricRelativisticStationMotion::getBodyFixedStationMotion(
        const Time& time,
        const std::shared_ptr< GroundStationState > station,
        const std::string& origin )
{
    if( !station ) throw std::invalid_argument( "Relativistic station correction requires a station." );
    Eigen::Vector6d nominal = Eigen::Vector6d::Zero( );
    nominal.head< 3 >( ) = station->getNominalCartesianPosition( );
    return correctCompleteLocalState( time, nominal, origin );
}

Eigen::Vector3d BodyCentricToBarycentricRelativisticStationMotion::getBodyFixedStationPositionCorrection(
        const Time& time,
        const std::shared_ptr< GroundStationState > station,
        const std::string& origin )
{
    if( !station ) throw std::invalid_argument( "Relativistic station correction requires a station." );
    return correctCompleteLocalPosition( time, station->getNominalCartesianPosition( ), origin );
}

CombinedStationMotionModel::CombinedStationMotionModel( const std::vector< std::shared_ptr< StationMotionModel > >& models )
{
    for( const auto& model : models ) appendModel( model );
}

void CombinedStationMotionModel::appendModel( const std::shared_ptr< StationMotionModel >& model )
{
    if( !model ) throw std::invalid_argument( "Null station motion model." );
    if( const auto combined = std::dynamic_pointer_cast< CombinedStationMotionModel >( model ) )
    {
        for( const auto& intrinsic : combined->modelList_ ) appendModel( intrinsic );
        if( combined->terminalModel_ ) appendModel( combined->terminalModel_ );
    }
    else if( const auto terminal = std::dynamic_pointer_cast< BodyCentricToBarycentricRelativisticStationMotion >( model ) )
    {
        if( terminalModel_ ) throw std::invalid_argument( "Only one terminal relativistic station transformation is permitted." );
        terminalModel_ = terminal;
    }
    else
        modelList_.push_back( model );
}

//! Sum intrinsic motion, then apply the terminal relativistic map to the complete intrinsic state.
Eigen::Vector6d CombinedStationMotionModel::getBodyFixedStationMotion( const Time& time,
                                                                       const std::shared_ptr< GroundStationState > station,
                                                                       const std::string& origin )
{
    Eigen::Vector6d motion = Eigen::Vector6d::Zero( );
    for( const auto& model : modelList_ )
    {
        motion += model->getBodyFixedStationMotion( time, station, origin );
    }
    if( terminalModel_ )
    {
        Eigen::Vector6d intrinsicState = motion;
        intrinsicState.head< 3 >( ) += station->getNominalCartesianPosition( );
        motion += terminalModel_->correctCompleteLocalState( time, intrinsicState, origin );
    }
    return motion;
}

//! Sum intrinsic position corrections, then apply the terminal relativistic map to the complete intrinsic position.
Eigen::Vector3d CombinedStationMotionModel::getBodyFixedStationPositionCorrection( const Time& time,
                                                                                   const std::shared_ptr< GroundStationState > station,
                                                                                   const std::string& origin )
{
    Eigen::Vector3d displacement = Eigen::Vector3d::Zero( );
    for( const auto& model : modelList_ )
    {
        displacement += model->getBodyFixedStationPositionCorrection( time, station, origin );
    }
    if( terminalModel_ )
    {
        displacement +=
                terminalModel_->correctCompleteLocalPosition( time, station->getNominalCartesianPosition( ) + displacement, origin );
    }
    return displacement;
}

}  // namespace ground_stations

}  // namespace tudat
