/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GROUNDSTATIONSTATE_H
#define TUDAT_GROUNDSTATIONSTATE_H

#include <memory>

#include "tudat/astro/basic_astro/bodyShapeModel.h"
#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"
#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace ground_stations
{

//! Function to generate unit vectors of topocentric frame.
/*!
 *  Function to generate unit vectors of topocentric frame. The set of unit vectors of topocentric frame is expressed in
 *  body-fixed frame, in ENU (Earth-North-Up) order.
 *  \param topocentricToPlanetFixedFrameMatrix Rotation matrix from topocentric to body-fixed frame.
 *  \return Unit vectors of topocentric frame, expressed in body-fixed frame.
 */
std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors( const Eigen::Matrix3d& topocentricToPlanetFixedFrameMatrix );

//! Function to generate unit vectors of topocentric frame.
/*!
 *  Function to generate unit vectors of topocentric frame. The set of unit vectors of topocentric frame is expressed in
 *  body-fixed frame, in ENU (Earth-North-Up) order.
 *  \param latitude Latitude of point for which topocentric frame is to be created.
 *  \param longitude Longitude of point for which  topocentric frame is to be created.
 *  \return Unit vectors of topocentric frame, expressed in body-fixed frame
 */
std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors( const double latitude, const double longitude );

class GroundStationState;

struct StationMotionModel {
public:
    StationMotionModel( ) {}

    virtual ~StationMotionModel( ) {}

    //! Function to compute the change in the body-fixed position and velocity of the station due to this type of motion.
    /*!
     *  Function to compute the change in the body-fixed position (in metres) and velocity (in metres per second) of the
     *  station due to this type of motion (for instance plate motion, or a type of tidal deformation). Each derived class
     *  implements this function.
     *  \param tdb TDB seconds since J2000 at which the motion is computed
     *  \param station Ground station, from which the model can retrieve the nominal position of the station
     *  \param origin Name of the origin of the frame in which the station state is used; only the relativistic correction
     *  (BodyCentricToBarycentricRelativisticStationMotion) depends on it
     *  \return Change in body-fixed position (first three entries) and velocity (last three entries)
     */
    virtual Eigen::Vector6d getBodyFixedStationMotion( const Time& tdb,
                                                       const std::shared_ptr< GroundStationState > station,
                                                       const std::string& origin ) = 0;

    //! Function to compute the change in the body-fixed position of the station due to this type of motion.
    /*!
     *  Function to compute the change in the body-fixed position (in metres) of the station due to this type of motion. By
     *  default, this is the position part of getBodyFixedStationMotion. A derived class overrides this function when it can
     *  compute the position in cases where it cannot compute the velocity (for instance at a data point next to a gap in the
     *  data), or when the position alone is much faster to compute.
     *  \param tdb TDB seconds since J2000 at which the motion is computed
     *  \param station Ground station, from which the model can retrieve the nominal position of the station
     *  \param origin Name of the origin of the frame in which the station position is used (see getBodyFixedStationMotion)
     *  \return Change in body-fixed position
     */
    virtual Eigen::Vector3d getBodyFixedStationPositionCorrection( const Time& tdb,
                                                                   const std::shared_ptr< GroundStationState > station,
                                                                   const std::string& origin )
    {
        return getBodyFixedStationMotion( tdb, station, origin ).head< 3 >( );
    }
};

//! Class storing and computing the (time-variable) state of a ground station in a body-fixed frame.
class GroundStationState : public std::enable_shared_from_this< GroundStationState >
{
public:
    //! Constructor
    /*!
     *  Constructor
     * \param stationPosition Position of ground station in body-fixed frame
     * \param inputElementType Element type (e.g. Cartesian, spherical, etc.) of groundStationPosition.
     * \param bodySurface Shape of body on which state is defined. If nullptr (default), no conversions to/from
     * geodetic position are possible.
     */
    explicit GroundStationState(
            const Eigen::Vector3d& stationPosition,
            const coordinate_conversions::PositionElementTypes inputElementType = coordinate_conversions::cartesian_position,
            const std::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface = nullptr,
            const std::shared_ptr< StationMotionModel > stationMotionModel = nullptr );

    virtual ~GroundStationState( ) {}

    //! Function to obtain the Cartesian state of the ground station in the local frame at a given TDB epoch.
    /*!
     *  Function to obtain the Cartesian state of the ground station in the local frame (body-fixed, not topocentric) at a
     *  given time.  Adds all position variations to the nominal state (at the requested time) and returns the state. NOTE:
     *  poisition variations are as yet not included.
     *  \param secondsSinceEpoch TDB seconds since J2000 at which the state or position is retrieved.
     *  \param inputReferenceEpoch Reference epoch julian day
     *  \return Cartesian state of station in local frame at requested time.
     */
    Eigen::Vector6d getCartesianStateInTime( const Time& tdb, const std::string& targetFrameOrigin );
    Eigen::Vector6d getCartesianStateInTime( const double secondsSinceEpoch, const std::string& targetFrameOrigin );

    //! Function to obtain the Cartesian position of the ground station in the local frame at a given TDB epoch.
    /*!
     *  Function to obtain the Cartesian position of the ground station in the local frame (body-fixed, not topocentric) at a
     *  given time.  Adds all position variations to the nominal state (at the requested time) and returns the state. NOTE:
     *  poisition variations are as yet not included.
     *  \param secondsSinceEpoch TDB seconds since J2000 at which the state or position is retrieved.
     *  \param inputReferenceEpoch Reference epoch julian day
     *  \return Cartesian position of station in local frame at requested time.
     */
    Eigen::Vector3d getCartesianPositionInTime( const Time& secondsSinceEpoch, const std::string& targetFrameOrigin )
    {
        Eigen::Vector3d position = cartesianPosition_;
        if( stationMotionModel_ != nullptr )
        {
            position +=
                    stationMotionModel_->getBodyFixedStationPositionCorrection( secondsSinceEpoch, shared_from_this( ), targetFrameOrigin );
        }
        return position;
    }

    Eigen::Vector3d getCartesianPositionInTime( const double secondsSinceEpoch, const std::string& targetFrameOrigin )
    {
        return getCartesianPositionInTime( Time( secondsSinceEpoch ), targetFrameOrigin );
    }

    //! Function to return the nominal (unperturbed) Cartesian position of the station
    /*!
     *  Function to return the nominal Cartesian (unperturbed, i.e. not including linear drift, eccentricity,
     *  tidal variations, etc.) position of the station in body-fixed system.
     *  \return Unperturbed Cartesian position of the station
     */
    Eigen::Vector3d getNominalCartesianPosition( )
    {
        return cartesianPosition_;
    }

    //! Function to return the nominal (unperturbed) spherical position of the station
    /*!
     *  Function to return the nominal spherical (unperturbed, i.e. not including linear drift, eccentricity,
     *  tidal variations, etc.) position of the station in body-fixed system.
     *  \return Unperturbed spherical position of the station
     */
    Eigen::Vector3d getNominalSphericalPosition( )
    {
        return sphericalPosition_;
    }

    //! Function to return the nominal (unperturbed) geodetic position of the station
    /*!
     *  Function to return the nominal geodetic (unperturbed, i.e. not including linear drift, eccentricity,
     *  tidal variations, etc.) position of the station in body-fixed system.
     *  NOTE: This vector is only defined in bodySurface_ was non-nullptr at last setting of nominal state
     *  (call to resetGroundStationPositionAtEpoch).
     *  \return Unperturbed geodetic position of the station
     */
    Eigen::Vector3d getNominalGeodeticPosition( )
    {
        if( !( geodeticPosition( 0 ) == geodeticPosition( 0 ) ) )
        {
            throw std::runtime_error( "Error, retrieving geodetic position from ground station state, but is not defined" );
        }
        return geodeticPosition;
    }

    //! Function to return the geocentric latitude of the station.
    /*!
     *  Function to return the geocentric latitude of the station.
     *  \return Geocentric latitude of the station.
     */
    double getNominalLatitude( )
    {
        return sphericalPosition_.y( );
    }

    //! Function to return the geocentric longtude of the station.
    /*!
     *  Function to return the geocentric longtude of the station.
     *  \return Geocentric longtude of the station.
     */
    double getNominalLongitude( )
    {
        return sphericalPosition_.z( );
    }

    //! Function to return current rotation from body-fixed to topocentric frame
    /*!
     *  Function to return current rotation from body-fixed to topocentric frame, with topocentric xyz-axes in ENU (Earth-North-Up) order.
     *  Note that no relativistic rescaling of vectors is performed, only a rotation is applied.
     *  \param time Time at which rotation quaternion is to be evaluated.
     *  \return Current rotation from body-fixed to topocentrix frame
     */
    Eigen::Quaterniond getRotationFromBodyFixedToTopocentricFrame( const double time )
    {
        return bodyFixedToTopocentricFrameRotation_;
    }

    Eigen::Matrix3d getRotationMatrixFromBodyFixedToTopocentricFrame( const double time )
    {
        return Eigen::Matrix3d( getRotationFromBodyFixedToTopocentricFrame( time ) );
    }

    Eigen::Matrix3d getConstantRotationMatrixFromBodyFixedToTopocentricFrame( )
    {
        return Eigen::Matrix3d( getRotationFromBodyFixedToTopocentricFrame( 0.0 ) );
    }

    //! Function to (re)set the nominal state of the station
    /*!
     *  Function to (re)set the nominal state of the station. Input may be in any type of elements defined in
     *  PositionElementTypes enum.
     * \param stationPosition Position of ground station in body-fixed frame
     * \param inputElementType Element type (e.g. Cartesian, spherical, etc.) of groundStationPosition.
     */
    void resetGroundStationPositionAtEpoch(
            const Eigen::Vector3d stationPosition,
            const coordinate_conversions::PositionElementTypes inputElementType = coordinate_conversions::cartesian_position );

    //! Function to retrieve the shape of body on which state is defined.
    /*!
     *  Function to retrieve the shape of body on which state is defined
     * \ return Shape of body on which state is defined
     */
    std::shared_ptr< basic_astrodynamics::BodyShapeModel > getBodySurface( )
    {
        return bodyShapeModel_;
    }

    std::vector< Eigen::Vector3d > getEnuGeocentricUnitVectors( )
    {
        if( geocentricUnitVectors_.size( ) == 0 )
        {
            try
            {
                setTransformationAndUnitVectors( );
            }
            catch( std::runtime_error& caughtException )
            {
                throw std::runtime_error( "Error when computing ground station unit vectors: " + std::string( caughtException.what( ) ) );
            }
        }
        return geocentricUnitVectors_;
    }

    void setSiteId( const std::string& siteId )
    {
        siteId_ = siteId;
    }

    std::string getSiteId( )
    {
        return siteId_;
    }

protected:
    //! Function to reset the rotation from the body-fixed to local topocentric frame, and associated unit vectors
    void setTransformationAndUnitVectors( );

    //! Cartesian position of station
    /*!
     *  Cartesian position of station, without variations (linear drift, eccentricity, tides, etc.), in the body-fixed frame.
     */
    Eigen::Vector3d cartesianPosition_;

    Eigen::Vector6d nominalCartesianState_;

    //! Spherical position of station
    /*!
     *  Spherical position of station, without variations (linear drift, eccentricity, tides, etc.), in the body-fixed frame.
     *  The order of the entries is: radius, geocentric latitude, longitude.
     */
    Eigen::Vector3d sphericalPosition_;

    //! Geodetic position of station
    /*!
     *  Geodetic position of station, without variations (linear drift, eccentricity, tides, etc.), in the body-fixed frame.
     *  The order of the entries is: altitude, geodetic latitude, longitude. This vector is only defined in bodySurface_ was
     *  non-nullptr at last setting of nominal state (call to resetGroundStationPositionAtEpoch).
     */
    Eigen::Vector3d geodeticPosition;

    //! Set of unit vectors of topocentric frame.
    /*!
     *  Set of unit vectors of topocentric frame, expressed in body-fixed frame, in ENU (Earth-North-Up) order.
     *  The up-vector is generated by drawing a line from the center of the body to the (nominal) position of the station.
     *  For a sphere, this coincides with the true topocentric frame, for which the East-North vectors span a local tangent plane.
     *  The transformation from body-fixed to 'true' topocentric frame is described by the
     *  bodyFixedToTopocentricFrameRotation_ variable.
     */
    std::vector< Eigen::Vector3d > geocentricUnitVectors_;

    //! Rotation from body-fixed to topocentrix frame
    /*!
     *  Rotation from body-fixed to topocentrix frame, with topocentric xyz-axes in ENU (Earth-North-Up) order.
     *  In this frame the East-North vectors span a local tangent plane. The frame is time-independent and is based on the
     *  cartesianPosition_ member, without variations.
     */
    Eigen::Quaterniond bodyFixedToTopocentricFrameRotation_;

    //! Shape of body on which state is defined
    std::shared_ptr< basic_astrodynamics::BodyShapeModel > bodyShapeModel_;

    std::shared_ptr< StationMotionModel > stationMotionModel_;

    std::string siteId_ = "";
};

//! Function to calculate the rotation from a body-fixed to a topocentric frame.
/*!
 *  Function to calculate the rotation from a body-fixed to a topocentric frame at a given point (localPoint), based on the geocentric
 *  latitude and longitude (i.e. angles based on line from center of body to localPoint. The topocentric calculation differs per type of
 *  body shape wrt which the point is located. The rotation is calculated by concatenating the unit vectors of the topocentric frame,
 *  expresse in the body-fixed frame.
 *  \param bodyShapeModel Shape model of body on which point is located.
 *  \param geocentricLatitude Geocentric latitude of point, i.e. angle between equatorial planet and line from center of body to localPoint
 *  \param geocentricLongitude Geocentric longitude of point, i.e. angle between positive x-axis and line from center of body to localPoint,
 *  projected on equatorial plane, measured positive in +y direction
 *  \param localPoint Cartesian position of point in body-fixed frame.
 *  \return The rotation from the body-fixed to the topocentric frame at localPoint.
 */
Eigen::Quaterniond getRotationQuaternionFromBodyFixedToTopocentricFrame(
        const std::shared_ptr< basic_astrodynamics::BodyShapeModel > bodyShapeModel,
        const double geocentricLatitude,
        const double geocentricLongitude,
        const Eigen::Vector3d localPoint );

struct LinearStationMotionModel : public StationMotionModel {
public:
    //! Reference epoch is in TDB seconds since J2000, matching station environment epochs. The default value is
    //! retained from earlier versions for compatibility (note that it is a Julian day number, not seconds).
    LinearStationMotionModel( const Eigen::Vector3d& linearVelocity,
                              const double referenceEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 ):
        linearVelocity_( linearVelocity ), referenceEpoch_( referenceEpoch )
    {}

    ~LinearStationMotionModel( ) {}

    Eigen::Vector6d getBodyFixedStationMotion( const Time& time,
                                               const std::shared_ptr< ground_stations::GroundStationState > groundStationState = nullptr,
                                               const std::string& targetFrameOrigin = "" )
    {
        return ( Eigen::Vector6d( ) << linearVelocity_ * static_cast< double >( time - Time( referenceEpoch_ ) ), linearVelocity_ )
                .finished( );
    }

protected:
    Eigen::Vector3d linearVelocity_;

    double referenceEpoch_;
};

struct PiecewiseConstantStationMotionModel : public StationMotionModel {
public:
    PiecewiseConstantStationMotionModel( const std::map< double, Eigen::Vector3d >& displacementList )
    {
        displacementTimes_ = utilities::createVectorFromMapKeys( displacementList );
        displacementVectors_ = utilities::createVectorFromMapValues( displacementList );

        firstDisplacementTime_ = displacementTimes_.at( 0 );
        finalDisplacementTime_ = displacementTimes_.at( displacementTimes_.size( ) - 1 );
    }

    ~PiecewiseConstantStationMotionModel( ) {}

    Eigen::Vector6d getBodyFixedStationMotion( const Time& time,
                                               const std::shared_ptr< ground_stations::GroundStationState > groundStationState = nullptr,
                                               const std::string& targetFrameOrigin = "" );

protected:
    double firstDisplacementTime_;

    double finalDisplacementTime_;

    std::vector< double > displacementTimes_;

    std::vector< Eigen::Vector3d > displacementVectors_;
};

//! Equation used to convert station positions from geocentric to barycentric coordinates (IERS Conventions 2010, Chapter 11).
/*!
 *  The conversion of a station position from geocentric to barycentric coordinates depends on the time scales with which the
 *  coordinates are consistent. Geocentric coordinates can be consistent with TCG or with TT, which differ by a constant scale
 *  factor; barycentric coordinates can likewise be consistent with TCB or with TDB. The ITRF uses TT-consistent coordinates.
 */
enum RelativisticStationCoordinateCompatibility {
    //! No equation selected; creating the relativistic station motion model with this value throws an exception.
    unspecified_station_coordinate_compatibility,

    //! TCG-consistent station coordinates to TCB-consistent coordinates (Eq. 11.18).
    unscaled_tcg_to_tcb_station_coordinates,

    //! TT-consistent station coordinates to TDB-consistent coordinates (Eq. 11.19), which adds a constant scale factor.
    tt_to_tdb_station_coordinates
};

//! Body whose gravitational potential at the centre of the station's central body is used in the relativistic station
//! correction. The central body itself must not be included.
struct RelativisticStationExternalBody {
    //! Name of the body, and functions returning its barycentric state at a TDB epoch and its gravitational parameter.
    std::string name;
    std::function< Eigen::Vector6d( const Time& ) > barycentricState;
    std::function< double( ) > gravitationalParameter;
};

//! Relativistic correction to the station position, for station states that are used in a barycentric frame.
/*!
 *  Correction from Eq. 11.18 or 11.19 of the IERS Conventions 2010 (see RelativisticStationCoordinateCompatibility). The
 *  station position relative to the centre of its central body (for instance the Earth) is scaled using the gravitational
 *  potential of the other bodies at that centre, and corrected for the barycentric velocity of the central body. The
 *  smaller terms of the full transformation, which depend on the acceleration of the central body, are not included. The
 *  velocity correction is the time derivative of the position correction. The correction is only applied when the station
 *  state is requested relative to the solar system barycentre (frame origin "SSB"); otherwise it is zero.
 */
struct BodyCentricToBarycentricRelativisticStationMotion : public StationMotionModel {
public:
    BodyCentricToBarycentricRelativisticStationMotion(
            const std::function< Eigen::Vector6d( const Time& ) >& bodyBarycentricStateFunction,
            const std::function< Eigen::Quaterniond( const Time& ) >& inertialToBodyFixedRotationFunction,
            const std::function< Eigen::Matrix3d( const Time& ) >& rotationDerivativeFunction,
            const std::vector< RelativisticStationExternalBody >& externalBodies,
            const RelativisticStationCoordinateCompatibility coordinateCompatibility,
            const double derivativeStep = 16.0 );

    Eigen::Vector6d getBodyFixedStationMotion( const Time& time,
                                               const std::shared_ptr< GroundStationState > groundStationState,
                                               const std::string& targetFrameOrigin ) override;
    Eigen::Vector3d getBodyFixedStationPositionCorrection( const Time& time,
                                                           const std::shared_ptr< GroundStationState > groundStationState,
                                                           const std::string& targetFrameOrigin ) override;

    //! Correction for a station state (or position) that already includes all other station motion. This is used by
    //! CombinedStationMotionModel, which applies this correction after summing all other motion models.
    Eigen::Vector6d correctCompleteLocalState( const Time& time, const Eigen::Vector6d& state, const std::string& targetFrameOrigin );
    Eigen::Vector3d correctCompleteLocalPosition( const Time& time, const Eigen::Vector3d& position, const std::string& targetFrameOrigin );
    const std::vector< RelativisticStationExternalBody >& getExternalBodies( ) const
    {
        return externalBodies_;
    }

private:
    Eigen::Matrix3d correctionMatrix( const Time& time );
    Eigen::Matrix3d correctionMatrixDerivative( const Time& time );
    double potential( const Time& time );
    double potentialDerivative( const Time& time );

    //! Functions returning the barycentric state of the central body, its rotation from the inertial to the body-fixed
    //! frame, and the time derivative of that rotation (if empty, the derivative is computed numerically).
    std::function< Eigen::Vector6d( const Time& ) > bodyBarycentricStateFunction_;
    std::function< Eigen::Quaterniond( const Time& ) > inertialToBodyFixedRotationFunction_;
    std::function< Eigen::Matrix3d( const Time& ) > rotationDerivativeFunction_;

    //! Bodies whose gravitational potential is included in the correction.
    std::vector< RelativisticStationExternalBody > externalBodies_;

    //! Equation used for the correction (see RelativisticStationCoordinateCompatibility).
    RelativisticStationCoordinateCompatibility coordinateCompatibility_;

    //! Time step (seconds) of the numerical derivatives used for the velocity correction.
    double derivativeStep_;
};

struct CustomStationMotionModel : public StationMotionModel {
public:
    CustomStationMotionModel( const std::function< Eigen::Vector6d( const double ) > customDisplacementModel ):
        customDisplacementModel_( customDisplacementModel )
    {}

    ~CustomStationMotionModel( ) {}

    Eigen::Vector6d getBodyFixedStationMotion( const Time& time,
                                               const std::shared_ptr< ground_stations::GroundStationState > groundStationState = nullptr,
                                               const std::string& targetFrameOrigin = "" )
    {
        return customDisplacementModel_( static_cast< double >( time ) );
    }

protected:
    std::function< Eigen::Vector6d( const double ) > customDisplacementModel_;
};

//! Combination of several station motion models. The changes in station state from all models are added, except for the
//! relativistic correction (at most one), which is applied last, to the station state that includes all other motion. The
//! position of the relativistic correction in the list of models therefore does not matter.
struct CombinedStationMotionModel : public StationMotionModel {
public:
    CombinedStationMotionModel( const std::vector< std::shared_ptr< StationMotionModel > >& modelList );

    Eigen::Vector6d getBodyFixedStationMotion( const Time& time,
                                               const std::shared_ptr< GroundStationState > groundStationState,
                                               const std::string& targetFrameOrigin ) override;
    Eigen::Vector3d getBodyFixedStationPositionCorrection( const Time& time,
                                                           const std::shared_ptr< GroundStationState > groundStationState,
                                                           const std::string& targetFrameOrigin ) override;

private:
    void appendModel( const std::shared_ptr< StationMotionModel >& model );

    //! Station motion models whose contributions are added, and the relativistic correction applied last (or nullptr).
    std::vector< std::shared_ptr< StationMotionModel > > modelList_;
    std::shared_ptr< BodyCentricToBarycentricRelativisticStationMotion > terminalModel_;
};

}  // namespace ground_stations

}  // namespace tudat
#endif  // TUDAT_GROUNDSTATIONSTATE_H
