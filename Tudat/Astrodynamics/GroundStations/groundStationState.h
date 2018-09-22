/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateRepresentationConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"


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
std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors(
            const Eigen::Matrix3d& topocentricToPlanetFixedFrameMatrix );

//! Function to generate unit vectors of topocentric frame.
/*!
 *  Function to generate unit vectors of topocentric frame. The set of unit vectors of topocentric frame is expressed in
 *  body-fixed frame, in ENU (Earth-North-Up) order.
 *  \param latitude Latitude of point for which topocentric frame is to be created.
 *  \param longitude Longitude of point for which  topocentric frame is to be created.
 *  \return Unit vectors of topocentric frame, expressed in body-fixed frame
 */
std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors(
        const double latitude,  const double longitude );

//! Class storing and computing the (time-variable) state of a ground station in a body-fixed frame.
class GroundStationState
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
    GroundStationState(
            const Eigen::Vector3d stationPosition,
            const coordinate_conversions::PositionElementTypes inputElementType = coordinate_conversions::cartesian_position,
            const std::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface =
            std::shared_ptr< basic_astrodynamics::BodyShapeModel >( ) );

    virtual ~GroundStationState( ){ }

    //! Function to obtain the Cartesian state of the ground station in the local frame at a given time.
    /*!
     *  Function to obtain the Cartesian state of the ground station in the local frame (body-fixed, not topocentric) at a
     *  given time.  Adds all position variations to the nominal state (at the requested time) and returns the state. NOTE:
     *  poisition variations are as yet not included.
     *  \param secondsSinceEpoch Seconds since reference epoch at which the position is to be retrieved.
     *  \param inputReferenceEpoch Reference epoch julian day
     *  \return Cartesian state of station in local frame at requested time.
     */
     Eigen::Vector6d getCartesianStateInTime(
            const double secondsSinceEpoch,
            const double inputReferenceEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

     //! Function to obtain the Cartesian position of the ground station in the local frame at a given time.
     /*!
      *  Function to obtain the Cartesian position of the ground station in the local frame (body-fixed, not topocentric) at a
      *  given time.  Adds all position variations to the nominal state (at the requested time) and returns the state. NOTE:
      *  poisition variations are as yet not included.
      *  \param secondsSinceEpoch Seconds since reference epoch at which the position is to be retrieved.
      *  \param inputReferenceEpoch Reference epoch julian day
      *  \return Cartesian position of station in local frame at requested time.
      */
     Eigen::Vector3d getCartesianPositionInTime(
            const double secondsSinceEpoch,
            const double inputReferenceEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 )
     {
         return getCartesianStateInTime( secondsSinceEpoch, inputReferenceEpoch ).segment( 0, 3 );
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
    double getLatitude( )
    {
        return sphericalPosition_.y( );
    }

    //! Function to return the geocentric longtude of the station.
    /*!
     *  Function to return the geocentric longtude of the station.
     *  \return Geocentric longtude of the station.
     */
    double getLongitude( )
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
        return bodySurface_;
    }

protected:


    //! Function to reset the rotation from the body-fixed to local topocentric frame, and associated unit vectors
    void setTransformationAndUnitVectors( );

    //! Cartesian position of station
    /*!
     *  Cartesian position of station, without variations (linear drift, eccentricity, tides, etc.), in the body-fixed frame.
     */
    Eigen::Vector3d cartesianPosition_;

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
    std::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface_;
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

} // namespace ground_stations

} // namespace tudat
#endif // TUDAT_GROUNDSTATIONSTATE_H
