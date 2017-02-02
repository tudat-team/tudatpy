/*    Copyright (c) 2010-2016, Delft University of Technology
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

//! Class storing and computing the (time-variable) state of a ground station in a body-fixed frame.
class GroundStationState
{
public:

    //! Constructor
    /*!
     *  Constructor
     * \param stationPosition Position of ground station in body-fixed frame
     * \param inputElementType Element type (e.g. Cartesian, spherical, etc.) of groundStationPosition.
     * \param bodySurface Shape of body on which state is defined. If NULL (default), no conversions to/from
     * geodetic position are possible.
     */
    GroundStationState(
            const Eigen::Vector3d stationPosition,
            const coordinate_conversions::PositionElementTypes inputElementType = coordinate_conversions::cartesian_position,
            const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface =
            boost::shared_ptr< basic_astrodynamics::BodyShapeModel >( ) );

    virtual ~GroundStationState( ){ }

    //! Function to obtain the Cartesian state of the ground station in the local frame at a given time.
    /*!
     *  Function to obtain the Cartesian state of the ground station in the local frame (body-fixed, not topocentric) at a
     *  given time.  Adds all position variations to the nominal state (at the requested time) and returns the state. NOTE:
     *  poisition variations are as yet not included.
     *  \param secondsSinceEpoch Secons since reference epoch at which the position is to be retrieved.
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
      *  \param secondsSinceEpoch Secons since reference epoch at which the position is to be retrieved.
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
     *  NOTE: This vector is only defined in bodySurface_ was non-NULL at last setting of nominal state
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
    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > getBodySurface( )
    {
        return bodySurface_;
    }

protected:

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
     *  non-NULL at last setting of nominal state (call to resetGroundStationPositionAtEpoch).
     */
    Eigen::Vector3d geodeticPosition;

    //! Shape of body on which state is defined
    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface_;
};

} // namespace ground_stations

} // namespace tudat
#endif // TUDAT_GROUNDSTATIONSTATE_H
