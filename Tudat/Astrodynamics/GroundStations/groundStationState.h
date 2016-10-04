#ifndef TUDAT_GROUNDSTATIONSTATE_H
#define TUDAT_GROUNDSTATIONSTATE_H

#include <vector>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"

#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateRepresentationConversions.h"


namespace tudat
{

namespace ground_stations
{

//! Class containing for storing of and calculations on the body-fixed position of a ground station
/*!
 *  Class containing for storing of and calculations on the body-fixed position of a ground station
 */
class GroundStationState
{
public:
    //! Constructor taking cartesian position data.
    /*!
     *  Constructor taking cartesian position data.
     */
    GroundStationState(
            const Eigen::Vector3d stationPosition,
            const coordinate_conversions::PositionElementTypes inputElementType = coordinate_conversions::cartesian_position,
            const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface =
            boost::shared_ptr< basic_astrodynamics::BodyShapeModel >( ) );

    virtual ~GroundStationState( ){ }

    //! Function to obtain the Cartesian position of the state in the local frame at a given time.
    /*!
     *  Function to obtain the Cartesian position of the state in the local (body-fixed, not topocentrix) frame at a given time.
     *  Adds all position variations to the nominal state (at the requested time) and returns the state.
     *  \param secondsSinceEpoch Secons since reference epoch at which the position is to be retrieved.
     *  \param inputReferenceEpoch Reference epoch julian day
     *  \return Cartiesian position of station in local frame at requested time.
     */
    Eigen::Vector3d getCartesianPositionInTime(
            const double secondsSinceEpoch,
            const double inputReferenceEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    //! Function to return the nominal cartesian (unperturbed) position of the station
    /*!
     *  Function to return the nominal cartesian (unperturbed, i.e. not including linear drift, eccentricity, tidal variations, etc.)
     *  position of the station in body-fixed system. Typically used for applications require moderate to low precision of ground station position.
     *  \return Unperturbed position of the station
     */
    Eigen::Vector3d getNominalCartesianPosition( )
    {
        return cartesianPosition_;
    }

    //! Function to return the nominal spherical (unperturbed) position of the station
    /*!
     *  Function to return the nominal spherical (unperturbed, i.e. not including linear drift, eccentricity, tidal variations, etc.)
     *  position of the station in body-fixed system.
     *  \return Unperturbed position of the station
     */
    Eigen::Vector3d getNominalSphericalPosition( )
    {
        return sphericalPosition_;
    }

    Eigen::Vector3d getNominalGeodeticPosition( )
    {
        return geodeticPosition;
    }

    virtual void resetGroundStationPositionAtEpoch(
            const Eigen::Vector3d stationPosition,
            const coordinate_conversions::PositionElementTypes inputElementType = coordinate_conversions::cartesian_position );

    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > getBodySurface( )
    {
        return bodySurface_;
    }

protected:

    //! Cartesian position of station
    /*!
     *  Cartesian position of station, without variations (linear drift, eccentricity, tides, etc. ), in the body-fixed frame.
     */
    Eigen::Vector3d cartesianPosition_; //(x,y,z)

    //! Spherical position of station
    /*!
     *  Spherical position of station, without variations (linear drift, eccentricity, tides, etc. ), in the body-fixed frame.
     *  The order of the entries is: radius, colatitude, longitude.
     */
    Eigen::Vector3d sphericalPosition_; //(radius, geocentric latitude, longitude)

    Eigen::Vector3d geodeticPosition; //(altitude, geodetic latitude, longitude)


    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface_;
};

}

}
#endif // TUDAT_GROUNDSTATIONSTATE_H
