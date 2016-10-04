#ifndef NOMINALGROUNDSTATIONSTATE_H
#define NOMINALGROUNDSTATIONSTATE_H

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


namespace tudat
{

namespace ground_stations
{

//! Class containing for storing of and calculations on the body-fixed position of a ground station
/*!
 *  Class containing for storing of and calculations on the body-fixed position of a ground station
 */
class NominalGroundStationState
{
public:
    //! Constructor taking cartesian position data.
    /*!
     *  Constructor taking cartesian position data.
     */
    NominalGroundStationState(
            const Eigen::Vector3d stationCartesianPosition,
            const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface = boost::shared_ptr< basic_astrodynamics::BodyShapeModel >( ),
            const double referenceJulianYear = basic_astrodynamics::JULIAN_DAY_ON_J2000 / physical_constants::JULIAN_YEAR_IN_DAYS,
            const bool setTransformation = 1 );

    virtual ~NominalGroundStationState( ){ }

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

    virtual void resetGroundStationPositionAtEpoch( const Eigen::Vector3d cartesianPosition );

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
    Eigen::Vector3d sphericalPosition_; //(radius, co-geocentric latitude, longitude)

    //! Name of ground station of which this object is a member
    /*!
     *  Name of ground station of which this object is a member, i.e. the station of which this object describes the state.
     */
    std::string siteId_;

    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface_;
};

}

}
#endif // NOMINALGROUNDSTATIONSTATE_H
