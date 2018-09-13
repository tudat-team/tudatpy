#ifndef TUDAT_OBLATESPHEROIDBODYSHAPEMODEL_H
#define TUDAT_OBLATESPHEROIDBODYSHAPEMODEL_H

/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */


#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! Body shape model for an oblate spheroid
/*!
 *  Body shape model for an oblate spheroid (flattened sphere), typically used as approximation for
 *  planets and large moons.
 */
class OblateSpheroidBodyShapeModel: public BodyShapeModel
{
public:

    //! Constructor
    /*!
     *  Constructor, sets the geomtric properties of the shape.
     *  \param equatorialRadius Equatorial radius of the oblate spheroid
     *  \param flattening Flattening of the oblate spheroid
     */
    OblateSpheroidBodyShapeModel( const double equatorialRadius, const double flattening ):
        equatorialRadius_( equatorialRadius ), flattening_( flattening )
    {
        // Calculate and set polar radius.
        polarRadius_ = equatorialRadius * ( 1.0 - flattening_ );
    }

    //! Destructor
    ~OblateSpheroidBodyShapeModel( ){ }

    //! Calculates the altitude above the oblate spheroid
    /*!
     *  Function to calculate the altitude above the oblate spheroid from a body fixed position.
     *  \param bodyFixedPosition Cartesian, body-fixed position of the point at which the altitude
     *  is to be determined.
     *  \return Altitude above the oblate spheroid.
     */
    double getAltitude( const Eigen::Vector3d& bodyFixedPosition )
    {
        return coordinate_conversions::calculateAltitudeOverOblateSpheroid(
                    bodyFixedPosition, equatorialRadius_, flattening_, 1.0E-4 );
    }

    //! Calculates the geodetic position w.r.t. the oblate spheroid.
    /*!
     *  Function to calculate the geodetic position w.r.t. the oblate spheroid.
     *  \sa convertCartesianToGeodeticCoordinates
     *  \param bodyFixedPosition Cartesian, body-fixed position of the point at which the geodetic
     *  position is to be determined.
     *  \param tolerance Convergence criterion for iterative algorithm that is employed. Represents
     *  the required change of position (in m) between two iterations.
     *  \return Geodetic coordinates at requested point.
     */
    Eigen::Vector3d getGeodeticPositionWrtShape( const Eigen::Vector3d& bodyFixedPosition,
                                        const double tolerance = 1.0E-4 )
    {
        return coordinate_conversions::convertCartesianToGeodeticCoordinates(
                    bodyFixedPosition, equatorialRadius_, flattening_, tolerance );
    }

    //! Calculates the geodetic latitude w.r.t. the oblate spheroid.
    /*!
     *  Function to calculate the geodetic latitude w.r.t. the oblate spheroid.
     *  \sa convertCartesianToGeodeticCoordinates
     *  \param bodyFixedPosition Cartesian, body-fixed position of the point at which the geodetic
     *  latitude is to be determined.
     *  \param tolerance Convergence criterion for iterative algorithm that is employed. Represents
     *  the required change of position (in m) between two iterations.
     *  \return Geodetic latitude at requested point.
     */
    double getGeodeticLatitude( const Eigen::Vector3d& bodyFixedPosition,
                                        const double tolerance = 1.0E-4 )
    {
        return coordinate_conversions::calculateGeodeticLatitude(
                    bodyFixedPosition, equatorialRadius_, flattening_, tolerance );
    }

    //! Function to return the mean radius of the oblate spheroid.
    /*!
     *  Function to return the mean radius of the oblate spheroid.
     *  \return Average radius of oblate spheroid.
     */
    double getAverageRadius( )
    {
        return ( ( 2.0 * equatorialRadius_ + polarRadius_ ) / 3.0 );
    }

    //! Function to obtain the equatorial radius
    /*!
     *  Function to obtain the equatorial radius
     *  \return Equatorial radius of the oblate spheroid
     */
    double getEquatorialRadius( )
    {
        return equatorialRadius_;
    }

    //! Function to obtain the flattening of the oblate spheroid
    /*!
     *  Function to obtain the flattening of the oblate spheroid
     *  \return Flattening radius of the oblate spheroid
     */
    double getFlattening( )
    {
        return flattening_;
    }

private:
    //! Equatorial radius of the oblate spheroid
    double equatorialRadius_;

    //! Polar radius of the oblate spheroid
    double polarRadius_;

    //! Flattening of the oblate spheroid
    double flattening_;
};

} // namespace basic_astrodynamics
} // namespace tudat


#endif // TUDAT_OBLATESPHEROIDBODYSHAPEMODEL_H
