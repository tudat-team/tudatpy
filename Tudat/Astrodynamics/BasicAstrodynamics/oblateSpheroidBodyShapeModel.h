#ifndef TUDAT_OBLATESPHEROIDBODYSHAPEMODEL_H
#define TUDAT_OBLATESPHEROIDBODYSHAPEMODEL_H

/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      150409    D. Dirkx          Migrated from personal code.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 *    Notes
 *
 */


#include <Eigen/Core>
#include <iostream>

#include <boost/lambda/lambda.hpp>

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
        return ( ( 2.0 * equatorialRadius_+ polarRadius_) / 3.0 );
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
