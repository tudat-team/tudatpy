#ifndef TUDAT_TRIAXIALELLIPSOIDBODYSHAPEMODEL_H
#define TUDAT_TRIAXIALELLIPSOIDBODYSHAPEMODEL_H

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/bodyShapeModel.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! Body shape model defined as a tri-axial ellipsoid
/*!
 *  Body shape model defined as a tri-axial ellipsoid, with semi-axes a, b and c aligned with the body-fixed x-, y- and
 *  z-axes, respectively. The altitude is computed along the line connecting the point to the center of the body (i.e.
 *  the distance from the point to the intersection of this line with the surface). Note that this is not the shortest
 *  distance to the surface: for elongated bodies, and close to the surface, it is larger than the shortest distance.
 */
class TriAxialEllipsoidBodyShapeModel : public BodyShapeModel
{
public:
    //! Constructor.
    /*!
     *  Constructor, defines body shape.
     *  \param semiAxisA Semi-axis of the ellipsoid along the body-fixed x-axis
     *  \param semiAxisB Semi-axis of the ellipsoid along the body-fixed y-axis
     *  \param semiAxisC Semi-axis of the ellipsoid along the body-fixed z-axis
     */
    TriAxialEllipsoidBodyShapeModel( const double semiAxisA, const double semiAxisB, const double semiAxisC ):
        semiAxisA_( semiAxisA ), semiAxisB_( semiAxisB ), semiAxisC_( semiAxisC )
    {
        if( !( semiAxisA_ > 0.0 ) || !( semiAxisB_ > 0.0 ) || !( semiAxisC_ > 0.0 ) )
        {
            throw std::runtime_error( "Error when creating tri-axial ellipsoid shape model, all semi-axes must be positive." );
        }
    }

    //! Destructor
    ~TriAxialEllipsoidBodyShapeModel( ) {}

    //! Calculates the altitude above the ellipsoid
    /*!
     *  Function to calculate the altitude above the ellipsoid from a body fixed position, measured along the line
     *  connecting the point to the center of the body. The altitude is negative for points inside the ellipsoid. At
     *  the center of the body, the direction is undefined, and minus the smallest semi-axis is returned.
     *  \param bodyFixedPosition Cartesian, body-fixed position of the point at which the altitude
     *  is to be determined.
     *  \return Altitude above the ellipsoid.
     */
    double getAltitude( const Eigen::Vector3d& bodyFixedPosition )
    {
        const double distanceToCenter = bodyFixedPosition.norm( );
        if( distanceToCenter == 0.0 )
        {
            return -std::min( { semiAxisA_, semiAxisB_, semiAxisC_ } );
        }
        return distanceToCenter - getSurfaceRadius( bodyFixedPosition );
    }

    //! Calculates the distance from the center to the surface, in the direction of the given position
    /*!
     *  Function to calculate the distance from the center of the body to the surface of the ellipsoid, along the
     *  direction of the given body-fixed position.
     *  \param bodyFixedPosition Cartesian, body-fixed position defining the direction (must be non-zero).
     *  \return Distance from the center to the surface in the given direction.
     */
    double getSurfaceRadius( const Eigen::Vector3d& bodyFixedPosition )
    {
        const double scaledSquaredDistance = std::pow( bodyFixedPosition.x( ) / semiAxisA_, 2 ) +
                std::pow( bodyFixedPosition.y( ) / semiAxisB_, 2 ) + std::pow( bodyFixedPosition.z( ) / semiAxisC_, 2 );
        return bodyFixedPosition.norm( ) / std::sqrt( scaledSquaredDistance );
    }

    //! Function to return the mean radius of the ellipsoid.
    /*!
     *  Function to return the mean radius of the ellipsoid, defined as the radius of a sphere with the same volume.
     *  \return Average (volumetric) radius of the ellipsoid.
     */
    double getAverageRadius( )
    {
        return std::cbrt( semiAxisA_ * semiAxisB_ * semiAxisC_ );
    }

    //! Function to obtain the semi-axis along the body-fixed x-axis
    double getSemiAxisA( )
    {
        return semiAxisA_;
    }

    //! Function to obtain the semi-axis along the body-fixed y-axis
    double getSemiAxisB( )
    {
        return semiAxisB_;
    }

    //! Function to obtain the semi-axis along the body-fixed z-axis
    double getSemiAxisC( )
    {
        return semiAxisC_;
    }

    //! Function to obtain the three semi-axes (a, b, c)
    Eigen::Vector3d getSemiAxes( )
    {
        return Eigen::Vector3d( semiAxisA_, semiAxisB_, semiAxisC_ );
    }

private:
    //! Semi-axis along the body-fixed x-axis
    double semiAxisA_;

    //! Semi-axis along the body-fixed y-axis
    double semiAxisB_;

    //! Semi-axis along the body-fixed z-axis
    double semiAxisC_;
};

}  // namespace basic_astrodynamics
}  // namespace tudat

#endif  // TUDAT_TRIAXIALELLIPSOIDBODYSHAPEMODEL_H
