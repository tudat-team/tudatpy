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
#include <limits>
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
 *  z-axes, respectively. The altitude is the shortest distance from the point to the surface.
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
     *  Function to calculate the altitude above the ellipsoid from a body fixed position, defined as the shortest
     *  distance to the surface. The altitude is negative for points inside the ellipsoid.
     *
     *  The closest surface point p of a position x satisfies p_i = a_i^2 x_i / ( a_i^2 + t ), with a_i the semi-axes,
     *  so that x - p = t n, with n_i = p_i / a_i^2 the (non-unit) surface normal. The value of t follows from requiring
     *  p to be on the surface, and is solved in terms of u = t + c_min^2, with c_min the smallest semi-axis:
     *  sum_i ( a_i x_i / ( a_i^2 - c_min^2 + u ) )^2 = 1. The left-hand side is convex and decreasing for u > 0, so
     *  Newton iterations from a lower bound of the root converge monotonically. If the left-hand side is below 1 at
     *  u = 0 (only possible for an interior point in the plane normal to the smallest axis), there is no root, and the
     *  closest point is out of this plane, with t = -c_min^2 and the remainder of the surface equation along the
     *  smallest axis.
     *  \param bodyFixedPosition Cartesian, body-fixed position of the point at which the altitude
     *  is to be determined.
     *  \return Altitude above the ellipsoid.
     */
    double getAltitude( const Eigen::Vector3d& bodyFixedPosition )
    {
        const Eigen::Array3d semiAxes( semiAxisA_, semiAxisB_, semiAxisC_ );
        const Eigen::Array3d scaledPosition = semiAxes * bodyFixedPosition.array( );
        const double smallestSquaredSemiAxis = semiAxes.square( ).minCoeff( );
        const Eigen::Array3d squaredSemiAxisOffsets = semiAxes.square( ) - smallestSquaredSemiAxis;

        // Lower bounds of the root: from the terms of the smallest semi-axis alone, and from replacing all offsets by
        // the largest one (which is close to the root far from the body).
        double u = std::max( std::sqrt( ( squaredSemiAxisOffsets == 0.0 ).select( scaledPosition.square( ), 0.0 ).sum( ) ),
                             std::sqrt( scaledPosition.square( ).sum( ) ) - squaredSemiAxisOffsets.maxCoeff( ) );

        Eigen::Array3d surfaceNormal;
        double residual = 0.0;
        for( int iteration = 0; iteration < 100; iteration++ )
        {
            // Denominator is only limited for u = 0 along the smallest axis, where the position component is zero
            const Eigen::Array3d denominators = ( squaredSemiAxisOffsets + u ).max( std::numeric_limits< double >::min( ) );
            surfaceNormal = bodyFixedPosition.array( ) / denominators;

            const Eigen::Array3d squaredTerms = ( semiAxes * surfaceNormal ).square( );
            residual = squaredTerms.sum( ) - 1.0;
            if( residual <= 0.0 )
            {
                break;
            }

            const double step = residual / ( 2.0 * ( squaredTerms / denominators ).sum( ) );
            if( u + step == u )
            {
                break;
            }
            u += step;
        }

        return ( u - smallestSquaredSemiAxis ) *
                std::sqrt( surfaceNormal.square( ).sum( ) + std::max( -residual, 0.0 ) / smallestSquaredSemiAxis );
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
