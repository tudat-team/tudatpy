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
 *      E.H. Hirschel and C. Weiland, Selected Aerothermodynamic Design Problems of Hypersonic
 *          Flight Vehicles (chapter 5), Springer/AIAA, 2009.
 *      D. Dirkx, Continuous Shape Optimization of Entry Vehicles, MSc thesis, Delft University
 *          of Technology, 2011 (Unpublished).
 *
 */

#include <cmath>
#include <iostream>

#include <boost/make_shared.hpp>
#include <memory>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Mathematics/GeometricShapes/capsule.h"
#include "Tudat/Mathematics/GeometricShapes/conicalFrustum.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"
#include "Tudat/Mathematics/GeometricShapes/torus.h"

namespace tudat
{
namespace geometric_shapes
{

//! Default constructor.
Capsule::Capsule( const double noseRadius,
                  const double middleRadius,
                  const double rearLength,
                  const double rearAngle,
                  const double sideRadius )
{
    using std::sin;
    using std::cos;
    using mathematical_constants::PI;

    // Call set functions for number of single and composite surface geometries
    // with predetermined values.
    setNumberOfCompositeSurfaceGeometries( 0 );
    setNumberOfSingleSurfaceGeometries( 4 );

    // Set member shape variables
    noseRadius_ = noseRadius;
    middleRadius_ = middleRadius;
    rearLength_ = rearLength;
    rearAngle_ = rearAngle;
    sideRadius_ = sideRadius;

    // Determine and set extent of spherical nose part.
    double noseSphereAngle_ = asin( ( middleRadius_ - sideRadius_ )
                                    / ( noseRadius_ - sideRadius_ ) );

    // Create nose sphere.
    std::shared_ptr< SphereSegment > noseSphere_ = std::make_shared< SphereSegment >(
                noseRadius_, 0, 2 * PI, 0, noseSphereAngle_ );

    // Declare translation vector.
    Eigen::VectorXd translationVector_ = Eigen::VectorXd( 3 );
    translationVector_( 2 ) = 0.0;
    translationVector_( 1 ) = 0.0;
    translationVector_( 0 ) = - noseRadius_ * cos( noseSphereAngle_ );

    // Set nose translation vector.
    noseSphere_->setOffset( translationVector_ );

    // Set noseSphere_ in singleSurfaceList_.
    setSingleSurfaceGeometry( noseSphere_, 0 );

    // Create rear cone, fully revolved.
    std::shared_ptr< ConicalFrustum > cone_ = std::make_shared< ConicalFrustum >(
                rearAngle_, middleRadius_ - sideRadius_ * ( 1.0 - cos( rearAngle_ ) ),
                rearLength_ );

    // Set translation vector of cone.
    translationVector_( 0 ) = -sideRadius_ * ( sin( PI / 2.0 - noseSphereAngle_ )
                                               + sin ( -rearAngle_ ) );
    cone_->setOffset( translationVector_ );

    // Set cone in singleSurfaceList_.
    setSingleSurfaceGeometry( cone_, 1 );

    // Calculate end radius of cone.
    double endRadius_ = cone_->getStartRadius( ) + rearLength_ * tan( rearAngle_ );

    // Calculate rear sphere radius.
    double rearNoseRadius_ = endRadius_ / cos( -rearAngle_ );

    // Create rear sphere ( "end cap" ), fully revolved.
    std::shared_ptr< SphereSegment > rearSphere_ = std::make_shared< SphereSegment >(
                rearNoseRadius_, 0.0, 2.0 * PI, PI / 2.0 - rearAngle_, PI );

    // Set translation vector of rear sphere.
    translationVector_( 0 ) =  ( rearNoseRadius_ * sin( -rearAngle_ ) ) - rearLength_
            - ( sideRadius_ * ( sin( PI / 2.0 - noseSphereAngle_ ) + sin ( -rearAngle_ ) ) );
    rearSphere_->setOffset( translationVector_ );
    setSingleSurfaceGeometry( rearSphere_, 2 );

    // Create torus section of capsule.
    double torusMajorRadius_ = ( noseRadius_ - sideRadius_ ) * sin( noseSphereAngle_ );
    std::shared_ptr< Torus > torus_ = std::make_shared< Torus >(
       torusMajorRadius_, sideRadius_, 0.0, 2.0 * PI, PI / 2.0 - noseSphereAngle_, rearAngle_ );

    // Set translation vector of rear sphere.
    translationVector_( 0 ) = -cos( noseSphereAngle_ ) * sideRadius_;
    torus_->setOffset( translationVector_ );
    setSingleSurfaceGeometry( torus_, 3 );

    // Set rotation matrix fo each part to be compatible with flow direction in
    // aerodynamic analysis.
    Eigen::MatrixXd rotationMatrix = Eigen::MatrixXd( 3, 3 );
    double angle_ = PI / 2.0;
    rotationMatrix( 0, 0 ) = cos( angle_ );
    rotationMatrix( 0, 1 ) = 0.0;
    rotationMatrix( 0, 2 ) = sin( angle_ );
    rotationMatrix( 1, 0 ) = 0.0;
    rotationMatrix( 1, 1 ) = 1.0;
    rotationMatrix( 1, 2 ) = 0.0;
    rotationMatrix( 2, 0 ) = -sin( angle_ );
    rotationMatrix( 2, 1 ) = 0.0;
    rotationMatrix( 2, 2 ) = cos( angle_ );

    // Set rotation matrix for single surface geometries.
    for ( unsigned i = 0; i < numberOfSingleSurfaceGeometries_ ; i++ )
    {
        singleSurfaceGeometryList_[ i ]->setRotationMatrix( rotationMatrix );
    }
}

//! Overload ostream to print class information.
std::ostream &operator << ( std::ostream &stream, Capsule& capsule )
{
    stream << "This is a capsule." << std::endl;
    stream << "The defining parameters are: " << std::endl
           << "Nose radius: " << capsule.getNoseRadius( ) << std::endl
           << "Mid radius: " << capsule.getMiddleRadius( ) << std::endl
           << "Rear length: " << capsule.getRearLength( ) << std::endl
           << "Rear angle: " << capsule.getRearAngle( ) << std::endl
           << "Side radius: " << capsule.getSideRadius( ) << std::endl;

    return stream;
}

} // namespace geometric_shapes
} // namespace tudat
