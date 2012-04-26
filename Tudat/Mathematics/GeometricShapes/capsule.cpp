/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      102511    D. Dirkx          First version of file.
 *      110120    D. Dirkx          Finalized for code check.
 *      110208    K. Kumar          Updated file header; corrected Doxygen comments; minor changes.
 *      110209    D. Dirkx          Minor changes.
 *      110209    K. Kumar          Minor changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120323    D. Dirkx          Removed set functions; moved functionality to constructor.
 *
 *    References
 *      E.H. Hirschel and C. Weiland, Selected Aerothermodynamic Design Problems
 *          of Hypersonic Flight Vehicles (chapter 5), Springer/AIAA, 2009.
 *      D. Dirkx, Continuous Shape Optimization of Entry Vehicles, MSc thesis,
 *          Delft University of Technology, 2011 (Unpublished).
 *
 */

#include <cmath>
#include <iostream>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/GeometricShapes/capsule.h"
#include "Tudat/Mathematics/GeometricShapes/conicalFrustum.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"
#include "Tudat/Mathematics/GeometricShapes/torus.h"

namespace tudat
{
namespace mathematics
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
    using tudat::mathematics::PI;

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
    boost::shared_ptr< SphereSegment > noseSphere_ = boost::make_shared< SphereSegment >(
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
    boost::shared_ptr< ConicalFrustum > cone_ = boost::make_shared< ConicalFrustum >(
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
    boost::shared_ptr< SphereSegment > rearSphere_ = boost::make_shared< SphereSegment >(
                rearNoseRadius_, 0.0, 2.0 * PI, PI / 2.0 - rearAngle_, PI );

    // Set translation vector of rear sphere.
    translationVector_( 0 ) =  ( rearNoseRadius_ * sin( -rearAngle_ ) ) - rearLength_
            - ( sideRadius_ * ( sin( PI / 2.0 - noseSphereAngle_ ) + sin ( -rearAngle_ ) ) );
    rearSphere_->setOffset( translationVector_ );
    setSingleSurfaceGeometry( rearSphere_, 2 );

    // Create torus section of capsule.
    double torusMajorRadius_ = ( noseRadius_ - sideRadius_ ) * sin( noseSphereAngle_ );
    boost::shared_ptr< Torus > torus_ = boost::make_shared< Torus >(
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
std::ostream &operator<<( std::ostream &stream, Capsule& capsule )
{
    using std::endl;

    stream << "This is a capsule." << endl;
    stream << "The defining parameters are: "<< endl
           << "Nose radius: " << capsule.getNoseRadius( ) << endl
           << "Mid radius: " << capsule.getMiddleRadius( ) << endl
           << "Rear length: " << capsule.getRearLength( ) << endl
           << "Rear angle: " << capsule.getRearAngle( ) << endl
           << "Side radius: " << capsule.getSideRadius( )<< endl;

    return stream;
}

} // namespace geometric_shapes
} // namespace mathematics
} // namespace tudat
