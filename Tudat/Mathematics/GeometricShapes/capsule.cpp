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
 *      110208    K. Kumar          Updated file header; corrected Doxygen
 *                                  comments; minor changes.
 *      110209    D. Dirkx          Minor changes.
 *      110209    K. Kumar          Minor changes.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *
 *    References
 *      E.H. Hirschel and C. Weiland, Selected Aerothermodynamic Design Problems
 *      of Hypersonic Flight Vehicles (chapter 5), Springer/AIAA, 2009.
 *
 *      D. Dirkx, Continuous Shape Optimization of Entry Vehicles, MSc thesis,
 *      Delft University of Technology, 2011 (Unpublished).
 *
 */

#include <cmath>
#include <Eigen/Core>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/GeometricShapes/capsule.h"
#include "Tudat/Mathematics/GeometricShapes/conicalFrustum.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"
#include "Tudat/Mathematics/GeometricShapes/torus.h"

namespace tudat
{

// Using declarations.
using tudat::mathematics::PI;
using std::cerr;
using std::endl;
using std::sin;
using std::cos;

//! Default constructor.
Capsule::Capsule( ) : middleRadius_( -0.0 ), noseRadius_( -0.0 ), rearLength_( -0.0 ),
    sideRadius_( -0.0 ), rearAngle_( -0.0 )
{
    // Call set functions for number of single and composite surface geometries
    // with predetermined values.
    setNumberOfCompositeSurfaceGeometries( 0 );
    setNumberOfSingleSurfaceGeometries( 4 );

    // Set single surface list contents to NULL to prevent destructor issues.
    for ( unsigned i = 0; i < numberOfSingleSurfaceGeometries_; i++ )
    {
        singleSurfaceGeometryList_[ i ] = NULL;
    }
}

//! Create capsule.
void Capsule::setCapsule( )
{
    // Create nose sphere.
    SphereSegment* noseSphere_ = new SphereSegment( );

    // Set predetermined bounds.
    noseSphere_->setMinimumAzimuthAngle( 0.0 );
    noseSphere_->setMaximumAzimuthAngle( 2.0 * PI );
    noseSphere_->setMinimumZenithAngle( 0.0 );

    // Determine and set extent of spherical nose part.
    double noseSphereAngle_ = asin( ( middleRadius_ - sideRadius_ )
                                    / ( noseRadius_ - sideRadius_ ) );
    noseSphere_->setMaximumZenithAngle( noseSphereAngle_ );

    // Set nose radius.
    noseSphere_->setRadius( noseRadius_ );

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
    ConicalFrustum* cone_ = new ConicalFrustum( );
    cone_->setMinimumAzimuthAngle( 0.0 );
    cone_->setMaximumAzimuthAngle( 2.0 * PI );

    // Set cone start radius.
    cone_->setStartRadius( middleRadius_ - sideRadius_ * ( 1.0 - cos( rearAngle_ ) ) );

    // Set cone length.
    cone_->setLength( rearLength_ );

    // Set cone half angle.
    cone_->setConeHalfAngle( rearAngle_ );

    // Set translation vector of cone.
    translationVector_( 0 ) = -sideRadius_ * ( sin( PI / 2.0 - noseSphereAngle_ )
                                               + sin ( -rearAngle_ ) );
    cone_->setOffset( translationVector_ );

    // Set cone in singleSurfaceList_.
    setSingleSurfaceGeometry( cone_, 1 );

    // Create rear sphere ( "end cap" ), fully revolved.
    SphereSegment* rearSphere_ = new SphereSegment( );
    rearSphere_->setMinimumAzimuthAngle( 0.0 );
    rearSphere_->setMaximumAzimuthAngle( 2.0 * PI );

    // Calculate end radius of cone.
    double endRadius_ = cone_->getStartRadius( ) + rearLength_ * tan( rearAngle_ );

    // Calculate rear sphere radius.
    double rearNoseRadius_ = endRadius_ / cos( -rearAngle_ );

    // Set extent of rear sphere.
    rearSphere_->setMinimumZenithAngle( PI / 2.0 - rearAngle_ );
    rearSphere_->setMaximumZenithAngle( PI );

    // Set rear sphere radius.
    rearSphere_->setRadius( rearNoseRadius_ );

    // Set translation vector of rear sphere.
    translationVector_( 0 ) =  ( rearNoseRadius_ * sin( -rearAngle_ ) ) - rearLength_
            - ( sideRadius_ * ( sin( PI / 2.0 - noseSphereAngle_ ) + sin ( -rearAngle_ ) ) );
    rearSphere_->setOffset( translationVector_ );
    setSingleSurfaceGeometry( rearSphere_, 2 );

    // Create torus section of capsule.
    Torus* torus_ = new Torus( );
    torus_->setMinimumIndependentVariable( 1, 0.0 );
    torus_->setMaximumIndependentVariable( 1, 2.0 * PI );
    torus_->setMinimumIndependentVariable( 2, PI / 2.0 - noseSphereAngle_ );
    torus_->setMaximumIndependentVariable( 2, rearAngle_ );

    // Calculate and set major torus radius.
    double torusMajorRadius_ = ( noseRadius_ - sideRadius_ ) * sin( noseSphereAngle_ );
    torus_->setParameter( 0, torusMajorRadius_ );
    torus_->setParameter( 1, sideRadius_ );

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
    rotationMatrix( 2, 1 ) = 0;
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
    stream << "This is a capsule." << endl;
    stream << "The defining parameters are: "<< endl
           << "Nose radius: " << capsule.getNoseRadius( ) << endl
           << "Mid radius: " << capsule.getMiddleRadius( ) << endl
           << "Rear length: " << capsule.getRearLength( ) << endl
           << "Rear angle: " << capsule.getRearAngle( ) << endl
           << "Side radius: " << capsule.getSideRadius( )<< endl;

    // Return stream.
    return stream;
}

} // namespace tudat
