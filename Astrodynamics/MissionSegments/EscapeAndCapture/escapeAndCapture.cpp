/*! \file escapeAndCapture.cpp
 *    This header file contains a base class for the implementation of required
 *    delta-V to escape from a parking orbit around the initial body into
 *    the interplanetary trajectory, and of the required delta-V needed to be
 *    captured in a parking orbit around the final body at the end of the
 *    interplanetary transfer trajectory.
 *
 *    Path              : /Astrodynamics/MissionSegments/EscapeAndCapture/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 29 January, 2011
 *    Last modified     : 14 February, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      110129    E. Iorfida        First creation of code.
 *      110131    E. Iorfida        Added pointerToCentralBody.
 *      110202    E. Iorfida        Modified structure of the code, unique
 *                                  base class for launch and capture paths.
 *      110208    E. Iorfida        Deleted inheritance from
 *                                  TrajectoryDesignMethod, and execute()
 *                                  function, too. Modified getDeltaV into
 *                                  computeDeltaV.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius,
 *                                  replaced by an element of GeometricShapes.
 */

// Include statements.
#include "escapeAndCapture.h"
#include "sphereSegment.h"

// Include directives.
using mathematics::raiseToIntegerPower;

//! Default constructor.
EscapeAndCapture::EscapeAndCapture( ) :
        semiMajorAxis_ ( -0.0 ),
        eccentricity_ ( -1.0 ),
        periapsisAltitude_ ( -0.0 ),
        apoapsisAltitude_( -0.0 ),
        hyperbolicExcessSpeed_( -1.0 ),
        deltaV_ ( -0.0 )

{
}

//! Default destructor.
EscapeAndCapture::~EscapeAndCapture( )
{
}

//! Set semi-major axis of parking orbit.
void EscapeAndCapture::setSemiMajorAxis( const double &semiMajorAxis )
{
    semiMajorAxis_ = semiMajorAxis;
}

//! Set eccentricity of parking orbit.
void EscapeAndCapture::setEccentricity( const double &eccentricity )
{
    eccentricity_ = eccentricity;
}

//! Set periapsis altitude of parking orbit.
void EscapeAndCapture::setPeriapsisAltitude(
        const double &periapsisAltitude )
{
    periapsisAltitude_ = periapsisAltitude;
}

//! Set apoapsis altitude of parking orbit.
void EscapeAndCapture::setApoapsisAltitude( const double &apoapsisAltitude )
{
    apoapsisAltitude_ = apoapsisAltitude;
}

//! Set central body of parking orbit.
void EscapeAndCapture::setCentralBody( CelestialBody *pointerToCentralBody )
{
    pointerToCentralBody_ = pointerToCentralBody;
}

//! Set hyperbolic excess speed at launch/capture phase.
void EscapeAndCapture::setHyperbolicExcessSpeed(
        const double &hyperbolicExcessSpeed )
{
    hyperbolicExcessSpeed_ = hyperbolicExcessSpeed;
}

//! Compute delta-V of launch/capture phase.
double& EscapeAndCapture::computeDeltaV( )
{
    // Local variables.
    double periapsisRadius_;
    double apoapsisRadius_;

    // Get shape model of central body.
    pointerToCentralBodySphere_ = static_cast< SphereSegment* >
                                  ( pointerToCentralBody_->getShapeModel( ) );

    // Check input parameters.
    // For the correct functioning of this routine, the periapsis
    // radius and the eccentricity will have to be known.
    // Eccentricity and semi-major axis set by user.
    if ( eccentricity_ != ( -1.0 ) && semiMajorAxis_ != ( -0.0 ) )
    {
        // Compute periapsis radius.
        periapsisRadius_ = semiMajorAxis_ * ( 1.0 - eccentricity_ );

    }
    // Periapsis and apoapsis altitudes set by user.
    else if ( periapsisAltitude_ != ( -0.0 ) &&
              apoapsisAltitude_ != ( -0.0 ) )
    {
        // Compute periapsis and apoapsis radii.
        periapsisRadius_ = periapsisAltitude_ +
                           pointerToCentralBodySphere_->getRadius( );
        apoapsisRadius_ = apoapsisAltitude_ +
                          pointerToCentralBodySphere_->getRadius( );

        // Compute eccentricity.
        eccentricity_ = ( apoapsisRadius_ - periapsisRadius_ ) /
                        ( apoapsisRadius_ + periapsisRadius_ ) ;

    }
    // Periapsis altitude and eccentricy set by user.
    else if ( periapsisAltitude_ != ( -0.0 ) && eccentricity_ != ( -1.0 ) )
    {
        // Compute periapsis radius.
        periapsisRadius_ = periapsisAltitude_ +
                           pointerToCentralBodySphere_->getRadius( );
    }

    // Compute escape velocity squared.
    double escapeVelocitySquared_;
    escapeVelocitySquared_ = 2.0 * pointerToCentralBody_->
                            getGravitationalParameter( ) / periapsisRadius_ ;

    // Compute delta-V.
    deltaV_ = sqrt( escapeVelocitySquared_ + raiseToIntegerPower(
                    hyperbolicExcessSpeed_ , 2 ) ) -
              sqrt( ( escapeVelocitySquared_ / 2.0 ) *
                    ( 1.0 + eccentricity_ ) );

    return deltaV_;
}

// End of file.
