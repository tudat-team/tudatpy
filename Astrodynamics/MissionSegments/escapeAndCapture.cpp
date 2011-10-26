/*! \file escapeAndCapture.cpp
 *    This source file contains a base class for the implementation of required delta-V to escape
 *    from a parking orbit around the initial body into the interplanetary trajectory, and of the
 *    required delta-V needed to be  captured in a parking orbit around the final body at the end
 *    of the interplanetary transfer trajectory.
 *
 *    Path              : /Astrodynamics/MissionSegments/
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
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110202    E. Iorfida        Modified structure of the code, unique base class for launch
 *                                  and capture paths.
 *      110208    E. Iorfida        Deleted inheritance from TrajectoryDesignMethod, and execute(),
 *                                  function too. Modified getDeltaV into computeDeltaV.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius, replaced by an element of
 *                                  GeometricShapes.
 */

// Include statements.
#include "Astrodynamics/MissionSegments/escapeAndCapture.h"
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/LinearAlgebra/linearAlgebra.h"

//! Compute delta-V of launch/capture phase.
double& EscapeAndCapture::computeDeltaV( )
{
    // Using declarations.
    using std::endl;
    using std::pow;
    using std::sqrt;

    // Local variables.
    double periapsisRadius_ = -0.0;
    double apoapsisRadius_ = -0.0;

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
        periapsisRadius_ = periapsisAltitude_ + pointerToCentralBodySphere_->getRadius( );
        apoapsisRadius_ = apoapsisAltitude_ + pointerToCentralBodySphere_->getRadius( );

        // Compute eccentricity.
        eccentricity_ = ( apoapsisRadius_ - periapsisRadius_ ) /
                        ( apoapsisRadius_ + periapsisRadius_ );

    }
    // Periapsis altitude and eccentricy set by user.
    else if ( periapsisAltitude_ != ( -0.0 ) && eccentricity_ != ( -1.0 ) )
    {
        // Compute periapsis radius.
        periapsisRadius_ = periapsisAltitude_ + pointerToCentralBodySphere_->getRadius( );
    }

    // Compute escape velocity squared.
    double escapeVelocitySquared_;
    escapeVelocitySquared_ = 2.0 * pointerToCentralBody_->getGravitationalParameter( )
            / periapsisRadius_;

    // Compute delta-V.
    deltaV_ = sqrt( escapeVelocitySquared_ + pow( hyperbolicExcessSpeed_, 2.0 ) ) -
              sqrt( ( escapeVelocitySquared_ / 2.0 ) * ( 1.0 + eccentricity_ ) );

    return deltaV_;
}

// End of file.
