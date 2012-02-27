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
 *      110129    E. Iorfida        First creation of code.
 *      110131    E. Iorfida        Added pointerToCentralBody.
 *      110202    E. Iorfida        Modified structure of the code, unique base class for launch
 *                                  and capture paths.
 *      110208    E. Iorfida        Deleted inheritance from TrajectoryDesignMethod, and execute( ),
 *                                  function too. Modified getDeltaV into computeDeltaV.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius, replaced by an element of
 *                                  GeometricShapes.
 *
 *    References
 *
 */

#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"

namespace tudat
{

//! Default constructor.
EscapeAndCapture::EscapeAndCapture( ) : 
    centralBodyGravityfield_( NULL ), 
    semiMajorAxis_ ( std::numeric_limits<double>::signaling_NaN( ) ), 
    eccentricity_ ( std::numeric_limits<double>::signaling_NaN( ) ),
    periapsisAltitude_ ( std::numeric_limits<double>::signaling_NaN( ) ), 
    apoapsisAltitude_( std::numeric_limits<double>::signaling_NaN( ) ), 
    hyperbolicExcessSpeed_( std::numeric_limits<double>::signaling_NaN( ) ),
    deltaV_ ( std::numeric_limits<double>::signaling_NaN( ) ), 
    parkingOrbitRadius_( std::numeric_limits<double>::signaling_NaN( ) ),
    pointerToCentralBodySphere_ ( NULL )
{
}

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
        periapsisRadius_ = periapsisAltitude_ + parkingOrbitRadius_;
        apoapsisRadius_ = apoapsisAltitude_ + parkingOrbitRadius_;

        // Compute eccentricity.
        eccentricity_ = ( apoapsisRadius_ - periapsisRadius_ ) /
                        ( apoapsisRadius_ + periapsisRadius_ );

    }
    // Periapsis altitude and eccentricy set by user.
    else if ( periapsisAltitude_ != ( -0.0 ) && eccentricity_ != ( -1.0 ) )
    {
        // Compute periapsis radius.
        periapsisRadius_ = periapsisAltitude_ + parkingOrbitRadius_;
    }

    // Compute escape velocity squared.
    double escapeVelocitySquared_;
    escapeVelocitySquared_ = 2.0 * centralBodyGravityfield_->getGravitationalParameter( )
            / periapsisRadius_;

    // Compute delta-V.
    deltaV_ = sqrt( escapeVelocitySquared_ + pow( hyperbolicExcessSpeed_, 2.0 ) ) -
              sqrt( ( escapeVelocitySquared_ / 2.0 ) * ( 1.0 + eccentricity_ ) );

    return deltaV_;
}

} // namespace tudat
