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
 *      110208    E. Iorfida        Deleted inheritance from TrajectoryDesignMethod, and
 *                                  execute( ), function too. Modified getDeltaV into
 *                                  computeDeltaV.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius, replaced by an element of
 *                                  GeometricShapes.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120416    T. Secretin       Corrected if-statements to detect NaN values.
 *
 *    References
 *
 */

#include <boost/shared_ptr.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"

namespace tudat
{
namespace astrodynamics
{
namespace mission_segments
{

//! Compute delta-V of launch/capture phase.
double EscapeAndCapture::computeDeltaV( )
{
    using std::endl;
    using std::pow;
    using std::sqrt;

    // Local variables.
    double periapsisRadius_ = TUDAT_NAN;
    double apoapsisRadius_ = TUDAT_NAN;

    // Check input parameters.
    // For the correct functioning of this routine, the periapsis
    // radius and the eccentricity will have to be known.
    // Eccentricity and semi-major axis set by user.
    if ( !( boost::math::isnan )( eccentricity_ ) && !( boost::math::isnan )( semiMajorAxis_ ) )
    {
        // Compute periapsis radius.
        periapsisRadius_ = semiMajorAxis_ * ( 1.0 - eccentricity_ );

    }
    // Periapsis and apoapsis altitudes set by user.
    else if ( !( boost::math::isnan )( periapsisAltitude_ ) &&
              !( boost::math::isnan )( apoapsisAltitude_ ) )
    {
        // Compute periapsis and apoapsis radii.
        periapsisRadius_ = periapsisAltitude_ + parkingOrbitRadius_;
        apoapsisRadius_ = apoapsisAltitude_ + parkingOrbitRadius_;

        // Compute eccentricity.
        eccentricity_ = ( apoapsisRadius_ - periapsisRadius_ ) /
                        ( apoapsisRadius_ + periapsisRadius_ );

    }
    // Periapsis altitude and eccentricy set by user.
    else if ( !( boost::math::isnan )( periapsisAltitude_ ) &&
              !( boost::math::isnan )( eccentricity_ ) )
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

} // namespace mission_segments
} // namespace astrodynamics
} // namespace tudat
