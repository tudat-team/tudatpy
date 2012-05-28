/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      110129    E. Iorfida        Creation of code.
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
