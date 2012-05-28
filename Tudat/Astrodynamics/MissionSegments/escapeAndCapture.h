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
 *      110131    E. Iorfida        Added comments and pointerToCelestialBody.
 *      110202    E. Iorfida        Modified structure of the code, unique base class for launch
 *                                  and capture paths.
 *      110206    E. Iorfida        Modified some comments and name of base class to
 *                                  EscapeAndCapture.
 *      110208    E. Iorfida        Deleted inheritance from TrajectoryDesignMethod, and
 *                                  execute( ), function too. Modified getDeltaV into
 *                                  computeDeltaV.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius, replaced by an element of
 *                                  GeometricShapes.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *
 *    At the moment the shape of the central body is a sphere segment, and the radius of the planet
 *    is set externally by the user. In the future it should be possible to get the radius of each
 *    planet directly from the CelestialBody class, by a link to GeometricShape class.
 */

#ifndef TUDAT_ESCAPE_AND_CAPTURE_H
#define TUDAT_ESCAPE_AND_CAPTURE_H

#include <boost/shared_ptr.hpp>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Bodies/celestialBody.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"

namespace tudat
{
namespace astrodynamics
{
namespace mission_segments
{

//! Escape and capture base class.
/*!
 * Escape and capture class.
 */
class EscapeAndCapture
{
public:

    //! Typedef for shared pointer to gravity field model.
    /*!
     * Typedef for shared pointer to gravity field model.
     */
    typedef boost::shared_ptr< astrodynamics::gravitation::GravityFieldModel >
    GravityFieldModelPointer;

    //! Typedef for shared pointer to sphere segment.
    /*!
     * Typedef for shared pointer to sphere segment.
     */
    typedef boost::shared_ptr< mathematics::geometric_shapes::SphereSegment > SphereSegmentPointer;

    //! Default constructor.
    /*!
     * Default constructor.
     */
    EscapeAndCapture( )
        : semiMajorAxis_( TUDAT_NAN ),
          eccentricity_( TUDAT_NAN ),
          periapsisAltitude_( TUDAT_NAN ),
          apoapsisAltitude_( TUDAT_NAN ),
          hyperbolicExcessSpeed_( TUDAT_NAN ),
          deltaV_ ( TUDAT_NAN ),
          parkingOrbitRadius_( TUDAT_NAN )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~EscapeAndCapture( ) { }

    //! Set central gravity field for the swing-by.
    /*!
     * Sets pointer to central body of the swing-by.
     * \param gravityField Central body of the swing-by.
     */
    void setCentralGravityField( GravityFieldModelPointer gravityField )
    {
        centralBodyGravityfield_ = gravityField;
    }

    //! Set semi-major axis of parking orbit.
    /*!
     * Sets semi-major axis of parking orbit.
     * \param semiMajorAxis Semi-major axis of parking orbit.
     */
    void setSemiMajorAxis( double semiMajorAxis ) { semiMajorAxis_ = semiMajorAxis; }

    //! Set eccentricity of parking orbit.
    /*!
     * Sets eccentricity of parking orbit.
     * \param eccentricity Eccentricity of parking orbit.
     */
    void setEccentricity( double eccentricity ) { eccentricity_ = eccentricity; }

    //! Set periapsis altitude of parking orbit.
    /*!
     * Sets periapsis altitude of parking orbit.
     * \param periapsisAltitude Periapsis altitude of parking orbit.
     */
    void setPeriapsisAltitude( double periapsisAltitude )
    {
        periapsisAltitude_ = periapsisAltitude;
    }

    //! Set apoapsis altitude of parking orbit.
    /*!
     * Sets apoapsis altitude of parking orbit.
     * \param apoapsisAltitude Apoapsis altitude of parking orbit.
     */
    void setApoapsisAltitude( double apoapsisAltitude )
    {
        apoapsisAltitude_ = apoapsisAltitude;
    }

    //! Set radius of the parking orbit.
    /*!
     * \param parkingOrbitRadius Set the radius of the departing/arriving parking orbit.
     */
    void setParkingOrbitRadius( double parkingOrbitRadius )
    {
        parkingOrbitRadius_ = parkingOrbitRadius;
    }

    //! Set hyperbolic excess speed at launch/capture phase.
    /*!
     * Sets hyperbolic excess speed at launch/capture phase.
     * \param hyperbolicExcessSpeed Hyperbolic excess speed at
     *            launch/capture phase in a parking orbit.
     */
    void setHyperbolicExcessSpeed( double hyperbolicExcessSpeed )
    {
        hyperbolicExcessSpeed_ = hyperbolicExcessSpeed;
    }

    //! Compute delta-V of launch/capture phase.
    /*!
     * Get Delta-V to get from the parking orbit to a hyperbolic escape
     * trajectory with the desired hyperbolic excess velocity.
     * \return Delta-V of parking orbit.
     */
    double computeDeltaV( );

protected:

    //! The gravity field produced by the CelestialBody.
    /*!
     * The gravity field in which the swing-by is performed.
     */
    GravityFieldModelPointer centralBodyGravityfield_;

    //! Semi-major axis of parking orbit.
    /*!
     * Semi-major axis of parking orbit.
     */
    double semiMajorAxis_;

    //! Eccentricity of parking orbit.
    /*!
     * Eccentricity of parking orbit.
     */
    double eccentricity_;

    //! Periapsis altitude of parking orbit.
    /*!
     * Periapsis altitude of parking orbit.
     */
    double periapsisAltitude_;

    //! Apoapsis altitude of parking orbit.
    /*!
     * Apoapsis altitude of parking orbit.
     */
    double apoapsisAltitude_;

    //! Hyperbolic excess speed at launch/capture phase in a parking orbit.
    /*!
     * Hyperbolic excess speed at launch/capture phase in a parking orbit.
     */
    double hyperbolicExcessSpeed_;

    //! Delta-V of parking orbit.
    /*!
     * Delta-V to get from parking orbit to escape speed.
     */
    double deltaV_;

    //! Radius of the parking orbit.
    /*!
     * Radius of the parking orbit.
     */
    double parkingOrbitRadius_;

    //! Shared pointer to SphereSegment class for central body.
    /*!
     * Shared pointer to SphereSegment class for central body.
     */
    SphereSegmentPointer centralBodySphere_;

private:
};

} // namespace mission_segments
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_ESCAPE_AND_CAPTURE_H
