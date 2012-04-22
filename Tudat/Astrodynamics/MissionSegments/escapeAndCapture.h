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
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// At the moment the shape of the central body is a sphere segment,
// and the radius of the planet is set externally by the user.
// In the future it should be possible to get the radius of each planet
// directly from the CelestialBody class, by a link to GeometricShape
// class.
// 

#ifndef TUDAT_ESCAPE_AND_CAPTURE_H
#define TUDAT_ESCAPE_AND_CAPTURE_H

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

    //! Default constructor.
    /*!
     * Default constructor.
     */
    EscapeAndCapture( );

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
    void setCentralGravityField( GravityFieldModel* gravityField )
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
    GravityFieldModel* centralBodyGravityfield_;

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

    //! Pointer to SphereSegment class for central body.
    /*!
     * Pointer to SphereSegment class for central body.
     */
    SphereSegment* pointerToCentralBodySphere_;

private:
};

} // namespace mission_segments
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_ESCAPE_AND_CAPTURE_H
