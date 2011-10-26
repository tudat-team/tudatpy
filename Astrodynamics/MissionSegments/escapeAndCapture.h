/*! \file escapeAndCapture.h
 *    This header file contains a base class for the implementation of required delta-V to escape
 *    from a parking orbit around the initial body into the interplanetary trajectory, and of the
 *    required delta-V needed to be  captured in a parking orbit around the final body at the end
 *    of the interplanetary transfer trajectory.
 *
 *    Path              : /Astrodynamics/MissionSegments/
 *    Version           : 6
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
 *    At the moment the shape of the central body is a sphere segment,
 *    and the radius of the planet is set externally by the user.
 *    In the future it should be possible to get the radius of each planet
 *    directly from the CelestialBody class, by a link to GeometricShape
 *    class.
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
 *      110131    E. Iorfida        Added comments and pointerToCelestialBody.
 *      110202    E. Iorfida        Modified structure of the code, unique base class for launch
 *                                  and capture paths.
 *      110206    E. Iorfida        Modified some comments and name of base class to
 *                                  EscapeAndCapture.
 *      110208    E. Iorfida        Deleted inheritance from TrajectoryDesignMethod, and execute(),
 *                                  function too. Modified getDeltaV into computeDeltaV.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius, replaced by an element of
 *                                  GeometricShapes.
 */

#ifndef ESCAPEANDCAPTURE_H
#define ESCAPEANDCAPTURE_H

// Include statements.
#include <iostream>
#include "Astrodynamics/Bodies/celestialBody.h"
#include "Mathematics/GeometricShapes/sphereSegment.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
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
    EscapeAndCapture( ) : semiMajorAxis_ ( -0.0 ), eccentricity_ ( -1.0 ),
        periapsisAltitude_ ( -0.0 ), apoapsisAltitude_( -0.0 ), hyperbolicExcessSpeed_( -1.0 ),
        deltaV_ ( -0.0 ), pointerToCentralBody_( NULL ), pointerToCentralBodySphere_ ( NULL ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~EscapeAndCapture( ) { }

    //! Set semi-major axis of parking orbit.
    /*!
     * Sets semi-major axis of parking orbit.
     * \param semiMajorAxis Semi-major axis of parking orbit.
     */
    void setSemiMajorAxis( const double &semiMajorAxis ) { semiMajorAxis_ = semiMajorAxis; }

    //! Set eccentricity of parking orbit.
    /*!
     * Sets eccentricity of parking orbit.
     * \param eccentricity Eccentricity of parking orbit.
     */
    void setEccentricity( const double &eccentricity ) { eccentricity_ = eccentricity; }

    //! Set periapsis altitude of parking orbit.
    /*!
     * Sets periapsis altitude of parking orbit.
     * \param periapsisAltitude Periapsis altitude of parking orbit.
     */
    void setPeriapsisAltitude( const double &periapsisAltitude )
    { periapsisAltitude_ = periapsisAltitude; }

    //! Set apoapsis altitude of parking orbit.
    /*!
     * Sets apoapsis altitude of parking orbit.
     * \param apoapsisAltitude Apoapsis altitude of parking orbit.
     */
    void setApoapsisAltitude( const double &apoapsisAltitude )
    { apoapsisAltitude_ = apoapsisAltitude; }

    //! Set central body of parking orbit.
    /*!
     * Sets pointer to central body of parking orbit.
     * \param pointerToCentralBody Central body of parking orbit.
     */
    void setCentralBody( CelestialBody *pointerToCentralBody )
    { pointerToCentralBody_ = pointerToCentralBody; }

    //! Set hyperbolic excess speed at launch/capture phase.
    /*!
     * Sets hyperbolic excess speed at launch/capture phase.
     * \param hyperbolicExcessSpeed Hyperbolic excess speed at
     *            launch/capture phase in a parking orbit.
     */
    void setHyperbolicExcessSpeed( const double &hyperbolicExcessSpeed )
    { hyperbolicExcessSpeed_ = hyperbolicExcessSpeed; }

    //! Compute delta-V of launch/capture phase.
    /*!
     * Get Delta-V to get from the parking orbit to a hyperbolic escape
     * trajectory with the desired hyperbolic excess velocity.
     * \return Delta-V of parking orbit.
     */
    double& computeDeltaV( );

protected:

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

    //! Pointer to CelestialBody class for parking orbit.
    /*!
     * Pointer to CelestialBody class for parking orbit.
     */
    CelestialBody* pointerToCentralBody_;

    //! Pointer to SphereSegment class for central body.
    /*!
     * Pointer to SphereSegment class for central body.
     */
    SphereSegment* pointerToCentralBodySphere_;

private:
};

}

#endif // ESCAPEANDCAPTURE_H

// End of file.
