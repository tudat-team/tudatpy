/*! \file keplerianElements.cpp
 *    This source file contains the Keplerian elements class included in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 22 Checked, 2010
 *    Last modified     : 02 December, 2010
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
 *      YYMMDD    author        comment
 *      101022    K. Kumar      First creation of code.
 *      101110    K. Kumar      Added get functions for auxilliary parameters.
 *      101130    E. Iorfida    Added set function for semi-latus rectum.
 *      101202    J. Melman     Only setting the state_ vector now, 6 private
 *                              Keplerian elements not used anymore.
 */

// Include statements.
#include "keplerianElements.h"

//! Default constructor.
KeplerianElements::KeplerianElements( )
{
}

//! Default destructor.
KeplerianElements::~KeplerianElements( )
{
}

//! Set semi-major axis.
void KeplerianElements::setSemiMajorAxis( const double& semiMajorAxis )
{
    state_( 0 ) = semiMajorAxis;
}

//! Set eccentricity.
void KeplerianElements::setEccentricity( const double& eccentricity )
{
    state_( 1 ) = eccentricity;
}

//! Set inclination.
void KeplerianElements::setInclination( const double& inclination )
{
    state_( 2 ) = inclination;
}

//! Set argument of periapsis.
void KeplerianElements::setArgumentOfPeriapsis( const double& argumentOfPeriapsis )
{
    state_( 3 ) = argumentOfPeriapsis;
}

//! Set right ascension of ascending node.
void KeplerianElements::setRightAscensionOfAscendingNode( const double&
                                                          rightAscensionOfAscendingNode )
{
    state_( 4 ) = rightAscensionOfAscendingNode;
}

//! Set true anomaly.
void KeplerianElements::setTrueAnomaly( const double& trueAnomaly )
{
    state_( 5 ) = trueAnomaly;
}

//! Set semi-latus rectum.
void KeplerianElements::setSemiLatusRectum(const double &semiLatusRectum)
{
    semiLatusRectum_ = semiLatusRectum;
}

//! Get semi-major axis.
double& KeplerianElements::getSemiMajorAxis( )
{
    return state_( 0 );
}

//! Get eccentricity.
double& KeplerianElements::getEccentricity( )
{
    return state_( 1 );
}

//! Get inclination.
double& KeplerianElements::getInclination()
{
    return state_( 2 );
}

//! Get argument of periapsis.
double& KeplerianElements::getArgumentOfPeriapsis( )
{
    return state_( 3 );
}

//! Get right ascension of ascending node.
double& KeplerianElements::getRightAscensionOfAscendingNode( )
{
    return state_( 4 );
}

//! Get true anomaly.
double& KeplerianElements::getTrueAnomaly( )
{
    return state_( 5 );
}

//! Get semi-latus rectum.
double KeplerianElements::getSemiLatusRectum( )
{
    // This only works if it has been set, which is necessary for a parabola.
    // It does not compute the semi-latus rectum from the semi-major axis and
    // the eccentricity, since in the case of a parabola, the semi-major axis
    // is not defined.
    return semiLatusRectum_;
}

//! Get longitude of periapsis.
double KeplerianElements::getLongitudeOfPeriapsis( )
{
    return state_( 3 ) + state_( 4 );
}

//! Get true longitude.
double KeplerianElements::getTrueLongitude( )
{
    return state_( 3 ) + state_( 4 ) + state_( 5 );
}

//! Get argument of latitude.
double KeplerianElements::getArgumentOfLatitude( )
{
    return state_( 3 ) + state_( 5 );
}

// End of file.
