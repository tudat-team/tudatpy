/*! \file keplerianElements.cpp
 *    This source file contains the Keplerian elements class included in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 22 October, 2010
 *    Last modified     : 10 March, 2011
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
 *      101022    K. Kumar          First creation of code.
 *      101110    K. Kumar          Added get functions for auxilliary
 *                                  parameters.
 *      101130    E. Iorfida        Added set function for semi-latus rectum.
 *      101202    J. Melman         Only setting the state_ vector now, 6
 *                                  private Keplerian elements not used anymore.
 *      110310    K. Kumar          Changed right ascension of ascending node
 *                                  to longitude of ascending node.
 */

// Include statements.
#include "keplerianElements.h"

//! Default constructor.
KeplerianElements::KeplerianElements( ) : semiLatusRectum_( -0.0 )
{
    // Initialize variables.
    state.setZero( 6 );
}

//! Default destructor.
KeplerianElements::~KeplerianElements( )
{
}

//! Set semi-major axis.
void KeplerianElements::setSemiMajorAxis( const double& semiMajorAxis )
{
    state( 0 ) = semiMajorAxis;
}

//! Set eccentricity.
void KeplerianElements::setEccentricity( const double& eccentricity )
{
    state( 1 ) = eccentricity;
}

//! Set inclination.
void KeplerianElements::setInclination( const double& inclination )
{
    state( 2 ) = inclination;
}

//! Set argument of periapsis.
void KeplerianElements::setArgumentOfPeriapsis( const double& argumentOfPeriapsis )
{
    state( 3 ) = argumentOfPeriapsis;
}

//! Set longitude of ascending node.
void KeplerianElements::setLongitudeOfAscendingNode( const double&
                                                          longitudeOfAscendingNode )
{
    state( 4 ) = longitudeOfAscendingNode;
}

//! Set true anomaly.
void KeplerianElements::setTrueAnomaly( const double& trueAnomaly )
{
    state( 5 ) = trueAnomaly;
}

//! Set semi-latus rectum ( for parabolic orbits ).
void KeplerianElements::setSemiLatusRectum( const double& semiLatusRectum )
{
    semiLatusRectum_ = semiLatusRectum;
}

//! Get semi-major axis.
double& KeplerianElements::getSemiMajorAxis( )
{
    return state( 0 );
}

//! Get eccentricity.
double& KeplerianElements::getEccentricity( )
{
    return state( 1 );
}

//! Get inclination.
double& KeplerianElements::getInclination()
{
    return state( 2 );
}

//! Get argument of periapsis.
double& KeplerianElements::getArgumentOfPeriapsis( )
{
    return state( 3 );
}

//! Get longitude of ascending node.
double& KeplerianElements::getLongitudeOfAscendingNode( )
{
    return state( 4 );
}

//! Get true anomaly.
double& KeplerianElements::getTrueAnomaly( )
{
    return state( 5 );
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
    return getArgumentOfPeriapsis( ) + getLongitudeOfAscendingNode( );
}

//! Get true longitude.
double KeplerianElements::getTrueLongitude( )
{
    return getArgumentOfPeriapsis( ) + getLongitudeOfAscendingNode( )
            + getTrueAnomaly( );
}

//! Get argument of latitude.
double KeplerianElements::getArgumentOfLatitude( )
{
    return getArgumentOfPeriapsis( ) + getTrueAnomaly( );
}

// End of file.
