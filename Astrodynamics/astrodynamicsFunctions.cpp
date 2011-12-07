/*! \file astrodynamicsFunctions.cpp
 *    This source file contains general astrodynamics functions.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 2
 *    Check status      : Unchecked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 11 November, 2011
 *    Last modified     : 11 November, 2011
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
 *      100906    K. Kumar          First creation of code.
 *      111115    K. Kumar          Added checker info.
 */

// Include statements.
#include <cmath>
#include <limits>
#include "Astrodynamics/astrodynamicsFunctions.h"
#include "Astrodynamics/physicalConstants.h"

//! Tudat library namespace.
namespace tudat
{

//! Astrodynamics namespace.
namespace astrodynamics
{

//! Compute two-body orbital period.
double computeTwoBodyOrbitalPeriod( double semiMajorAxis,
                                    double gravitationalParameterOfCentralBody,
                                    double massOfOrbitingBody )
{
    return 2.0 * M_PI * std::sqrt( std::pow( semiMajorAxis, 3.0 )
                                   /  ( ( PhysicalConstants::GRAVITATIONAL_CONSTANT
                                          * massOfOrbitingBody )
                                        + gravitationalParameterOfCentralBody ) );
}

//! Compute two-body angular momentum.
double computeTwoBodyAngularMomentum( double semiMajorAxis, double eccentricity,
                                      double gravitationalParameterOfCentralBody,
                                      double massOfOrbitingBody )
{
    return massOfOrbitingBody * std::sqrt( gravitationalParameterOfCentralBody * semiMajorAxis
                                           * ( 1.0 - std::pow( eccentricity, 2.0 ) ) );
}

//! Compute two-body mean motion.
double computeTwoBodyMeanMotion( double semiMajorAxis, double gravitationalParameterOfCentralBody,
                                 double massOfOrbitingBody )
{
    return std::sqrt( ( ( PhysicalConstants::GRAVITATIONAL_CONSTANT * massOfOrbitingBody )
                      + gravitationalParameterOfCentralBody ) / std::pow( semiMajorAxis, 3.0 ) );
}

//! Compute two-body orbital energy.
double computeTwoBodyOrbitalEnergy( double semiMajorAxis,
                                    double gravitationalParameterOfCentralBody,
                                    double massOfOrbitingBody )
{  return -massOfOrbitingBody * gravitationalParameterOfCentralBody / ( 2.0 * semiMajorAxis ); }

//! Compute synodic period.
double computeSynodicPeriod( double orbitalPeriodBody1, double orbitalPeriodBody2 )
{ return 1.0 / std::fabs( 1.0 / orbitalPeriodBody1 - 1.0 / orbitalPeriodBody2 ); }

}

}

// End of file.
