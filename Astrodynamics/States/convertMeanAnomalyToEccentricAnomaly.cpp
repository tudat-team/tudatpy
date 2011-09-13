/*! \file convertMeanAnomalyToEccentricAnomaly.cpp
 *    This source file contains a class to convert mean anomly to eccentric
 *    anomaly for elliptical and hyperbolic orbits.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 10 February, 2011
 *    Last modified     : 10 February, 2011
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series,
 *          VA, 2002.
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
 *      110210    K. Kumar          First creation of code.
 */

// Include statement.
#include "convertMeanAnomalyToEccentricAnomaly.h"

// Using declarations.
using mathematics::raiseToIntegerPower;

//! Orbital element conversions namespace.
namespace orbital_element_conversions
{

//! Default constructor.
ConvertMeanAnomalyToEccentricAnomaly::ConvertMeanAnomalyToEccentricAnomaly( )
    : eccentricAnomaly_( -1.0 )
{
}

//! Default destructor.
ConvertMeanAnomalyToEccentricAnomaly::~ConvertMeanAnomalyToEccentricAnomaly( )
{
}

//! Convert mean anomaly to eccentric anomaly.
double ConvertMeanAnomalyToEccentricAnomaly::convert( )
{
    // Set the class that contains the functions needed for Newton-Raphson.
    newtonRaphsonAdaptor_.setClass( this );

    // Set NewtonRaphson adaptor class.
    pointerToNewtonRaphson_->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptor_ );

    // Check if orbit is elliptical, and not near-parabolic.
    if ( eccentricity_ < 0.8 )
    {
        // Set mathematical functions.
        newtonRaphsonAdaptor_.setPointerToFunction(
                &ConvertMeanAnomalyToEccentricAnomaly::
                computeKeplersFunctionForEllipticalOrbits_ );
        newtonRaphsonAdaptor_
                .setPointerToFirstDerivativeFunction(
                        &ConvertMeanAnomalyToEccentricAnomaly::
                   computeFirstDerivativeKeplersFunctionForEllipticalOrbits_ );

        // Set initial guess of eccentric anomaly to the mean anomaly.
        pointerToNewtonRaphson_->setInitialGuessOfRoot( meanAnomaly_ );

        // Execute Newton-Raphon method.
        pointerToNewtonRaphson_->execute( );

        // Set eccentric anomaly based on result of Newton-Raphson
        // root-finding algorithm.
        eccentricAnomaly_ = pointerToNewtonRaphson_
                            ->getComputedRootOfFunction( );
    }

    // Check if orbit is near-parabolic.
    else if ( eccentricity_ >= 0.8 )
    {
        cerr << "Orbit is near-parabolic and, at present conversion, "
             << "between eccentric anomaly and mean anomaly is not "
             << "possible for eccentricities in the range: "
             << "0.8 < eccentricity < 1.2." << endl;

        // Set eccentric anomaly to error value.
        eccentricAnomaly_ = -1.0;
    }

    // Return eccentric anomaly.
    return eccentricAnomaly_;
}

//! Compute Kepler's function for elliptical orbits.
double ConvertMeanAnomalyToEccentricAnomaly::
        computeKeplersFunctionForEllipticalOrbits_( double& eccentricAnomaly )
{
    // Return value of Kepler's function for elliptical orbits.
    return eccentricAnomaly - eccentricity_ * sin( eccentricAnomaly )
            - meanAnomaly_;
}

//! Compute first-derivative of Kepler's function for elliptical orbits.
double ConvertMeanAnomalyToEccentricAnomaly::
        computeFirstDerivativeKeplersFunctionForEllipticalOrbits_(
                double& eccentricAnomaly )
{
    // Return value of first-derivative of Kepler's function for
    // elliptical orbits.
    return 1.0 - eccentricity_ * cos( eccentricAnomaly );
}

}

// End of file.
