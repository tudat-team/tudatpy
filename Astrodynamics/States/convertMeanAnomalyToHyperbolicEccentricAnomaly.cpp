/*! \file convertMeanAnomalyToHyperbolicEccentricAnomaly.cpp
 *    This source file contains a class to convert mean anomaly to hyperbolic
 *    eccentric anomaly for hyperbolic orbits.
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
 *    Date created      : 14 February, 2011
 *    Last modified     : 14 February, 2011
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series,
 *          VA, 2002.
 *      http://www.cdeagle.com/omnum/pdf/demokep1.pdf, last accessed:
 *          16th February, 2011.
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
 *      110214    K. Kumar          First creation of code.
 */

// Include statements.
#include <iostream>
#include "convertMeanAnomalyToHyperbolicEccentricAnomaly.h"

// Using declarations.
using std::cerr;
using std::endl;

//! Orbital element conversions namespace.
namespace orbital_element_conversions
{

//! Default constructor.
ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
        ConvertMeanAnomalyToHyperbolicEccentricAnomaly( )
            : hyperbolicEccentricAnomaly_( -1.0 )
{
}

//! Default destructor.
ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
        ~ConvertMeanAnomalyToHyperbolicEccentricAnomaly( )
{
}

//! Convert mean anomaly to hyperbolic eccentric anomaly.
double ConvertMeanAnomalyToHyperbolicEccentricAnomaly::convert( )
{
    // Set the class that contains the functions needed for Newton-Raphson.
    newtonRaphsonAdaptor_.setClass( this );

    // Set NewtonRaphson adaptor class.
    pointerToNewtonRaphson_->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptor_ );

    // Check if orbit is hyperbolic, and not near-parabolic.
    if ( eccentricity_ > 1.2 )
    {
        // Set mathematical functions.
        newtonRaphsonAdaptor_.setPointerToFunction(
                &ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
                computeKeplersFunctionForHyperbolicOrbits_ );
        newtonRaphsonAdaptor_
                .setPointerToFirstDerivativeFunction(
                    &ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
                   computeFirstDerivativeKeplersFunctionForHyperbolicOrbits_ );

        // Set initial guess of hyperbolic eccentric anomaly to the mean anomaly.
        pointerToNewtonRaphson_
                ->setInitialGuessOfRoot( 2.0 * meanAnomaly_ / eccentricity_
                                         - 1.8 );

        // Execute Newton-Raphon method.
        pointerToNewtonRaphson_->execute( );

        // Set hyperbolic eccentric anomaly based on result of Newton-Raphson
        // root-finding algorithm.
        hyperbolicEccentricAnomaly_ = pointerToNewtonRaphson_
                                      ->getComputedRootOfFunction( );
    }

    // Check if orbit is near-parabolic.
    else if ( eccentricity_ <= 1.2 )
    {
        cerr << "Orbit is near-parabolic and, at present conversion, "
             << "between hyperbolic eccentric anomaly and hyperbolic mean "
             << "anomaly is not possible for eccentricities in the range: "
             << "0.8 < eccentricity < 1.2." << endl;

        // Set hyperbolic eccentric anomaly to error value.
        hyperbolicEccentricAnomaly_ = -1.0;
    }

    // Return hyperbolic eccentric anomaly.
    return hyperbolicEccentricAnomaly_;
}

//! Compute Kepler's function for hyperbolic orbits.
double ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
        computeKeplersFunctionForHyperbolicOrbits_(
                double& hyperbolicEccentricAnomaly )
{
    // Return value of Kepler's function for hyperbolic orbits.
    return eccentricity_ * sinh( hyperbolicEccentricAnomaly )
            - hyperbolicEccentricAnomaly - meanAnomaly_;
}

//! Compute first-derivative of Kepler's function for hyperbolic orbits.
double ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
        computeFirstDerivativeKeplersFunctionForHyperbolicOrbits_(
                double& hyperbolicEccentricAnomaly )
{
    // Return value of first-derivative of Kepler's function for
    // hyperbolic orbits.
    return eccentricity_ * cosh( hyperbolicEccentricAnomaly ) - 1.0;
}

}

// End of file.
