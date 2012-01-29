/*! \file convertMeanAnomalyToEccentricAnomaly.cpp
 *    This source file contains a class to convert mean anomly to eccentric anomaly for elliptical
 *    orbits. It makes use of the simple iterative method to solve Kepler's equation.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : T. Secretin
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : T.A.LeitePintoSecretin@student.tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : Simon Billemont
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : simon@angelcorp.be
 *
 *    Date created      : 10 February, 2011
 *    Last modified     : 21 December, 2011
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *
 *    Notes
 *      Currently, this conversion is only valid for eccentricities up to 0.97 due to the
 *      difficulties of the iterative method to converge for eccentricities close to 1.
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
 *      111209    T. Secretin       Relaxed constraints on near-parabolic check.
 *      111221    T. Secretin       Added zero eccentricity case and check for negative
 *                                  eccentricities.
 */

// Include statements.
#include <iostream>
#include "Astrodynamics/States/convertMeanAnomalyToEccentricAnomaly.h"

//! Tudat library namespace.
namespace tudat
{

// Using declarations.
using std::cerr;
using std::endl;

//! Orbital element conversions namespace.
namespace orbital_element_conversions
{

//! Convert mean anomaly to eccentric anomaly.
double ConvertMeanAnomalyToEccentricAnomaly::convert( )
{
    // Set the class that contains the functions needed for Newton-Raphson.
    newtonRaphsonAdaptor_.setClass( this );

    // Set NewtonRaphson adaptor class.
    pointerToNewtonRaphson_->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptor_ );

    // Check if orbit is circular.
    if ( eccentricity_ == 0.0 )
    {
        // If orbit is circular mean anomaly and eccentric anomaly are equal.
        eccentricAnomaly_ = meanAnomaly_;
    }
    // Check if orbit is elliptical, and not near-parabolic.
    else if ( eccentricity_ < 0.98 && eccentricity_ > 0.0 )
    {
        // Set mathematical functions.
        newtonRaphsonAdaptor_.setPointerToFunction( &ConvertMeanAnomalyToEccentricAnomaly::
                                                    computeKeplersFunctionForEllipticalOrbits_ );
        newtonRaphsonAdaptor_.setPointerToFirstDerivativeFunction(
                    &ConvertMeanAnomalyToEccentricAnomaly::
                    computeFirstDerivativeKeplersFunctionForEllipticalOrbits_ );

        // Set initial guess of eccentric anomaly to the mean anomaly.
        pointerToNewtonRaphson_->setInitialGuessOfRoot( meanAnomaly_ );

        // Execute Newton-Raphon method.
        pointerToNewtonRaphson_->execute( );

        // Set eccentric anomaly based on result of Newton-Raphson root-finding algorithm.
        eccentricAnomaly_ = pointerToNewtonRaphson_->getComputedRootOfFunction( );
    }
    // Check if orbit is near-parabolic, i.e. eccentricity >= 0.98.
    else
    {
        cerr << "Orbit is near-parabolic and, at present conversion, between eccentric anomaly "
             << "and mean anomaly is not possible for eccentricities larger than: 0.98" << endl;

        // Set eccentric anomaly to error value.
        eccentricAnomaly_ = std::numeric_limits<double>::signaling_NaN() ;
    }

    // Return eccentric anomaly.
    return eccentricAnomaly_;
}

} // Namespace orbital_element_conversions.

} // Namespace tudat.

// End of file.
