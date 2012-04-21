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
 *      110210    K. Kumar          First creation of code.
 *      111209    T. Secretin       Relaxed constraints on near-parabolic check.
 *      111221    T. Secretin       Added zero eccentricity case and check for negative
 *                                  eccentricities.
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *
 */

// Temporary notes (move to class/function doxygen):
// Currently, this conversion is only valid for eccentricities up to 0.97 due to the
// difficulties of the iterative method to converge for eccentricities close to 1.
// 

#include <iostream>
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"

namespace tudat
{
namespace orbital_element_conversions
{

// Using declarations.
using std::cerr;
using std::endl;

//! Convert mean anomaly to eccentric anomaly.
double ConvertMeanAnomalyToEccentricAnomaly::convert( )
{
    // Declare eccentric anomaly.
    double eccentricAnomaly_ = std::numeric_limits<double>::signaling_NaN( );

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
    }

    // Return eccentric anomaly.
    return eccentricAnomaly_;
}

} // namespace orbital_element_conversions
} // namespace tudat
