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
 *      110214    K. Kumar          First creation of code.
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *      http://www.cdeagle.com/omnum/pdf/demokep1.pdf, last accessed: 16th February, 2011.
 *
 */

#include <iostream>
#include <limits>

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"

namespace tudat
{
namespace orbital_element_conversions
{

// Using declarations.
using std::cerr;
using std::endl;

//! Convert mean anomaly to hyperbolic eccentric anomaly.
double ConvertMeanAnomalyToHyperbolicEccentricAnomaly::convert( )
{
    // Declare hyperbolic eccentric anomaly.
    double hyperbolicEccentricAnomaly_ = std::numeric_limits< double >::signaling_NaN( );

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
        newtonRaphsonAdaptor_.setPointerToFirstDerivativeFunction(
                    &ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
                    computeFirstDerivativeKeplersFunctionForHyperbolicOrbits_ );

        // Set initial guess of hyperbolic eccentric anomaly to the mean anomaly.
        pointerToNewtonRaphson_->setInitialGuessOfRoot(
                    2.0 * hyperbolicMeanAnomaly_ / eccentricity_ - 1.8 );

        // Execute Newton-Raphon method.
        pointerToNewtonRaphson_->execute( );

        // Set hyperbolic eccentric anomaly based on result of Newton-Raphson root-finding
        // algorithm
        hyperbolicEccentricAnomaly_ = pointerToNewtonRaphson_ ->getComputedRootOfFunction( );
    }

    // Check if orbit is near-parabolic.
    else if ( eccentricity_ <= 1.2 )
    {
        cerr << "Orbit is near-parabolic and, at present conversion, between hyperbolic eccentric "
             << "anomaly and hyperbolic mean anomaly is not possible for eccentricities in the "
             << "range: 0.8 < eccentricity < 1.2." << endl;
    }

    // Return hyperbolic eccentric anomaly.
    return hyperbolicEccentricAnomaly_;
}

} // namespace orbital_element_conversions
} // tudat
