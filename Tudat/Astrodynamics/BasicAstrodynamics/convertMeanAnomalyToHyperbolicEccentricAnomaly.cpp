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
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *      http://www.cdeagle.com/omnum/pdf/demokep1.pdf, last accessed: 16th February, 2011.
 *
 */

#include <iostream>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"

namespace tudat
{
namespace orbital_element_conversions
{

//! Convert mean anomaly to hyperbolic eccentric anomaly.
double ConvertMeanAnomalyToHyperbolicEccentricAnomaly::convert( )
{
    // Declare hyperbolic eccentric anomaly.
    double hyperbolicEccentricAnomaly_ = TUDAT_NAN;

    // Set the class that contains the functions needed for Newton-Raphson.
    newtonRaphsonAdaptor_.setClass( this );

    // Set NewtonRaphson adaptor class.
    newtonRaphson_->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptor_ );

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
        newtonRaphson_->setInitialGuessOfRoot( 2.0 * hyperbolicMeanAnomaly_
                                               / eccentricity_ - 1.8 );

        // Execute Newton-Raphon method.
        newtonRaphson_->execute( );

        // Set hyperbolic eccentric anomaly based on result of Newton-Raphson root-finding
        // algorithm
        hyperbolicEccentricAnomaly_ = newtonRaphson_ ->getComputedRootOfFunction( );
    }

    // Check if orbit is near-parabolic.
    else if ( eccentricity_ <= 1.2 )
    {
        std::cerr << "Orbit is near-parabolic and, at present conversion, between hyperbolic "
                  << "eccentric anomaly and hyperbolic mean anomaly is not possible for "
                  << "eccentricities in the range: 0.8 < eccentricity < 1.2." << std::endl;
    }

    // Return hyperbolic eccentric anomaly.
    return hyperbolicEccentricAnomaly_;
}

} // namespace orbital_element_conversions
} // namespace tudat
