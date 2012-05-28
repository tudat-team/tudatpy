/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110214    K. Kumar          Creation of code.
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
