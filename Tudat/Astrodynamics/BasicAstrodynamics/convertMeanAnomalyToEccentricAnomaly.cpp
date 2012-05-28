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
 *      110210    K. Kumar          Creation of code.
 *      111209    T. Secretin       Relaxed constraints on near-parabolic check.
 *      111221    T. Secretin       Added zero eccentricity case and check for negative
 *                                  eccentricities.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *
 *    Currently, this conversion is only valid for eccentricities up to 0.97 due to the
 *    difficulties of the iterative method to converge for eccentricities close to 1.
 *
 */

#include <iostream>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"

namespace tudat
{
namespace orbital_element_conversions
{

//! Convert mean anomaly to eccentric anomaly.
double ConvertMeanAnomalyToEccentricAnomaly::convert( )
{
    // Declare eccentric anomaly.
    double eccentricAnomaly_ = TUDAT_NAN;

    // Set the class that contains the functions needed for Newton-Raphson.
    newtonRaphsonAdaptor_.setClass( this );

    // Set NewtonRaphson adaptor class.
    newtonRaphson_->setNewtonRaphsonAdaptor( &newtonRaphsonAdaptor_ );

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
        newtonRaphson_->setInitialGuessOfRoot( meanAnomaly_ );

        // Execute Newton-Raphon method.
        newtonRaphson_->execute( );

        // Set eccentric anomaly based on result of Newton-Raphson root-finding algorithm.
        eccentricAnomaly_ = newtonRaphson_->getComputedRootOfFunction( );
    }

    // Check if orbit is near-parabolic, i.e. eccentricity >= 0.98.
    else
    {
        std::cerr << "Orbit is near-parabolic and, at present conversion, between eccentric "
                  << "anomaly and mean anomaly is not possible for eccentricities larger than: "
                  << "0.98" << std::endl;
    }

    // Return eccentric anomaly.
    return eccentricAnomaly_;
}

} // namespace orbital_element_conversions
} // namespace tudat
