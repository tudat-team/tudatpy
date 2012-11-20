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
 *      120813    P. Musegaas       Changed code to new root finding structure.
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *
 *    Notes
 *      Currently, this conversion is only valid for eccentricities up to 0.97 due to the
 *      difficulties of the iterative method to converge for eccentricities close to 1.
 *
 */

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <boost/exception/all.hpp>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"
#include "Tudat/Mathematics/BasicMathematics/functionProxy.h"

namespace tudat
{
namespace basic_astrodynamics
{
namespace orbital_element_conversions
{

using namespace root_finders;
using namespace root_finders::termination_conditions;

//! Construct converter with eccentricity and mean anomaly.
ConvertMeanAnomalyToEccentricAnomaly::ConvertMeanAnomalyToEccentricAnomaly( 
        const double anEccentricity, 
        const double aMeanAnomaly,
        RootFinderPointer aRootFinder )
    : eccentricity( anEccentricity ),
      meanAnomaly( aMeanAnomaly ),
      rootFinder( aRootFinder )
{
    // Required because the make_shared in the function definition gives problems for MSVC.
    if ( !rootFinder.get( ) )
    {
        rootFinder = boost::make_shared< NewtonRaphson >(
            boost::bind( &RootAbsoluteToleranceTerminationCondition::checkTerminationCondition,
            boost::make_shared< RootAbsoluteToleranceTerminationCondition >( 5.0e-15, 1000 ),
            _1, _2, _3, _4, _5 ) );
    }
}

//! Convert mean anomaly to eccentric anomaly.
double ConvertMeanAnomalyToEccentricAnomaly::convert( )
{
    // Declare eccentric anomaly.
    double eccentricAnomaly = TUDAT_NAN;

    // Check if orbit is circular.
    if ( std::fabs( eccentricity ) < std::numeric_limits< double >::min( ) )
    {
        // If orbit is circular mean anomaly and eccentric anomaly are equal.
        eccentricAnomaly = meanAnomaly;
    }

    // Check if eccentricity is non-negative.
    else if ( eccentricity < 0.0 )
    {
        boost::throw_exception( boost::enable_error_info(
                                    std::runtime_error( "Orbit eccentricity is negative!" ) ) );
    }

    // Check if orbit is near-parabolic, i.e. eccentricity >= 0.98.
    else if ( eccentricity > 0.98 )
    {
        std::stringstream errorMessage;
        errorMessage << "Orbit is near-parabolic and, at present conversion, between eccentric "
                     << "anomaly and mean anomaly is not possible for eccentricities larger than: "
                     << "0.98" << std::endl;

        boost::throw_exception( boost::enable_error_info(
                                    std::runtime_error( errorMessage.str( ) ) ) );
    }

    // Check if orbit is elliptical, and not near-parabolic.
    else
    {
        // Create an object containing the function of which we whish to obtain the root from.
        basic_mathematics::UnivariateProxyPtr rootFunction
                = boost::make_shared< basic_mathematics::UnivariateProxy >(
                    boost::bind( &ConvertMeanAnomalyToEccentricAnomaly::
                                 computeKeplersFunctionForEllipticalOrbits, this, _1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding( -1, boost::bind( &ConvertMeanAnomalyToEccentricAnomaly::
                computeFirstDerivativeKeplersFunctionForEllipticalOrbits, this, _1 ) );

        // Set eccentric anomaly based on result of Newton-Raphson root-finding algorithm.
        eccentricAnomaly = rootFinder->execute( rootFunction, meanAnomaly );
    }

    // Return eccentric anomaly.
    return eccentricAnomaly;
}

} // namespace orbital_element_conversions
} // namespace basic_astrodynamics
} // namespace tudat
