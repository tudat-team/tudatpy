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
 *      120813    P. Musegaas       Changed code to new root finding structure.
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *      http://www.cdeagle.com/omnum/pdf/demokep1.pdf, last accessed: 16th February, 2011.
 *
 */

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <boost/exception/all.hpp>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"
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
ConvertMeanAnomalyToHyperbolicEccentricAnomaly::ConvertMeanAnomalyToHyperbolicEccentricAnomaly( 
        const double anEccentricity, 
        const double aHyperbolicMeanAnomaly,
        RootFinderPointer aRootFinder )
    : eccentricity( anEccentricity ),
      hyperbolicMeanAnomaly( aHyperbolicMeanAnomaly ),
      rootFinder( aRootFinder )
{
    // Required because the make_shared in the function definition gives problems for MSVC.
    if ( !rootFinder.get( ) )
    {
        rootFinder = boost::make_shared< NewtonRaphson >(
                    boost::bind(
                        &RootAbsoluteOrRelativeToleranceTerminationCondition::
                        checkTerminationCondition,
                        boost::make_shared< RootAbsoluteOrRelativeToleranceTerminationCondition >(
                            5.0e-15, 1000 ), _1, _2, _3, _4, _5 ) );
    }
}

//! Convert mean anomaly to hyperbolic eccentric anomaly.
double ConvertMeanAnomalyToHyperbolicEccentricAnomaly::convert( )
{
    // Declare hyperbolic eccentric anomaly.
    double hyperbolicEccentricAnomaly = TUDAT_NAN;

    // Check if orbit is near-parabolic, i.e. eccentricity <= 1.2.
    if ( eccentricity <= 1.2 )
    {
        std::stringstream errorMessage;
        errorMessage << "Orbit is near-parabolic and, at present conversion, between hyperbolic "
                     << "eccentric anomaly and hyperbolic mean anomaly is not possible for "
                     << "eccentricities in the range: 0.8 < eccentricity < 1.2." << std::endl;

        boost::throw_exception( boost::enable_error_info(
                                    std::runtime_error( errorMessage.str( ) ) ) );
    }

    // Check orbit is hyperbolic, and not near-parabolic.
    else
    {
        // Create an object containing the function of which we whish to obtain the root from.
        basic_mathematics::UnivariateProxyPtr rootFunction
                = boost::make_shared< basic_mathematics::UnivariateProxy >(
                    boost::bind( &ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
                                 computeKeplersFunctionForHyperbolicOrbits, this, _1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding(
                    -1, boost::bind( &ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
                                     computeFirstDerivativeKeplersFunctionForHyperbolicOrbits,
                                     this, _1 ) );

        // Set hyperbolic eccentric anomaly based on result of Newton-Raphson root-finding
        // algorithm.
        hyperbolicEccentricAnomaly = rootFinder->execute( rootFunction,
                2.0 * hyperbolicMeanAnomaly / eccentricity - 1.8 );
    }

    // Return hyperbolic eccentric anomaly.
    return hyperbolicEccentricAnomaly;
}

} // namespace orbital_element_conversions
} // namespace basic_astrodynamics
} // namespace tudat
