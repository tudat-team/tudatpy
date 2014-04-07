/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      110214    K. Kumar          File created.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120813    P. Musegaas       Changed code to new root finding structure.
 *      120822    P. Musegaas       Tested and improved initial guess of hyperbolic eccentric
 *                                  anomaly. Added functionality for near-parabolic cases. Added
 *                                  option for user to specifiy initial guess. Changed error
 *                                  message to a runtime error.
 *      120903    P. Musegaas       Added additional warning message.
 *      121205    P. Musegaas       Updated code to final version of rootfinders.
 *
 *    References
 *      Regarding method in general, including starter values used:
 *          Wakker, K. F. Astrodynamics I + II. Lecture Notes AE4-874, Delft University of
 *              Technology, Delft, Netherlands.
 *      Regarding the performance of different starter values and performance near-parabolic.
 *          Musegaas, P., Optimization of Space Trajectories Including Multiple Gravity Assists and
 *              Deep Space Maneuvers, MSc thesis report, Delft University of Technology, 2012.
 *              [unpublished so far]. Section available on tudat website (tudat.tudelft.nl)
 *              under issue #539.
 *      Regarding older versions of the code (120813 and before):
 *          Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *          http://www.cdeagle.com/omnum/pdf/demokep1.pdf, last accessed: 16th February, 2011.
 *
 *    Notes
 *
 */

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <boost/math/special_functions/asinh.hpp>

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
        const bool useDefaultInitialGuess_,
        const double userSpecifiedInitialGuess_,
        RootFinderPointer aRootFinder )
    : eccentricity( anEccentricity ),
      hyperbolicMeanAnomaly( aHyperbolicMeanAnomaly ),
      useDefaultInitialGuess( useDefaultInitialGuess_ ),
      userSpecifiedInitialGuess( userSpecifiedInitialGuess_ ),
      rootFinder( aRootFinder )
{
    // Required because the make_shared in the function definition gives problems for MSVC.
    if ( !rootFinder.get( ) )
    {
        rootFinder = boost::make_shared< NewtonRaphson >(
                    boost::bind(
                        &RootAbsoluteToleranceTerminationCondition::
                        checkTerminationCondition,
                        boost::make_shared< RootAbsoluteToleranceTerminationCondition >(
                            5.0e-15, 1000 ), _1, _2, _3, _4, _5 ) );
    }
}

//! Convert mean anomaly to hyperbolic eccentric anomaly.
double ConvertMeanAnomalyToHyperbolicEccentricAnomaly::convert( )
{
    // Declare hyperbolic eccentric anomaly.
    double hyperbolicEccentricAnomaly = TUDAT_NAN;

    // Check if orbit is hyperbolic.
    if ( eccentricity > 1.0 )
    {
        // Create an object containing the function of which we whish to obtain the root from.
        basic_mathematics::UnivariateProxyPointer rootFunction
                = boost::make_shared< basic_mathematics::UnivariateProxy >(
                    boost::bind( &ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
                                 computeKeplersFunctionForHyperbolicOrbits, this, _1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding(
                    -1, boost::bind( &ConvertMeanAnomalyToHyperbolicEccentricAnomaly::
                                     computeFirstDerivativeKeplersFunctionForHyperbolicOrbits,
                                     this, _1 ) );

        // Declare initial guess.
        double initialGuess = TUDAT_NAN;

        // Set the initial guess. Check if the default scheme is to be used or a user specified
        // value should be used. See [Wakker, 2007] for derivations of the default values.
        // Note that an error was detected in these starter values, as is discussed in
        // [Musegaas,2012].
        // !!!!!!!!!!!!!     IMPORTANT     !!!!!!!!!!!!!
        // If this scheme is changed, please run a very extensive test suite. The root finder
        // function tends to be chaotic for some very specific combinations of mean anomaly and
        // eccentricity. Various random tests of 100.000.000 samples were done to verify the
        // functionality of this one. [Musegaas,2012]
        if ( useDefaultInitialGuess )
        {
            if ( std::abs( hyperbolicMeanAnomaly ) < 6.0 * eccentricity )
            {
                initialGuess =
                        std::sqrt( 8.0 * ( eccentricity - 1.0 ) / eccentricity ) *
                        std::sinh( 1.0 / 3.0 * boost::math::asinh( 3.0 * hyperbolicMeanAnomaly /
                                ( std::sqrt( 8.0 * ( eccentricity - 1.0 ) / eccentricity ) *
                                             ( eccentricity - 1.0 ) ) ) );
            }
            else if ( hyperbolicMeanAnomaly > 6.0 * eccentricity )
            {
                initialGuess = ( std::log( 2.0 * hyperbolicMeanAnomaly / eccentricity ) );
            }
            else
            {
                initialGuess = ( - std::log( -2.0 * hyperbolicMeanAnomaly / eccentricity ) );
            }
        }
        else
        {
            initialGuess = userSpecifiedInitialGuess;
        }

        // Set hyperbolic eccentric anomaly based on result of Newton-Raphson root-finding
        // algorithm.
        hyperbolicEccentricAnomaly = rootFinder->execute( rootFunction, initialGuess );
    }

    // In this case the orbit is not hyperbolic.
    else
    {
        boost::throw_exception( std::runtime_error( boost::str( boost::format(
                "Invalid eccentricity. Valid range is e > 1.0. Eccentricity was: '%f'." )
                    % eccentricity ) ) );
    }

    // Return hyperbolic eccentric anomaly.
    return hyperbolicEccentricAnomaly;
}

} // namespace orbital_element_conversions
} // namespace basic_astrodynamics
} // namespace tudat
