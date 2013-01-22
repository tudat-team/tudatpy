/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      110210    K. Kumar          File created.
 *      111209    T. Secretin       Relaxed constraints on near-parabolic check.
 *      111221    T. Secretin       Added zero eccentricity case and check for negative
 *                                  eccentricities.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120813    P. Musegaas       Changed code to new root finding structure.
 *      120822    P. Musegaas       Tested and improved initial guess of eccentric anomaly. Added
 *                                  functionality for near-parabolic cases. Added option for user
 *                                  to specifiy initial guess. Changed error message to a runtime
 *                                  error.
 *      120903    P. Musegaas       Added additional warning message.
 *      121205    P. Musegaas       Updated code to final version of rootfinders.
 *      130116    E. Heeren         Changed RootAbsoluteToleranceTerminationCondition to work on
 *                                  all systems.
 *
 *    References
 *      Regarding method in general:
 *          Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *          Wakker, K. F. Astrodynamics I + II. Lecture Notes AE4-874, Delft University of
 *              Technology, Delft, Netherlands.
 *      Regarding the choice of initial guess:
 *          Musegaas, P., Optimization of Space Trajectories Including Multiple Gravity Assists and
 *              Deep Space Maneuvers, MSc thesis report, Delft University of Technology, 2012.
 *              [unpublished so far]. Section available on tudat website (tudat.tudelft.nl)
 *              under issue #539.
 *
 *    Notes
 *
 */

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <boost/exception/all.hpp>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>

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
        const double anEccentricity, const double aMeanAnomaly,
        const bool useDefaultInitialGuess_,
        const double userSpecifiedInitialGuess_,
        root_finders::RootFinderPointer aRootFinder )
    : eccentricity( anEccentricity ),
      useDefaultInitialGuess( useDefaultInitialGuess_ ),
      userSpecifiedInitialGuess( userSpecifiedInitialGuess_ ),
      rootFinder( aRootFinder )

{
    // Set mean anomaly to region between 0 and 2 PI.
    meanAnomaly = mathematics::computeModulo( aMeanAnomaly, 2.0 * mathematics::PI );

    // Required because the make_shared in the function definition gives problems for MSVC.
    if ( !rootFinder.get( ) )
    {
        rootFinder = boost::make_shared< NewtonRaphson >(
            boost::bind( &RootAbsoluteToleranceTerminationCondition::checkTerminationCondition,
            boost::make_shared< RootAbsoluteToleranceTerminationCondition >( 1.0e-13, 1000 ),
            _1, _2, _3, _4, _5 ) );
    }
}

//! Convert mean anomaly to eccentric anomaly.
double ConvertMeanAnomalyToEccentricAnomaly::convert( )
{
    // Declare eccentric anomaly.
    double eccentricAnomaly = TUDAT_NAN;

    // Check if orbit is elliptical, and not near-parabolic.
    if ( eccentricity < 1.0 && eccentricity >= 0.0 )
    {
        // Create an object containing the function of which we whish to obtain the root from.
        basic_mathematics::UnivariateProxyPointer rootFunction
                = boost::make_shared< basic_mathematics::UnivariateProxy >(
                    boost::bind( &ConvertMeanAnomalyToEccentricAnomaly::
                                 computeKeplersFunctionForEllipticalOrbits, this, _1 ) );

        // Add the first derivative of the root function.
        rootFunction->addBinding( -1, boost::bind( &ConvertMeanAnomalyToEccentricAnomaly::
                computeFirstDerivativeKeplersFunctionForEllipticalOrbits, this, _1 ) );

        // Declare initial guess.
        double initialGuess = TUDAT_NAN;

        // Set the initial guess. Check if the default scheme is to be used or a user specified
        // value should be used.
        // !!!!!!!!!!!!!     IMPORTANT     !!!!!!!!!!!!!
        // If this scheme is changed, please run a very extensive test suite. The root finder
        // function tends to be chaotic for some very specific combinations of mean anomaly and
        // eccentricity. Various random tests of 100.000.000 samples were done to verify the
        // functionality of this one, and of another option for the starter: PI. [Musegaas,2012]
        if ( useDefaultInitialGuess )
        {
            if ( meanAnomaly > mathematics::PI )
            {
                initialGuess =  meanAnomaly - eccentricity;
            }
            else
            {
                initialGuess = meanAnomaly + eccentricity;
            }
        }
        else
        {
            initialGuess = userSpecifiedInitialGuess;
        }

        // Set eccentric anomaly based on result of Newton-Raphson root-finding algorithm.
        eccentricAnomaly = rootFinder->execute( rootFunction, initialGuess );
    }

    //  Eccentricity is invalid: eccentricity < 0.0 or eccentricity >= 1.0.
    else
    {
        boost::throw_exception( std::runtime_error( boost::str( boost::format(
                "Invalid eccentricity. Valid range is 0.0 <= e < 1.0. Eccentricity was: '%f'." )
                    % eccentricity ) ) );
    }

    // Return eccentric anomaly.
    return eccentricAnomaly;
}

} // namespace orbital_element_conversions
} // namespace basic_astrodynamics
} // namespace tudat
