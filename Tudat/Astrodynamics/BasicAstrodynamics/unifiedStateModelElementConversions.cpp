/*    Copyright (c) 2010-2016, Delft University of Technology
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
 *      160406    M. Van den Broeck Creation
 *      <YYMMDD>  <author name>     <comment>
 *
 *    References
 *      <First reference>
 *      <Second reference>
 *
 *    Notes
 *
 */

#include <cmath>

#include <boost/exception/all.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian to Unified State Model elements.
basic_mathematics::Vector6d convertKeplerianToModifiedEquinoctialElements(
        const basic_mathematics::Vector6d& keplerianElements,
        const double centralBodyGravitationalParameter )
// Based on Vittaldev, 2010.
{

    // Declaring eventual output vector.
    basic_mathematics::Vector6d UnifiedStateModelState( 7 );

    // Compute the C hodograph element of the Unified State Model
    UnifiedStateModelState( CHodographIndex ) =
            std::sqrt( centralBodyGravitationalParameter / ( keplerianElements( semiMajorAxisIndex )
                                                  * ( 1 - keplerianElements( eccentricityIndex ) *
                                                      keplerianElements( eccentricityIndex ) ) ) );

    // Calculate the additional R hodograph parameter
    double RHodographElement = keplerianElements( eccentricityIndex ) * UnifiedStateModelState( CHodographIndex );

    // Compute the Rf1 hodograph element of the Unified State Model
    UnifiedStateModelState( Rf1HodographIndex ) =
            - RHodographElement * std::sin( keplerianElements( longitudeOfAscendingNodeIndex )
                                            + keplerianElements( argumentOfPeriapsisIndex ) );

    // Compute the Rf2 hodograph element of the Unified State Model
    UnifiedStateModelState( Rf2HodographIndex ) =
              RHodographElement * std::cos( keplerianElements( longitudeOfAscendingNodeIndex )
                                            + keplerianElements( argumentOfPeriapsisIndex ) );

    // Calculate the additional argument of longitude u
    double argumentOfLongitude = keplerianElements( argumentOfPeriapsisIndex ) +
            keplerianElements( trueAnomalyIndex );

    // Compute the epsilon1 quaternion of the Unified State Model
    UnifiedStateModelState( epsilon1QuaternionIndex ) =
            std::sin( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::cos( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) - argumentOfLongitude ) );

    // Compute the epsilon2 quaternion of the Unified State Model
    UnifiedStateModelState( epsilon2QuaternionIndex ) =
            std::sin( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::sin( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) - argumentOfLongitude ) );

    // Compute the epsilon3 quaternion of the Unified State Model
    UnifiedStateModelState( epsilon3QuaternionIndex ) =
            std::cos( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::sin( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) + argumentOfLongitude ) );

    // Compute the eta quaternion of the Unified State Model
    UnifiedStateModelState( etaQuaternionIndex ) =
            std::cos( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::cos( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) + argumentOfLongitude ) );

    // Give back result
    return UnifiedStateModelState;

}
