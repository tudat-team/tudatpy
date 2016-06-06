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
 *      Vittaldev, V. (2010). The Unified State Model: Derivation and application in astrodynamics
 *          and navigation. Master's thesis, Delft University of Technology.
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
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian elements to Unified State Model elements.
Eigen::Matrix< double, 7, 1 > convertKeplerianToUnifiedStateModelElements(
        const basic_mathematics::Vector6d& keplerianElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Matrix< double, 7, 1 > convertedUnifiedStateModelElements = Eigen::Matrix< double, 7, 1 >::Zero( );

    // Define the tolerance of a singularity
    double singularityTolerance = 1.0e-15; // Based on tolerance chosen in
                                           // orbitalElementConversions.cpp in Tudat Core.

    // If eccentricity is outside range [0,inf)
    if ( keplerianElements( eccentricityIndex ) < 0.0 )
    {
        //Define the error message
        std::stringstream errorMessage;
        errorMessage << "Eccentricity is expected in range [0,inf)\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << std::endl;

        // Throw exception
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }

    // If inclination is outside range [0,PI]
    if ( ( keplerianElements( inclinationIndex ) < 0.0 ) || ( keplerianElements( inclinationIndex ) > PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "Inclination is expected in range [0," << PI << "]\n"
                     << "Specified inclination: " << keplerianElements( inclinationIndex ) << " rad." << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }

    // If argument of pericenter is outside range [0,2.0 * PI]
    if ( ( keplerianElements( argumentOfPeriapsisIndex ) < 0.0 ) || ( keplerianElements( argumentOfPeriapsisIndex ) >
                                                                      2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "RAAN is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( argumentOfPeriapsisIndex ) << " rad." << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }

    // If right ascension of ascending node is outside range [0,2.0 * PI]
    if ( ( keplerianElements( longitudeOfAscendingNodeIndex ) < 0.0 ) ||
         ( keplerianElements( longitudeOfAscendingNodeIndex ) > 2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "RAAN is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( longitudeOfAscendingNodeIndex ) << " rad."
                     << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }

    // If true anomaly is outside range [0,2.0 * PI]
    if ( ( keplerianElements( trueAnomalyIndex ) < 0.0 ) || ( keplerianElements( trueAnomalyIndex ) > 2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "RAAN is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( trueAnomalyIndex ) << " rad." << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }

    // If inclination is zero and the right ascension of ascending node is non-zero
    if ( ( std::fabs( keplerianElements( inclinationIndex ) ) < singularityTolerance ) &&
         ( std::fabs( keplerianElements( longitudeOfAscendingNodeIndex ) ) > singularityTolerance ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the inclination is zero, the right ascending node should be zero by definition\n"
                     << "Specified right ascension of ascending node: " <<
                        keplerianElements( longitudeOfAscendingNodeIndex ) << " rad." << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }

    // If eccentricity is zero and the argument of pericenter is non-zero
    if ( ( std::fabs( keplerianElements( eccentricityIndex ) ) < singularityTolerance ) &&
         ( std::fabs( keplerianElements( argumentOfPeriapsisIndex ) ) > singularityTolerance ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the eccentricity is zero, the argument of pericenter should be zero by definition\n"
                     << "Specified argument of pericenter: " <<
                        keplerianElements( argumentOfPeriapsisIndex ) << " rad." << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }

    // If semi-major axis is negative and the eccentricity is smaller or equal to one
    if ( ( keplerianElements( semiMajorAxisIndex ) < 0.0 ) && ( keplerianElements( eccentricityIndex ) <= 1.0 ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the semi-major axis is negative, the eccentricity should be larger than one\n"
                     << "Specified semi-major axis: " << keplerianElements( semiMajorAxisIndex ) << " m.\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << " rad." << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }

    // If semi-major axis is positive and the eccentricity is larger than one
    if ( ( keplerianElements( semiMajorAxisIndex ) > 0.0 ) && ( keplerianElements( eccentricityIndex ) > 1.0 ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the semi-major axis is positive, the eccentricity should be smaller than or equal to one\n"
                     << "Specified semi-major axis: " << keplerianElements( semiMajorAxisIndex ) << " m.\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << " rad." << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }
    //Else, nothing wrong and continue

    // Compute the C hodograph element of the Unified State Model
    if ( std::fabs( keplerianElements( eccentricityIndex ) - 1.0) < singularityTolerance )
            // parabolic orbit -> semi-major axis is not defined
    {
        convertedUnifiedStateModelElements( CHodographIndex ) =
                std::sqrt( centralBodyGravitationalParameter / keplerianElements( semiLatusRectumIndex ) );
    }
    else
    {
        convertedUnifiedStateModelElements( CHodographIndex ) =
                std::sqrt( centralBodyGravitationalParameter / ( keplerianElements( semiMajorAxisIndex )
                                                  * ( 1 - keplerianElements( eccentricityIndex ) *
                                                      keplerianElements( eccentricityIndex ) ) ) );
    }

    // Calculate the additional R hodograph parameter
    double RHodographElement = keplerianElements( eccentricityIndex ) *
            convertedUnifiedStateModelElements( CHodographIndex );

    // Compute the Rf1 hodograph element of the Unified State Model
    convertedUnifiedStateModelElements( Rf1HodographIndex ) =
            - RHodographElement * std::sin( keplerianElements( longitudeOfAscendingNodeIndex )
                                            + keplerianElements( argumentOfPeriapsisIndex ) );

    // Compute the Rf2 hodograph element of the Unified State Model
    convertedUnifiedStateModelElements( Rf2HodographIndex ) =
              RHodographElement * std::cos( keplerianElements( longitudeOfAscendingNodeIndex )
                                            + keplerianElements( argumentOfPeriapsisIndex ) );

    // Calculate the additional argument of longitude u
    double argumentOfLongitude = keplerianElements( argumentOfPeriapsisIndex ) +
            keplerianElements( trueAnomalyIndex );

    // Compute the epsilon1 quaternion of the Unified State Model
    convertedUnifiedStateModelElements( epsilon1QuaternionIndex ) =
            std::sin( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::cos( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) - argumentOfLongitude ) );

    // Compute the epsilon2 quaternion of the Unified State Model
    convertedUnifiedStateModelElements( epsilon2QuaternionIndex ) =
            std::sin( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::sin( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) - argumentOfLongitude ) );

    // Compute the epsilon3 quaternion of the Unified State Model
    convertedUnifiedStateModelElements( epsilon3QuaternionIndex ) =
            std::cos( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::sin( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) + argumentOfLongitude ) );

    // Compute the eta quaternion of the Unified State Model
    convertedUnifiedStateModelElements( etaQuaternionIndex ) =
            std::cos( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::cos( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) + argumentOfLongitude ) );

    // Give back result
    return convertedUnifiedStateModelElements;

}

//! Convert Unified State Model elements to Keplerian elements.
basic_mathematics::Vector6d convertUnifiedStateModelToKeplerianElements(
        const Eigen::Matrix< double, 7, 1 >& unifiedStateModelElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    basic_mathematics::Vector6d convertedKeplerianElements = basic_mathematics::
            Vector6d::Zero( 6 );

    // Define the tolerance of a singularity
    double singularityTolerance = 1.0e-15; // Based on tolerance chosen in
                                           // orbitalElementConversions.cpp in Tudat Core.

    // Declare auxiliary parameters before using them in the if loop
    double cosineLambda = 0.0;
    double sineLambda = 0.0;
    double lambdaFromSineAndCosine = 0.0;

    // Check whether the Unified State Model elements are within expected limits
    // If inclination is zero and the right ascension of ascending node is non-zero
    const double normOfQuaternionElements = std::sqrt( std::pow( unifiedStateModelElements( epsilon1QuaternionIndex ), 2 ) +
                                                       std::pow( unifiedStateModelElements( epsilon2QuaternionIndex ), 2 ) +
                                                       std::pow( unifiedStateModelElements( epsilon3QuaternionIndex ), 2 ) +
                                                       std::pow( unifiedStateModelElements( etaQuaternionIndex ), 2 ) );

    if ( std::fabs( normOfQuaternionElements - 1.0 ) > singularityTolerance )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "The norm of the quaternion should be equal to one.\n"
                     << "Norm of the specified quaternion is: " << normOfQuaternionElements << " ." << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }
    //Else, nothing wrong and continue

    // Compute auxiliary parameters cosineLambda, sineLambda and Lambda
    if ( ( std::fabs( unifiedStateModelElements( epsilon3QuaternionIndex ) ) < singularityTolerance )
        && ( std::fabs( unifiedStateModelElements( etaQuaternionIndex ) ) < singularityTolerance ) )
            // pure-retrograde orbit -> inclination  = pi
    {
        //Define the error message
        std::stringstream errorMessage;
        errorMessage << "Pure-retrograde orbit (inclination = pi).\n"
                     << "Unified State Model elements cannot be transformed to Kepler elements." << std::endl;

        // Throw exception
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }
    else
    {
        cosineLambda = ( unifiedStateModelElements( etaQuaternionIndex ) *
                                unifiedStateModelElements( etaQuaternionIndex ) -
                                unifiedStateModelElements( epsilon3QuaternionIndex ) *
                                unifiedStateModelElements( epsilon3QuaternionIndex ) )
                / ( unifiedStateModelElements( epsilon3QuaternionIndex ) *
                    unifiedStateModelElements( epsilon3QuaternionIndex ) +
                    unifiedStateModelElements( etaQuaternionIndex ) *
                    unifiedStateModelElements( etaQuaternionIndex ) );
        sineLambda = ( 2.0 *
                              unifiedStateModelElements( epsilon3QuaternionIndex ) *
                              unifiedStateModelElements( etaQuaternionIndex ) )
                / ( unifiedStateModelElements( epsilon3QuaternionIndex ) *
                    unifiedStateModelElements( epsilon3QuaternionIndex ) +
                    unifiedStateModelElements( etaQuaternionIndex ) *
                    unifiedStateModelElements( etaQuaternionIndex ) );
        lambdaFromSineAndCosine = std::atan2( sineLambda, cosineLambda );
    }

    // Compute auxiliary parameters auxiliaryParameter1 and auxiliaryParameter2
    double auxiliaryParameter1 = unifiedStateModelElements( Rf1HodographIndex ) * cosineLambda +
            unifiedStateModelElements( Rf2HodographIndex ) * sineLambda;
    double auxiliaryParameter2 = unifiedStateModelElements( CHodographIndex ) -
            unifiedStateModelElements( Rf1HodographIndex ) * sineLambda +
            unifiedStateModelElements( Rf2HodographIndex ) * cosineLambda;

    // Compute auxiliary R hodograph parameter
    double RHodographElement = std::sqrt( unifiedStateModelElements( Rf1HodographIndex )
                                          * unifiedStateModelElements( Rf1HodographIndex )
                                          + unifiedStateModelElements( Rf2HodographIndex )
                                          * unifiedStateModelElements( Rf2HodographIndex ));

    // Compute eccentricity
    convertedKeplerianElements( eccentricityIndex ) =
            RHodographElement / unifiedStateModelElements( CHodographIndex );

    // Compute semi-major axis or, in case of a parabolic orbit, the semi-latus rectum.
    if ( std::fabs( convertedKeplerianElements( eccentricityIndex ) - 1.0 ) < singularityTolerance )
            // parabolic orbit -> semi-major axis is not defined. Use semi-latus rectum instead.
    {
        convertedKeplerianElements( semiLatusRectumIndex ) = centralBodyGravitationalParameter /
                ( unifiedStateModelElements( CHodographIndex ) * unifiedStateModelElements( CHodographIndex ) );
    }
    else
    {
        convertedKeplerianElements( semiMajorAxisIndex ) =
                centralBodyGravitationalParameter /
                ( 2.0 * unifiedStateModelElements( CHodographIndex ) * auxiliaryParameter2 -
                    ( auxiliaryParameter1 * auxiliaryParameter1 + auxiliaryParameter2 * auxiliaryParameter2 ) );
    }

    // Compute inclination
    convertedKeplerianElements( inclinationIndex ) =
            std::acos( 1.0 - 2.0 * ( unifiedStateModelElements( epsilon1QuaternionIndex ) *
                                     unifiedStateModelElements( epsilon1QuaternionIndex ) +
                                     unifiedStateModelElements( epsilon2QuaternionIndex ) *
                                     unifiedStateModelElements( epsilon2QuaternionIndex ) ) );
        // This acos is always defined correctly because the inclination is always below pi rad.

    // Compute longitude of ascending node
    if ( ( ( std::fabs( unifiedStateModelElements( epsilon1QuaternionIndex ) ) < singularityTolerance )
           && ( std::fabs( unifiedStateModelElements( epsilon2QuaternionIndex ) ) < singularityTolerance ) ) ||
         ( ( std::fabs( unifiedStateModelElements( epsilon3QuaternionIndex ) ) < singularityTolerance )
         && ( std::fabs( unifiedStateModelElements( etaQuaternionIndex ) ) < singularityTolerance ) ) )
            // pure-prograde or pure-retrograde orbit
    {
        convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0; // by definition
    }
    else
    {
        convertedKeplerianElements( longitudeOfAscendingNodeIndex ) =
                std::atan2( ( ( unifiedStateModelElements( epsilon1QuaternionIndex ) *
                                unifiedStateModelElements( epsilon3QuaternionIndex ) +
                                unifiedStateModelElements( epsilon2QuaternionIndex ) *
                                unifiedStateModelElements( etaQuaternionIndex ) )
                           / ( std::sqrt( ( unifiedStateModelElements( epsilon1QuaternionIndex ) *
                                               unifiedStateModelElements( epsilon1QuaternionIndex ) +
                                               unifiedStateModelElements( epsilon2QuaternionIndex ) *
                                               unifiedStateModelElements( epsilon2QuaternionIndex ) ) *
                                          ( unifiedStateModelElements( etaQuaternionIndex ) *
                                               unifiedStateModelElements( etaQuaternionIndex ) +
                                               unifiedStateModelElements( epsilon3QuaternionIndex ) *
                                               unifiedStateModelElements( epsilon3QuaternionIndex ) ) ) ) ),
                            ( ( unifiedStateModelElements( epsilon1QuaternionIndex ) *
                             unifiedStateModelElements( etaQuaternionIndex ) -
                             unifiedStateModelElements( epsilon2QuaternionIndex ) *
                             unifiedStateModelElements( epsilon3QuaternionIndex ) )
                        / ( std::sqrt( ( unifiedStateModelElements( epsilon1QuaternionIndex ) *
                                            unifiedStateModelElements( epsilon1QuaternionIndex ) +
                                            unifiedStateModelElements( epsilon2QuaternionIndex ) *
                                            unifiedStateModelElements( epsilon2QuaternionIndex ) ) *
                                       ( unifiedStateModelElements( etaQuaternionIndex ) *
                                            unifiedStateModelElements( etaQuaternionIndex ) +
                                            unifiedStateModelElements( epsilon3QuaternionIndex ) *
                                            unifiedStateModelElements( epsilon3QuaternionIndex ) ) ) ) ) );

        // Round off small values of the right ascension of ascending node to zero
        if ( std::fabs( convertedKeplerianElements( longitudeOfAscendingNodeIndex ) ) < singularityTolerance )
        {
            convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
        }
        // Ensure the longitude of ascending node is positive
        while ( convertedKeplerianElements( longitudeOfAscendingNodeIndex ) < 0.0 )
                // Because of the previous if loop, if the longitude of ascending node is smaller than 0, it will
                // always be smaller than -singularityTolerance
        {
            convertedKeplerianElements( longitudeOfAscendingNodeIndex ) =
                    convertedKeplerianElements( longitudeOfAscendingNodeIndex ) + 2.0 * PI;
        }
    }

    // Compute true anomaly and argument of periapsis
    if ( std::fabs( RHodographElement ) < singularityTolerance ) // circular orbit
    {
        convertedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0; // by definition
        convertedKeplerianElements( trueAnomalyIndex ) =
                lambdaFromSineAndCosine - convertedKeplerianElements( longitudeOfAscendingNodeIndex );

        // Round off small theta to zero
        if ( std::fabs( convertedKeplerianElements( trueAnomalyIndex ) ) < singularityTolerance )
        {
            convertedKeplerianElements( trueAnomalyIndex ) = 0.0;
        }

        // Ensure the true anomaly is positive
        while ( convertedKeplerianElements( trueAnomalyIndex ) < 0.0 )
                // Because of the previous if loop, if the true anomaly is smaller than zero, it will always be smaller than
                // -singularityTolerance
        {
            convertedKeplerianElements( trueAnomalyIndex ) =
                    convertedKeplerianElements( trueAnomalyIndex ) + 2.0 * PI;
        }
    }
    else
    {
        convertedKeplerianElements( trueAnomalyIndex ) =
                std::atan2( ( auxiliaryParameter1 / RHodographElement ),
                            ( ( auxiliaryParameter2 - unifiedStateModelElements( CHodographIndex ) )
                           / RHodographElement ) );

        // Round off small theta to zero
        if ( std::fabs( convertedKeplerianElements( trueAnomalyIndex ) ) < singularityTolerance )
        {
            convertedKeplerianElements( trueAnomalyIndex ) = 0.0;
        }

        // Ensure the true anomaly is positive
        while ( convertedKeplerianElements( trueAnomalyIndex ) < 0.0 )
            // Because of the previous if loop, if the true anomaly is smaller than zero, it will always
            // be smaller than -singularityTolerance
        {
            convertedKeplerianElements( trueAnomalyIndex ) =
                    convertedKeplerianElements( trueAnomalyIndex ) + 2.0 * PI;
        }

        convertedKeplerianElements( argumentOfPeriapsisIndex ) =
                lambdaFromSineAndCosine -
                convertedKeplerianElements( longitudeOfAscendingNodeIndex ) -
                convertedKeplerianElements( trueAnomalyIndex );

        // Round off small omega to zero
        if ( std::fabs( convertedKeplerianElements( argumentOfPeriapsisIndex ) ) < singularityTolerance )
        {
            convertedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;
        }

        // Ensure the argument of periapsis is positive
        while ( convertedKeplerianElements( argumentOfPeriapsisIndex ) < 0.0 )
            // Because of the previous if loop, if the argument of pericenter is smaller than zero,
            // it will be smaller than -singularityTolerance
        {
            convertedKeplerianElements( argumentOfPeriapsisIndex ) =
                    convertedKeplerianElements( argumentOfPeriapsisIndex ) + 2.0 * PI;
        }
    }

    // Give back result
    return convertedKeplerianElements;
}

} // close namespace orbital_element_conversions

} // close namespace tudat
