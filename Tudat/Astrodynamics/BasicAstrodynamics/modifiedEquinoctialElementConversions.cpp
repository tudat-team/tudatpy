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
 *      121128    R.C.A. Boon       File created (in progress).
 *      130124    R.C.A. Boon       Removed retrograde function, reworked Kepler to MEE conversion
 *                                  to properly convert inclination to range [0,180) and use
 *                                  retrograde factor internally
 *      130131    R.C.A. Boon       Added Cartesian conversion cases, added avoidSingularityAtPi~
 *                                  boolean flag to Kepler conversions, optimized computation
 *                                  (partially)
 *      130225    D. Dirkx          Added overloaded function for Kepler to MEE that determines
 *                                  retrogradeness based on Kepler state
 *      130301    R.C.A. Boon       Updated use of mathematics::PI to basic_mathematics::
 *                                  mathematical_constants::PI, minor textual changes.
 *      130305    R.C.A. Boon       Replaced Eigen::VectorXd by tudat::basic_mathematics::Vector6d
 *
 *    References
 *      Verified Interval Propagation, Bart Rˆmgens; Delft (2011). Code archive.
 *      Modified Equinoctial Orbital Elements, author unknown;
 *          http://www.cdeagle.com/pdf/mee.pdf (2010?).
 *      Survey of Orbital Element Sets, Gerald R. Hintz; Journal of Guidance, Control and
 *          Dynamics (2008, Vol. 31 - Nr. 3).
 *      Code archive, E. Heeren (fellow Tudat developer).
 *
 *    Notes
 *
 */

#include <cmath>

#include <boost/exception/all.hpp>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{
namespace basic_astrodynamics
{
namespace orbital_element_conversions
{

//! Convert Keplerian to modified equinoctial orbital elements using implicit MEE equation set.
tudat::basic_mathematics::Vector6d convertKeplerianToModifiedEquinoctialElements(
        const tudat::basic_mathematics::Vector6d& keplerianElements )
// Based on Hintz, 2008.
{
    using tudat::basic_mathematics::mathematical_constants::PI;

    // Check if orbit is retrograde
    bool avoidSingularityAtPiInclination =
            mission_geometry::isOrbitRetrograde( keplerianElements );

    // Convert to modified equinoctial elements.
    return convertKeplerianToModifiedEquinoctialElements( keplerianElements,
                                                          avoidSingularityAtPiInclination );
}

//! Convert Keplerian to modified equinoctial orbital elements using MEE explicit equation set.
tudat::basic_mathematics::Vector6d convertKeplerianToModifiedEquinoctialElements(
        const tudat::basic_mathematics::Vector6d& keplerianElements,
        const bool avoidSingularityAtPiInclination )
// Based on Hintz, 2008.
{
    using basic_mathematics::mathematical_constants::PI;

    // Declaring eventual output vector.
    tudat::basic_mathematics::Vector6d modifiedEquinoctialState( 6 );

    // Compute semi-latus rectum.
    double singularityTolerance = 1.0e-15; // Based on tolerance chosen in
                                           // orbitalElementConversions.cpp in Tudat Core.

    // Extracting eccentricity for ease of referencing.
    double eccentricity = keplerianElements( eccentricityIndex );

    // If e is (very near) one, then semi-major axis is undefined and thus the first kepler
    // element is the semi-latus rectum.
    if ( std::fabs( eccentricity - 1.0 ) < singularityTolerance )
    {
        modifiedEquinoctialState( semiLatusRectumIndex )
                = keplerianElements( semiLatusRectumIndex );
    }
    else // eccentricity is significantly away from singularity and can be computed with p=a(1-e≤).
    {
        modifiedEquinoctialState( semiLatusRectumIndex ) = keplerianElements( semiMajorAxisIndex )
                                                           * ( 1.0 - eccentricity * eccentricity );
    }

    // Determine prograde-ness, as other five parameters depend on that fact.
    double inclination = keplerianElements( inclinationIndex );

    // If inclination is outside range [0,PI].
    if ( ( inclination < 0.0 ) || ( inclination > PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "Inclination is expected in range [0," << PI << "]\n"
                     << "Specified inclination: " << inclination << " rad." << std::endl;

        // Throw exception.
        boost::throw_exception( std::runtime_error( errorMessage.str( ) ) );
    }
    //Else, nothing wrong and continue.

    // Compute set dependant helper parameters.
    double argumentOfPeriapsisAndAscendingNode = 0.0;
    double tangentOfHalfInclination = 0.0;

    // If normal set of equations is to be used (i.e. not inverse set for singularity).
    if ( !avoidSingularityAtPiInclination )
    {
        // Add (+) argument of periapsis and longitude of ascending node.
        argumentOfPeriapsisAndAscendingNode = keplerianElements( argumentOfPeriapsisIndex )
                + keplerianElements( longitudeOfAscendingNodeIndex );

        // Take the tangent of half the inclination.
        tangentOfHalfInclination = std::tan( inclination / 2.0 );
    }
    else // the orbit is retrograde.
    {
        // Subtract (-) longitude of ascending node from argument of periapsis.
        argumentOfPeriapsisAndAscendingNode = keplerianElements( argumentOfPeriapsisIndex )
                - keplerianElements( longitudeOfAscendingNodeIndex );

        // Take the inverse of the tangent of half the inclination to avoid singularity at tan
        // PI/2.
        tangentOfHalfInclination = 1.0 / std::tan( inclination / 2.0 );
    }

    // Compute f-element.
    modifiedEquinoctialState( fElementIndex ) = eccentricity
            * std::cos( argumentOfPeriapsisAndAscendingNode );

    // Compute g-element.
    modifiedEquinoctialState( gElementIndex ) = eccentricity
            * std::sin( argumentOfPeriapsisAndAscendingNode );

    // Compute h-element.
    modifiedEquinoctialState( hElementIndex ) = tangentOfHalfInclination
            * std::cos( keplerianElements( longitudeOfAscendingNodeIndex ) );

    // Compute k-element.
    modifiedEquinoctialState( kElementIndex ) = tangentOfHalfInclination
            * std::sin( keplerianElements( longitudeOfAscendingNodeIndex ) );

    // Compute true longitude (modulo 2 PI to keep within interval -2PI to 2PI).
    modifiedEquinoctialState( trueLongitudeIndex )
            = tudat::basic_mathematics::computeModulo(
                argumentOfPeriapsisAndAscendingNode
                + keplerianElements( trueAnomalyIndex ), 2.0 * PI );

    // Give back result.
    return modifiedEquinoctialState;
}

//! Convert modified equinoctial to Keplerian orbital elements.
tudat::basic_mathematics::Vector6d convertModifiedEquinoctialToKeplerianElements(
        const tudat::basic_mathematics::Vector6d& modifiedEquinoctialElements,
        const bool avoidSingularityAtPiInclination )
// Using unknown source pdf, code archive E. Heeren and personal derivation based on Hintz 2008.
{
    using tudat::basic_mathematics::mathematical_constants::PI;

    // Declaration of output vector.
    tudat::basic_mathematics::Vector6d convertedKeplerianElements = tudat::basic_mathematics::
            Vector6d::Zero( 6 );

    // for ease of referencing, almost all modified equinoctial elements.
    double fElement = modifiedEquinoctialElements( fElementIndex );
    double gElement = modifiedEquinoctialElements( gElementIndex );
    double hElement = modifiedEquinoctialElements( hElementIndex );
    double kElement = modifiedEquinoctialElements( kElementIndex );

    // Tolerance for singularities.
    double singularityTolerance = 1.0e-15;

    // Compute eccentricity.
    double eccentricity = std::sqrt( fElement * fElement + gElement * gElement );
    convertedKeplerianElements( eccentricityIndex ) = eccentricity;

    // Compute semi-major axis.
    // If eccentricity is not near-parabolic.
    if ( std::fabs( eccentricity - 1.0 ) > singularityTolerance )
    {
        // Use semi-latus rectum and eccentricity to calculate semi-major axis with a=p/(1-e^2).
        convertedKeplerianElements( semiMajorAxisIndex ) =
                modifiedEquinoctialElements( semiLatusRectumIndex ) /
                ( 1.0 - eccentricity * eccentricity );
    }
    // This is (almost) a parabola and the semimajor axis will tend to infinity and be useless.
    else
    {
        // Give semi-latus rectum instead.
        convertedKeplerianElements( semiLatusRectumIndex ) =
                modifiedEquinoctialElements( semiLatusRectumIndex );
    }

    // Compute longitude of ascending node.

    // Compute solution in [-PI,PI] interval with atan2.
    double longitudeOfAscendingNode = std::atan2( kElement, hElement );

    // Store longitude of ascending node.
    convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = tudat::basic_mathematics::
            computeModulo( longitudeOfAscendingNode, 2.0 * PI );

    // Compute inclination.

    // If prograde factor must be positive, otherwise negative.
    double retrogradeFactor = avoidSingularityAtPiInclination ? -1.0 : 1.0;

    // Calculate magnitude of inclination with retrogradeFactor (derived based on Hintz, 2008).
    double hSquaredPlusKSquared = std::pow( hElement * hElement + kElement * kElement ,
                                            retrogradeFactor );

    convertedKeplerianElements( inclinationIndex ) =
            2.0 * std::atan( std::sqrt( hSquaredPlusKSquared ) );
    // Was: std::atan2(2.0 * std::sqrt(hSquaredPlusKSquared) , (1.0 - hSquaredPlusKSquared));

    // Compute argument of periapsis.

    // Helper quantity (omega + I * OMEGA ), for argument of periapsis and true anomaly.
    double argumentOfPeriapsisAndLongitude;

    // If eccentricity is (near) circular
    if ( eccentricity < singularityTolerance )
        // Then the composite argument of periapsis and longitude together cannot be determined
        // because both arguments of the atan2 functions will be zero.
    {
        // Set argument of periapsis to zero.
        convertedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;

        // Set composite argument to longitude only, since argument of periapsis is now zero.
        argumentOfPeriapsisAndLongitude = retrogradeFactor * longitudeOfAscendingNode;
    }
    else
        // Argument of periapsis can be found from composite argument.
    {
        // Compute composite argument.
        argumentOfPeriapsisAndLongitude = std::atan2( gElement, fElement );

        // Compute argument of periapsis.
        convertedKeplerianElements( argumentOfPeriapsisIndex ) = tudat::basic_mathematics::
                computeModulo( argumentOfPeriapsisAndLongitude
                               - retrogradeFactor * longitudeOfAscendingNode, 2.0 * PI );
    }

    // Compute true anomaly.
    convertedKeplerianElements( trueAnomalyIndex ) = tudat::basic_mathematics::
            computeModulo( modifiedEquinoctialElements( trueLongitudeIndex )
                           - argumentOfPeriapsisAndLongitude, 2.0 * PI );

    // Return converted elements.
    return convertedKeplerianElements;
}

//! Convert Cartesian to modified equinoctial orbital elements using implicit MEE equation set.
tudat::basic_mathematics::Vector6d convertCartesianToModifiedEquinoctialElements(
        const tudat::basic_mathematics::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter )
{
    using basic_mathematics::mathematical_constants::PI;

    // Convert to keplerian elements.
    tudat::basic_mathematics::Vector6d keplerianElements = convertCartesianToKeplerianElements(
                cartesianElements, centralBodyGravitationalParameter );

    // Check whether orbit is retrograde.
    bool avoidSingularityAtPiInclination =
            mission_geometry::isOrbitRetrograde( keplerianElements );

    // Convert to modified equinoctial elements.
    return convertKeplerianToModifiedEquinoctialElements(
                keplerianElements, avoidSingularityAtPiInclination );
}

//! Convert Cartesian to modified equinoctial orbital elements using explicit MEE equation set.
tudat::basic_mathematics::Vector6d convertCartesianToModifiedEquinoctialElements(
        const tudat::basic_mathematics::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter,
        const bool avoidSingularityAtPiInclination )
{
    // Convert Cartesian to Keplerian to modified equinoctial elements.
    return convertKeplerianToModifiedEquinoctialElements(
                convertCartesianToKeplerianElements(
                    cartesianElements, centralBodyGravitationalParameter ),
                avoidSingularityAtPiInclination );
}

//! Convert Modified Equinoctial Elements to Cartesian Elements.
tudat::basic_mathematics::Vector6d convertModifiedEquinoctialToCartesianElements(
        const tudat::basic_mathematics::Vector6d& modifiedEquinoctialElements,
        const double centralBodyGravitationalParameter,
        const bool avoidSingularityAtPiInclination )
// Using unnamed pdf and code archive Bart Rˆmgens.
{
    // Creating output vector.
    tudat::basic_mathematics::Vector6d convertedCartesianElements = tudat::basic_mathematics::
            Vector6d::Zero( 6 );

    // If the prograde equations are to be used.
    if ( !avoidSingularityAtPiInclination )
    {
        // Local storage of elements for ease of access.
        double semiLatusRectum = modifiedEquinoctialElements( semiLatusRectumIndex );
        double fElement = modifiedEquinoctialElements( fElementIndex );
        double gElement = modifiedEquinoctialElements( gElementIndex );
        double hElement = modifiedEquinoctialElements( hElementIndex );
        double kElement = modifiedEquinoctialElements( kElementIndex );
        double trueLongitude = modifiedEquinoctialElements( trueLongitudeIndex );

        // Computing helper parameters.

        // Square-root
        double squareRootOfGravitationalParameterOverSemiLatusRectum =
                std::sqrt( centralBodyGravitationalParameter / semiLatusRectum );

        // Cosines and sine.
        double cosineTrueLongitude = std::cos( trueLongitude );
        double sineTrueLongitude = std::sin( trueLongitude );

        // s≤, a≤, w, r from references.
        double sSquaredParameter = 1.0 + hElement * hElement + kElement * kElement;
        double aSquaredParameter = hElement * hElement - kElement*kElement;
        double wParameter = 1.0 + fElement * cosineTrueLongitude
                + gElement * sineTrueLongitude;
        double radius = semiLatusRectum / wParameter;

        // Computing position and storing (using code archive Bart Rˆmgens and unnamed pdf).
        convertedCartesianElements( xPositionIndex )
                = radius / sSquaredParameter * ( cosineTrueLongitude + aSquaredParameter
                                                 * cosineTrueLongitude + 2.0 * hElement * kElement
                                                 * sineTrueLongitude );
        convertedCartesianElements( yPositionIndex )
                = radius / sSquaredParameter * ( sineTrueLongitude - aSquaredParameter
                                                 * sineTrueLongitude + 2.0 * hElement * kElement
                                                 * cosineTrueLongitude );
        convertedCartesianElements( zPositionIndex ) = 2.0 * radius / sSquaredParameter
                * ( hElement * sineTrueLongitude - kElement * cosineTrueLongitude );

        // Computing velocities (these can probably use more optimal computation).
        convertedCartesianElements( xVelocityIndex ) = -1.0 / sSquaredParameter
                * squareRootOfGravitationalParameterOverSemiLatusRectum
                * ( sineTrueLongitude + aSquaredParameter * sineTrueLongitude
                    - 2.0 * hElement * kElement * cosineTrueLongitude + gElement
                    - 2.0 * fElement * hElement * kElement + aSquaredParameter * gElement );
        convertedCartesianElements( yVelocityIndex ) = -1.0 / sSquaredParameter
                * squareRootOfGravitationalParameterOverSemiLatusRectum
                * ( -cosineTrueLongitude + aSquaredParameter * cosineTrueLongitude
                    + 2.0 * hElement * kElement * sineTrueLongitude - fElement
                    + 2.0 * gElement * hElement * kElement + aSquaredParameter * fElement );
        convertedCartesianElements( zVelocityIndex ) = 2.0 / sSquaredParameter
                * squareRootOfGravitationalParameterOverSemiLatusRectum
                * ( hElement * cosineTrueLongitude + kElement * sineTrueLongitude
                    + fElement * hElement + gElement * kElement );
    }
    else // use the indirect transformation via an intermediate Kepler state.
    {
        // Use the Keplerian state as an intermediary step, as there is no direct transformation
        // available in literature (personal derivation forgone for now due to time constraints).
        convertedCartesianElements =
                convertKeplerianToCartesianElements(
                    convertModifiedEquinoctialToKeplerianElements(
                        modifiedEquinoctialElements, avoidSingularityAtPiInclination ),
                    centralBodyGravitationalParameter );
    }

    // Return converted set of elements.
    return convertedCartesianElements;
}

} // namespace orbital_element_conversions
} // namespace basic_astrodynamics
} // namespace tudat
