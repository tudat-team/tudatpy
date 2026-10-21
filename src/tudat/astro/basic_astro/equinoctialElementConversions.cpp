/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 */

#include "tudat/astro/basic_astro/equinoctialElementConversions.h"

#include <cmath>
#include <stdexcept>

#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace orbital_element_conversions
{

Eigen::Vector6d convertOrbitalElementsWithMeanMotionAndMeanAnomalyToEquinoctialElements( const Eigen::Vector6d& orbitalElements,
                                                                                         const bool useRetrogradeElements )
{
    if( !orbitalElements.allFinite( ) || orbitalElements( 0 ) <= 0.0 || orbitalElements( 1 ) < 0.0 )
    {
        throw std::invalid_argument( "Orbital elements must be finite, with positive mean motion and non-negative eccentricity." );
    }

    const double meanMotion = orbitalElements( 0 );
    const double eccentricity = orbitalElements( 1 );
    const double inclination = orbitalElements( 2 );
    const double ascendingNode = orbitalElements( 3 );
    const double argumentOfPeriapsis = orbitalElements( 4 );
    const double meanAnomaly = orbitalElements( 5 );

    // The retrograde branch replaces tan(i/2) by cot(i/2) and uses the corresponding longitude of perigee.
    const double inclinationParameter = useRetrogradeElements ? 1.0 / std::tan( inclination / 2.0 ) : std::tan( inclination / 2.0 );
    const double longitudeOfPeriapsis = useRetrogradeElements ? argumentOfPeriapsis - ascendingNode : argumentOfPeriapsis + ascendingNode;

    Eigen::Vector6d equinoctialElements;
    equinoctialElements << std::log( meanMotion ), eccentricity * std::cos( longitudeOfPeriapsis ),
            eccentricity * std::sin( longitudeOfPeriapsis ), inclinationParameter * std::cos( ascendingNode ),
            inclinationParameter * std::sin( ascendingNode ), longitudeOfPeriapsis + meanAnomaly;
    return equinoctialElements;
}

Eigen::Vector6d convertEquinoctialToOrbitalElementsWithMeanMotionAndMeanAnomaly( const Eigen::Vector6d& equinoctialElements,
                                                                                 const bool useRetrogradeElements )
{
    if( !equinoctialElements.allFinite( ) )
    {
        throw std::invalid_argument( "Equinoctial orbital parameters must be finite." );
    }

    const double twoPi = 2.0 * mathematical_constants::PI;
    const double eccentricity = std::hypot( equinoctialElements( 1 ), equinoctialElements( 2 ) );
    const double inclinationParameterNorm = std::hypot( equinoctialElements( 3 ), equinoctialElements( 4 ) );
    const double ascendingNode =
            basic_mathematics::computeModulo( std::atan2( equinoctialElements( 4 ), equinoctialElements( 3 ) ), twoPi );
    const double longitudeOfPeriapsis =
            basic_mathematics::computeModulo( std::atan2( equinoctialElements( 2 ), equinoctialElements( 1 ) ), twoPi );
    const double inclination =
            useRetrogradeElements ? 2.0 * std::atan2( 1.0, inclinationParameterNorm ) : 2.0 * std::atan( inclinationParameterNorm );
    const double argumentOfPeriapsis = basic_mathematics::computeModulo(
            useRetrogradeElements ? longitudeOfPeriapsis + ascendingNode : longitudeOfPeriapsis - ascendingNode, twoPi );
    const double meanAnomaly = basic_mathematics::computeModulo( equinoctialElements( 5 ) - longitudeOfPeriapsis, twoPi );

    Eigen::Vector6d orbitalElements;
    orbitalElements << std::exp( equinoctialElements( 0 ) ), eccentricity, inclination, ascendingNode, argumentOfPeriapsis, meanAnomaly;
    return orbitalElements;
}

}  // namespace orbital_element_conversions
}  // namespace tudat
