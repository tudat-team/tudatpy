#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/mission_segments/lambertTargeter.h"

namespace tudat
{

namespace mission_segments
{

Eigen::Vector6d getLambertTargeterInitialKeplerianState(
        const LambertTargeter& lambertTargeter )
{
    return orbital_element_conversions::convertCartesianToKeplerianElements(
                lambertTargeter.getDepartureState( ), lambertTargeter.getCentralBodyGravitationalParameter( ) );
}

Eigen::Vector6d getLambertTargeterKeplerianStateDuringTransfer(
        const LambertTargeter& lambertTargeter,
        const double timeAfterDeparture )
{
    Eigen::Vector6d initialKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements(
                lambertTargeter.getDepartureState( ), lambertTargeter.getCentralBodyGravitationalParameter( ) );
    return orbital_element_conversions::propagateKeplerOrbit(
                initialKeplerianState, timeAfterDeparture, lambertTargeter.getCentralBodyGravitationalParameter( ) );
}

Eigen::Vector6d getLambertTargeterCartesianStateDuringTransfer(
        const LambertTargeter& lambertTargeter,
        const double timeAfterDeparture )
{
    return orbital_element_conversions::convertKeplerianToCartesianElements(
                getLambertTargeterKeplerianStateDuringTransfer( lambertTargeter, timeAfterDeparture ),
                lambertTargeter.getCentralBodyGravitationalParameter( ) );
}

} // namespace mission_segments

} // namespace tudat

