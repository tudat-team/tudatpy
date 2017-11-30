#include "Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h"

namespace tudat
{

namespace simulation_setup
{
std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getExtendedSinglePointMassInteractions(
        const int maximumDegreeOfBodyUndergoingAcceleration,
        const int maximumOrderOfBodyUndergoingAcceleration,
        const int maximumDegreeOfBodyExertingAcceleration,
        const int maximumOrderOfBodyExertingAcceleration )
{
    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse;

    for( unsigned int i = 1; i <= maximumDegreeOfBodyUndergoingAcceleration; i++ )
    {
        for( unsigned int j = 0; ( j <= maximumOrderOfBodyUndergoingAcceleration && j <= i ); j++ )
        {
            coefficientCombinationsToUse.push_back( boost::make_tuple( i, j, 0, 0 ) );
        }
    }

    for( unsigned int k = 1; k <= maximumDegreeOfBodyExertingAcceleration; k++ )
    {
        for( unsigned int l = 0; ( l <= maximumOrderOfBodyExertingAcceleration && l <= k ); l++ )
        {
            coefficientCombinationsToUse.push_back( boost::make_tuple( 0, 0, k, l ) );
        }
    }

    coefficientCombinationsToUse.push_back( boost::make_tuple( 0, 0, 0, 0 ) );

    return coefficientCombinationsToUse;
}

} // namespace simulation_setup

} // namespace tudat
