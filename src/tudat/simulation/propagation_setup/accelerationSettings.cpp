#include "tudat/simulation/propagation_setup/accelerationSettings.h"

namespace tudat
{

namespace simulation_setup
{
std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getExtendedSinglePointMassInteractions(
        const int maximumDegreeOfBodyUndergoingAcceleration,
        const int maximumOrderOfBodyUndergoingAcceleration,
        const int maximumDegreeOfBodyExertingAcceleration,
        const int maximumOrderOfBodyExertingAcceleration )
{
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse;

    for( int i = 1; i <= maximumDegreeOfBodyUndergoingAcceleration; i++ )
    {
        for( int j = 0; ( j <= maximumOrderOfBodyUndergoingAcceleration && j <= i ); j++ )
        {
            coefficientCombinationsToUse.push_back( std::make_tuple( i, j, 0, 0 ) );
        }
    }

    for( int k = 1; k <= maximumDegreeOfBodyExertingAcceleration; k++ )
    {
        for( int l = 0; ( l <= maximumOrderOfBodyExertingAcceleration && l <= k ); l++ )
        {
            coefficientCombinationsToUse.push_back( std::make_tuple( 0, 0, k, l ) );
        }
    }

    coefficientCombinationsToUse.push_back( std::make_tuple( 0, 0, 0, 0 ) );

    return coefficientCombinationsToUse;
}

} // namespace simulation_setup

} // namespace tudat
