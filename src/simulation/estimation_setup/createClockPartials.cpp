#include "tudat/simulation/estimation_setup/createClockPartials.h"

namespace tudat
{

namespace observation_partials
{

std::map< int, double > getTimingPartialMultipliers( const observation_models::ObservableType observable )
{
    std::map< int, double > multipliers;
    switch( observable )
    {
    case observation_models::one_way_range:
        multipliers[ 1 ] = 1.0;
        multipliers[ 0 ] = -1.0;
        break;
    case observation_models::n_way_range:
        multipliers[ 3 ] = 1.0;
        multipliers[ 2 ] = -1.0;
        multipliers[ 1 ] = 1.0;
        multipliers[ 0 ] = -1.0;
        break;
    default:
        throw std::runtime_error( "Error, observable " + observation_models::getObservableName( observable ) +
            " not implemented when getting timing partial multipliers" );

    }
    return multipliers;

}

}

}

