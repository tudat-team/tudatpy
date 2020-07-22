#include <iostream>

#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/pso.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/bee_colony.hpp"
#include "pagmo/algorithms/sga.hpp"
#include "pagmo/algorithms/sea.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/bee_colony.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/nlopt.hpp"
#include "pagmo/algorithms/ihs.hpp"
#include "pagmo/algorithms/xnes.hpp"
#include "pagmo/algorithms/nsga2.hpp"
#include "pagmo/algorithms/moead.hpp"
#include "pagmo/algorithms/ihs.hpp"
#include "pagmo/algorithms/gaco.hpp"

pagmo::algorithm getMultiObjectiveAlgorithm( const int index )
{
    switch( index )
    {
    case 0:
    {
        pagmo::algorithm algo{ pagmo::nsga2( ) };
        return algo;
        break;
    }
    case 1:
    {
        pagmo::algorithm algo{ pagmo::moead( ) };
        return algo;
        break;
    }
    case 2:
    {
        pagmo::algorithm algo{ pagmo::ihs( ) };
        return algo;
        break;
    }
    default:
    {
        throw std::runtime_error( "Error, multi-objective pagmo algorithm " + std::to_string( index ) + " was not found." );
    }
    }
}
pagmo::algorithm getAlgorithm( const int index )
{
    switch( index )
    {
    case 0:
    {
        pagmo::algorithm algo{ pagmo::de( ) };
        return algo;
        break;
    }
    case 1:
    {
        pagmo::algorithm algo{ pagmo::sade( ) };
        return algo;
        break;
    }
    case 2:
    {
        pagmo::algorithm algo{ pagmo::de1220( ) };
        return algo;
        break;
    }
    case 3:
    {
        pagmo::algorithm algo{ pagmo::pso( ) };
        return algo;
        break;
    }
    case 4:
    {
        pagmo::algorithm algo{ pagmo::sea( ) };
        return algo;
        break;
    }
    case 5:
    {
        pagmo::algorithm algo{ pagmo::sga( ) };
        return algo;
        break;
    }
    case 6:
    {
        pagmo::algorithm algo{ pagmo::simulated_annealing( ) };
        return algo;
        break;
    }
    case 7:
    {
        pagmo::algorithm algo{ pagmo::bee_colony( ) };
        return algo;
        break;
    }
    case 8:
    {
        pagmo::algorithm algo{ pagmo::cmaes( ) };
        return algo;
        break;
    }
    case 9:
    {
        pagmo::algorithm algo{ pagmo::ihs( ) };
        return algo;
        break;
    }
    case 10:
    {
        pagmo::algorithm algo{ pagmo::xnes( ) };
        return algo;
        break;
    }
    default:
    {
        throw std::runtime_error( "Error, sinlge-objective pagmo algorithm " + std::to_string( index ) + " was not found." );
    }
    }
}

