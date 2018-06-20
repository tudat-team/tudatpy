/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/EstimationSetup/createDifferencedOneWayRangeRatePartials.h"

namespace tudat
{

namespace observation_partials
{

//! Function to split the total list of light-time corrections for one-way differenced range rate into list for either arc.
std::pair< PerLinkEndPerLightTimeSolutionCorrections, PerLinkEndPerLightTimeSolutionCorrections >
splitOneWayRangeRateLightTimeCorrectionsBetweenArcs(
        const PerLinkEndPerLightTimeSolutionCorrections& combinedCorrections )
{
    PerLinkEndPerLightTimeSolutionCorrections arcStartCorrections;
    PerLinkEndPerLightTimeSolutionCorrections arcEndCorrections;

    std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > currentArcStartCorrections;
    std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > currentArcEndCorrections;

    // Iterate over all link ends
    for( PerLinkEndPerLightTimeSolutionCorrections::const_iterator correctionIterator = combinedCorrections.begin( );
         correctionIterator != combinedCorrections.end( ); correctionIterator++ )
    {
        currentArcStartCorrections.clear( );
        currentArcEndCorrections.clear( );

        // if light-time corrections exist; split and put into separate lists.
        if( correctionIterator->second.size( ) > 0 )
        {
            if( correctionIterator->second.size( ) != 2 )
            {
               throw std::runtime_error(
                            "Error when splitting one-way range rate light time corrections, size is " +
                            std::to_string( correctionIterator->second.size( ) ) );
            }
            else
            {
                currentArcStartCorrections.push_back( correctionIterator->second.at( 0 ) );
                currentArcEndCorrections.push_back( correctionIterator->second.at( 1 ) );
            }
        }
        arcStartCorrections[ correctionIterator->first ] = currentArcStartCorrections;
        arcEndCorrections[ correctionIterator->first ] = currentArcEndCorrections;
    }
    return std::make_pair( arcStartCorrections, arcEndCorrections );
}

}

}
