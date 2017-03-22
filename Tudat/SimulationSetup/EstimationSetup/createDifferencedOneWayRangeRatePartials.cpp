#include "Tudat/SimulationSetup/EstimationSetup/createDifferencedOneWayRangeRatePartials.h"

namespace tudat
{

namespace observation_partials
{

typedef std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > > SingleLinkObservationPartialList;
void removeObservationBiasPartialFromList( SingleLinkObservationPartialList& observationPartialList )
{
    std::vector< std::pair< int, int > > entriesToDelete;
    for( SingleLinkObservationPartialList::iterator partialIterator = observationPartialList.begin( );
         partialIterator != observationPartialList.end( ); partialIterator++ )
    {
        if( false )
                //boost::dynamic_pointer_cast< ObservationPartialWrtBias< 1 > >( partialIterator->second ) )
        {
            entriesToDelete.push_back( partialIterator->first );
        }
    }

    for( unsigned int i = 0; i < entriesToDelete.size( ); i++ )
    {
        observationPartialList.erase( entriesToDelete[ i ] );
    }
}

std::pair< PerLinkEndPerLightTimeSolutionCorrections, PerLinkEndPerLightTimeSolutionCorrections > splitOneWayRangeRateLightTimeCorrectionsBetweenArcs(
        const PerLinkEndPerLightTimeSolutionCorrections& combinedCorrections )
{
    PerLinkEndPerLightTimeSolutionCorrections arcStartCorrections;
    PerLinkEndPerLightTimeSolutionCorrections arcEndCorrections;

    std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > currentArcStartCorrections;
    std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > currentArcEndCorrections;

    for( PerLinkEndPerLightTimeSolutionCorrections::const_iterator correctionIterator = combinedCorrections.begin( );
         correctionIterator != combinedCorrections.end( ); correctionIterator++ )
    {
        currentArcStartCorrections.clear( );
        currentArcEndCorrections.clear( );
        if( correctionIterator->second.size( ) > 0 )
        {
            if( correctionIterator->second.size( ) != 2 )
            {
                std::cerr<<"Error when splitting one-way range rate light time corrections, size is "<<correctionIterator->second.size( )<<std::endl;
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
