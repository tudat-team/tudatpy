#include <boost/assign/list_of.hpp>

#include "Astrodynamics/ObservationModels/observationBias.h"
#include "Astrodynamics/Bodies/getHardwareModels.h"

namespace tudat
{

namespace observation_models
{

std::map< ObservableType, std::map< LinkEnds, boost::shared_ptr< ObservationBiasInterface > > > observationBiasInterfaceList;

Eigen::VectorXd RangeObservationBiasInterface::getObservationBias( const std::vector< double >& observationTimes )
{
    double timingError = 0.0;

    for( timingSystemIterator_ = timingSystems_.begin( ); timingSystemIterator_ != timingSystems_.end( ); timingSystemIterator_++ )
    {
        double currentTimingError = 0.0;
        if( useStochasticError_ )
        {
            currentTimingError = timingSystemIterator_->second->getCompleteClockError( observationTimes[ timingSystemIterator_->first ] );
        }
        else
        {
            currentTimingError = timingSystemIterator_->second->getDeterministicClockError( observationTimes[ timingSystemIterator_->first ] );
        }

        if( addContributions_[ timingSystemIterator_->first ] == 1 )
        {
            timingError += currentTimingError;
        }
        else
        {
            timingError -= currentTimingError;
        }
    }

    Eigen::VectorXd totalBias = observationBias_;
    totalBias.x( ) += physical_constants::SPEED_OF_LIGHT * timingError;

    return totalBias;
}


void setObservationBiases( const std::map< ObservableType, std::map< LinkEnds, double > >& biasesForDoubleObservables )
{
    for( std::map< ObservableType, std::map< LinkEnds, double > >::const_iterator observableIterator = biasesForDoubleObservables.begin( );
         observableIterator != biasesForDoubleObservables.end( ); observableIterator++ )
    {
        for( std::map< LinkEnds, double >::const_iterator linkIterator = observableIterator->second.begin( );
             linkIterator != observableIterator->second.end( ); linkIterator++ )
        {
            observationBiasInterfaceList[ observableIterator->first ][ linkIterator->first ] =
                    boost::make_shared< ObservationBiasInterface >( Eigen::VectorXd::Constant( 1, linkIterator->second ) );
        }
    }
}

void setObservationBiases( const std::map< ObservableType, std::vector< LinkEnds > >& linkEndList,
                           const std::map< ObservableType, double >& biasPerObservables )
{
    for( std::map< ObservableType, std::vector< LinkEnds > >::const_iterator observableIterator = linkEndList.begin( );
         observableIterator != linkEndList.end( ); observableIterator++ )
    {
        for( unsigned int i = 0; i < observableIterator->second.size( ); i++ )
        {
            observationBiasInterfaceList[ observableIterator->first ][ observableIterator->second[ i ] ] =
                    boost::make_shared< ObservationBiasInterface >(
                        Eigen::VectorXd::Constant( 1, biasPerObservables.at( observableIterator->first ) ) );
        }
    }
}

void setObservationBiases( const std::map< ObservableType, std::vector< LinkEnds > >& linkEndList,
                           const double biasForAllObservations )
{
    for( std::map< ObservableType, std::vector< LinkEnds > >::const_iterator observableIterator = linkEndList.begin( );
         observableIterator != linkEndList.end( ); observableIterator++ )
    {
        for( unsigned int i = 0; i < observableIterator->second.size( ); i++ )
        {
            observationBiasInterfaceList[ observableIterator->first ][ observableIterator->second[ i ] ] =
                    boost::make_shared< ObservationBiasInterface >( Eigen::VectorXd::Constant( 1, biasForAllObservations ) );
        }
    }
}

std::vector< bool > getClockErrorAdditions( const ObservableType observableType, const LinkEnds& linkEnds )
{
    std::vector< bool > clockAdditionVector;
    switch( observableType )
    {
    case oneWayRange:
    {
        clockAdditionVector.push_back( 0 );
        clockAdditionVector.push_back( 1 );
        break;
    }
//    case twoWayRange:
//    {
//        clockAdditionVector.push_back( 0 );
//        clockAdditionVector.push_back( 1 );
//        clockAdditionVector.push_back( 0 );
//        break;
//    }
    case nWayRange:
    {
        clockAdditionVector.push_back( 0 );
        for( unsigned int i = 0; i < linkEnds.size( ) - 2; i++ )
        {
            clockAdditionVector.push_back( 0 );
            clockAdditionVector.push_back( 1 );
        }
        clockAdditionVector.push_back( 1 );
        break;
    }
    default:
        std::cerr<<"Error, clock addition vector not found for observable "<<observableType<<std::endl;
    }

    return clockAdditionVector;
}

void setObservationBiasesWithTimingErrors( const std::map< ObservableType, std::vector< LinkEnds > >& linkEndList,
                                           const NamedBodyMap& bodyMap,
                                           const bool includeStochasticNoise,
                                           const std::map< ObservableType, std::map< LinkEnds, double > >& biasesForDoubleObservables )
{
    for( std::map< ObservableType, std::vector< LinkEnds > >::const_iterator observableIterator = linkEndList.begin( );
         observableIterator != linkEndList.end( ); observableIterator++ )
    {
        std::map< LinkEnds, double > currentObservableTypeConstantBiases;
        if( biasesForDoubleObservables.count( observableIterator->first ) > 0 )
        {
            currentObservableTypeConstantBiases = biasesForDoubleObservables.at( observableIterator->first );
        }

        for( unsigned int i = 0; i < observableIterator->second.size( ); i++ )
        {
            double currentConstantBias = 0.0;
            if( currentObservableTypeConstantBiases.count( observableIterator->second[ i ] ) > 0 )
            {
                currentConstantBias = currentObservableTypeConstantBiases[ observableIterator->second[ i ] ];
            }

            std::map< int, boost::shared_ptr< hardware_models::TimingSystem > > linkEndTimingSystems;

            for( LinkEnds::const_iterator linkEndIterator = observableIterator->second[ i ].begin( );
                 linkEndIterator != observableIterator->second[ i ].end( ); linkEndIterator++ )
            {

                boost::shared_ptr< hardware_models::TimingSystem > currentTimingSystem =
                        bodies::getTimingSystem( linkEndIterator->second, bodyMap );
                if( currentTimingSystem != NULL )
                {
                    std::vector< int > currentLinkEndIndices = getLinkEndIndicesForLinkEndTypeAtObservable(
                                observableIterator->first, linkEndIterator->first, observableIterator->second[ i ].size( ) );
                    for( unsigned int j = 0; j < currentLinkEndIndices.size( ); j++ )
                    {
                        if( linkEndTimingSystems.count( currentLinkEndIndices.at( j ) ) > 0 )
                        {
                            std::cerr<<"Error when setting timing error observation bias for type "<<observableIterator->first<<
                                       ", overriding timing system with link end "<<linkEndIterator->first<<std::endl;
                        }
                        else
                        {
                            linkEndTimingSystems[ currentLinkEndIndices.at( j ) ] = currentTimingSystem;
                        }
                    }
                }

            }

            observationBiasInterfaceList[ observableIterator->first ][ observableIterator->second[ i ] ] =
                    boost::make_shared< RangeObservationBiasInterface >(
                        linkEndTimingSystems, getClockErrorAdditions( observableIterator->first, observableIterator->second.at( i ) ),
                        includeStochasticNoise, Eigen::VectorXd::Constant( 1, currentConstantBias ) );

        }
    }

}

}

}
