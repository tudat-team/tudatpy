#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tudat
{

namespace simulation_setup
{

int noiseSeed = 0;

std::function< Eigen::VectorXd( const double ) > getNoiseFunctionForObservable(
        const std::function< double( const double ) > singleNoiseFunction,
        const observation_models::ObservableType observableType )
{
    int observableSize = observation_models::getObservableSize( observableType );
    std::function< Eigen::VectorXd( const double ) > noiseFunction;
    if( observableSize != 1 )
    {
        noiseFunction = [=]( const double time )
        {
            Eigen::VectorXd noise = Eigen::VectorXd::Zero( observableSize );
            for( int i = 0; i < observableSize; i++ )
            {
                noise( i ) = singleNoiseFunction( time );
            }
            return noise;
        };
    }
    else
    {
        noiseFunction = [=]( const double time ){ return ( Eigen::VectorXd( 1 )<<singleNoiseFunction( time ) ).finished( ); };
    }
    return noiseFunction;
}


Eigen::VectorXd getIdenticallyAndIndependentlyDistributedNoise(
        const std::function< double( const double ) > noiseFunction,
        const int observationSize,
        const double evaluationTime )
{
    Eigen::VectorXd noiseValues = Eigen::VectorXd( observationSize );
    for( int i = 0; i < observationSize; i++ )
    {
        noiseValues( i ) = noiseFunction( evaluationTime );
    }
    return noiseValues;
}



}

}
