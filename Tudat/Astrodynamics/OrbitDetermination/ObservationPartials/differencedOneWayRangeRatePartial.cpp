#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/differencedOneWayRangeRatePartial.h"

namespace tudat
{

namespace observation_partials
{

std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > DifferencedOneWayRangeRatePartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime )
{

    using namespace observation_partials;

    std::vector< Eigen::Vector6d > arcStartStates;
    arcStartStates.push_back( states.at( 0 ) );
    arcStartStates.push_back( states.at( 1 ) );
    std::vector< Eigen::Vector6d > arcEndStates;
    arcEndStates.push_back( states.at( 2 ) );
    arcEndStates.push_back( states.at( 3 ) );

    std::vector< double > arcStartTimes;
    arcStartTimes.push_back( times.at( 0 ) );
    arcStartTimes.push_back( times.at( 1 ) );
    std::vector< double > arcEndTimes;
    arcEndTimes.push_back( times.at( 2 ) );
    arcEndTimes.push_back( times.at( 3 ) );


    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > arcStartPartials =
            arcStartRangePartial_->calculatePartial( arcStartStates, arcStartTimes, linkEndOfFixedTime );

    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > arcEndPartials =
            arcEndRangePartial_->calculatePartial( arcEndStates, arcEndTimes, linkEndOfFixedTime );

    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > differencedRangeRatePartials;

    if( arcStartPartials.size( ) != arcEndPartials.size( ) )
    {
        std::cerr<<"Error when making differenced one way range rate partials, arc start and end partials inconsistent"<<std::endl;
    }

    //throw std::runtime_error( "Error, must compute arc duration" );
    double arcDuration = 60.0;

    for( unsigned int i = 0; i < arcStartPartials.size( ); i++ )
    {
        //std::cout<<arcStartPartials[ i ].first - arcEndPartials[ i ].first<<" "<<arcStartPartials[ i ].second -  arcEndPartials[ i ].second <<std::endl;
        differencedRangeRatePartials.push_back(
                    std::make_pair( -arcStartPartials[ i ].first / arcDuration, arcStartPartials[ i ].second ) );
    }

    for( unsigned int i = 0; i < arcEndPartials.size( ); i++ )
    {
        //std::cout<<arcEndPartials[ i ].first<<std::endl;
        differencedRangeRatePartials.push_back(
                    std::make_pair( arcEndPartials[ i ].first / arcDuration, arcEndPartials[ i ].second ) );
    }

    //std::cout<<( differencedRangeRatePartials[ 0 ].first+differencedRangeRatePartials[ 1 ].first ).cwiseQuotient( differencedRangeRatePartials[ 1 ].first  )<<std::endl;

    return differencedRangeRatePartials;
}

}

}

