#ifndef POSITIONOBSERVATIONMODEL_H
#define POSITIONOBSERVATIONMODEL_H

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/ObservationModels/observationModel.h"

namespace tudat
{

namespace observation_models
{


template< typename ObservationScalarType = double, typename TimeType = double >
class PositionObservationModel: public ObservationModel< 3, ObservationScalarType, TimeType, ObservationScalarType >
{
public:
    PositionObservationModel( const std::pair< std::string, std::string > bodyWithState,
                              const std::map< std::string, boost::shared_ptr< bodies::Body > > bodyMap ):
        ObservationModel< 3, ObservationScalarType, TimeType, ObservationScalarType >( positionObservable )
    {
        if( bodyWithState.second != "" )
        {
            std::cerr<<"Error when making state observation model, can only observe full body states"<<std::endl;
        }
        else
        {
            stateFunction_ = boost::bind(
                        &ephemerides::Ephemeris::getTemplatedStateFromEphemeris< ObservationScalarType, TimeType >,
                        bodyMap.at( bodyWithState.first )->getEphemeris( ), _1 );
        }
    }


    PositionObservationModel( const boost::function<  Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType ) > stateFunction ):
        stateFunction_( stateFunction ){ }

    ~PositionObservationModel( ) { }

    Eigen::Matrix< ObservationScalarType, 3, 1 > computeObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime = observed_body ) const
    {
        return stateFunction_( time ).segment( 0, 3 );
    }


    Eigen::Matrix< ObservationScalarType, 3, 1 > computeObservationsAndFullPrecisionLinkEndData(
                const TimeType time,
                const LinkEndType linkEndAssociatedWithTime,
                std::vector< TimeType >& linkEndTimes,
                std::vector< Eigen::Matrix< ObservationScalarType, 6, 1 > >& linkEndStates ) const
    {
        Eigen::Matrix< ObservationScalarType, 3, 1 >  observation = computeObservations( time, linkEndAssociatedWithTime );

        linkEndTimes.clear( );
        linkEndTimes.push_back( static_cast< Time >( time ) );

        linkEndStates.clear( );
        linkEndStates.push_back( stateFunction_( time ).template cast< ObservationScalarType >( ) );

        return observation;
    }


private:
    boost::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType ) > stateFunction_;
};

}

}

#endif // POSITIONOBSERVATIONMODEL_H
