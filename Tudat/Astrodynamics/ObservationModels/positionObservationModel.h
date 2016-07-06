#ifndef POSITIONOBSERVATIONMODEL_H
#define POSITIONOBSERVATIONMODEL_H

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"

namespace tudat
{

namespace observation_models
{


template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
class PositionObservationModel: public ObservationModel< 3, ObservationScalarType, TimeType, StateScalarType >
{
public:

    PositionObservationModel(
            const boost::function<  Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType& ) > stateFunction,
            const boost::shared_ptr< ObservationBias< 3 > > observationBiasCalculator = NULL ):
        ObservationModel< 3, ObservationScalarType, TimeType, StateScalarType >(
            position_observable, observationBiasCalculator ), stateFunction_( stateFunction ){ }

    ~PositionObservationModel( ) { }

    Eigen::Matrix< ObservationScalarType, 3, 1 > computeUnbiasedObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime = observed_body ) const
    {
        if( linkEndAssociatedWithTime != observed_body )
        {
            throw std::runtime_error(
                        "Error when computing position observable, associated link end must be observed_body " );
        }

        return stateFunction_( time ).segment( 0, 3 );
    }


    Eigen::Matrix< StateScalarType, 3, 1 > computeUnbiasedObservationsWithLinkEndData(
                const TimeType time,
                const LinkEndType linkEndAssociatedWithTime,
                std::vector< TimeType >& linkEndTimes,
                std::vector< Eigen::Matrix< StateScalarType, 6, 1 > >& linkEndStates ) const
    {
        if( linkEndAssociatedWithTime != observed_body )
        {
            throw std::runtime_error(
                        "Error when computing position observable, associated link end must be observed_body " );
        }

        Eigen::Matrix< ObservationScalarType, 3, 1 >  observation = stateFunction_( time ).segment( 0, 3 );

        linkEndTimes.clear( );
        linkEndTimes.push_back( static_cast< TimeType >( time ) );

        linkEndStates.clear( );
        linkEndStates.push_back( stateFunction_( time ).template cast< StateScalarType >( ) );

        return observation;
    }


private:
    boost::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType& ) > stateFunction_;
};

}

}

#endif // POSITIONOBSERVATIONMODEL_H
