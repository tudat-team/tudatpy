#ifndef DIFFERENCEDONEWAYRANGERATEPARTIAL_H
#define DIFFERENCEDONEWAYRANGERATEPARTIAL_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayRangePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"

namespace tudat
{

namespace observation_partials
{

class OneWayRangeRateScaling: public PositionPartialScaling
{
public:
    OneWayRangeRateScaling(
            const boost::shared_ptr< OneWayRangeScaling > oneWayRangeScalerArcStart,
            const boost::shared_ptr< OneWayRangeScaling > oneWayRangeScalerArcEnd ):
        oneWayRangeScalerArcStart_( oneWayRangeScalerArcStart ), oneWayRangeScalerArcEnd_( oneWayRangeScalerArcEnd ){ }

    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd )
    {
        oneWayRangeScalerArcStart_->update(
                    std::vector< Eigen::Vector6d >(
                        linkEndStates.begin( ), linkEndStates.begin( ) + 2 ), times, fixedLinkEnd );
        oneWayRangeScalerArcEnd_->update(
                    std::vector< Eigen::Vector6d >(
                        linkEndStates.begin( ) + 2, linkEndStates.begin( ) + 4 ), times, fixedLinkEnd );
    }

    boost::shared_ptr< OneWayRangeScaling > oneWayRangeScalerArcStart_;
    boost::shared_ptr< OneWayRangeScaling > oneWayRangeScalerArcEnd_;
};

class DifferencedOneWayRangeRatePartial: public ObservationPartial< 1 >
{

public:

    DifferencedOneWayRangeRatePartial(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const boost::shared_ptr< ObservationPartial< 1 > > arcStartRangePartial,
            const boost::shared_ptr< ObservationPartial< 1 > > arcEndRangePartial):
        ObservationPartial< 1 >( parameterIdentifier ),
        arcStartRangePartial_( arcStartRangePartial ),
        arcEndRangePartial_( arcEndRangePartial ){ }

    ~DifferencedOneWayRangeRatePartial( ) { }

    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Zero( ) );

protected:

    boost::shared_ptr< ObservationPartial< 1 > > arcStartRangePartial_;

    boost::shared_ptr< ObservationPartial< 1 > > arcEndRangePartial_;
};

}

}


#endif // DIFFERENCEDONEWAYRANGERATEPARTIAL_H
