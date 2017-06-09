#ifndef NWAYRANGEPARTIAL_H
#define NWAYRANGEPARTIAL_H

#include <iostream>

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayRangePartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{

class NWayRangeScaling: public PositionPartialScaling
{
public:
    NWayRangeScaling( const std::map< int, boost::shared_ptr< OneWayRangeScaling > >& constituentRangeScalings ):
        constituentRangeScalings_( constituentRangeScalings )
    {
        numberOfLinkEnds_ = constituentRangeScalings.size( ) + 1;
    }

    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation );

    double getProjectedRelativeVelocityRatio( const int linkIndex )
    {
        return projectedRelativeVelocityRatios_.at( linkIndex );
    }

private:

    std::map< int, boost::shared_ptr< OneWayRangeScaling > > constituentRangeScalings_;

    std::map< int, boost::shared_ptr< OneWayRangeScaling > >::iterator constituentRangeScalingIterator_;

    std::map< int, double > projectedRelativeVelocityRatios_;

    int numberOfLinkEnds_;

};

class NWayRangePartial: public ObservationPartial< 1 >
{

public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > NWayRangePartialReturnType;

    NWayRangePartial( const boost::shared_ptr< NWayRangeScaling > nWayRangeScaler,
                      const std::map< int, boost::shared_ptr< ObservationPartial< 1 > > >& rangePartialList,
                      const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
                      const int numberOfLinkEnds ):
        ObservationPartial< 1 >( parameterIdentifier ), nWayRangeScaler_( nWayRangeScaler ), rangePartialList_( rangePartialList ),
        numberOfLinkEnds_( numberOfLinkEnds ){ }

    ~NWayRangePartial( ) { }

    NWayRangePartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector1d& currentObservation = Eigen::Vector1d::Constant( TUDAT_NAN ) );

    estimatable_parameters::EstimatebleParameterIdentifier getParameterIdentifier( )
    {
        return parameterIdentifier_;
    }

    boost::shared_ptr< NWayRangeScaling > getNWayRangeScaler( )
    {
        return nWayRangeScaler_;
    }


protected:

    boost::shared_ptr< NWayRangeScaling > nWayRangeScaler_;

    std::map< int, boost::shared_ptr< ObservationPartial< 1 > > > rangePartialList_;

    std::map< int, boost::shared_ptr< ObservationPartial< 1 > > >::iterator rangePartialIterator_;

    estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier_;

    int minimumDependentLinkEnd_;

    int maximumDependentLinkEnd_;

    int numberOfLinkEnds_;
};

}

}


#endif // NWAYRANGEPARTIAL_H
