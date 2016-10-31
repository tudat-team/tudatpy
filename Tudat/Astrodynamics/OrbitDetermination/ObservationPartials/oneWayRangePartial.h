#ifndef RANGEPARTIAL_H
#define RANGEPARTIAL_H

#include <iostream>

#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/positionPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/LightTimeCorrectionPartials/lightTimeCorrectionPartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{

class OneWayRangeScaling: public PositionPartialScaling
{
public:
    ~OneWayRangeScaling( ){ }

    void update( const std::vector< basic_mathematics::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd );

    Eigen::Matrix< double, 1, 3 > getScalingFactor( const observation_models::LinkEndType linkEndType, const observation_models::LinkEndType referenceTimeLinkEnd  );

private:
    Eigen::Matrix< double, 1, 3 > baseScalingFactor_;

    Eigen::Matrix< double, 1, 3 > transmitterReferenceScalingFactor_;

    Eigen::Matrix< double, 1, 3 > receiverReferenceScalingFactor_;

};

class OneWayRangePartial: public ObservationPartial< 1 >
{

public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > OneWayRangePartialReturnType;
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleOneWayRangePartialReturnType;

    OneWayRangePartial( const boost::shared_ptr< OneWayRangeScaling > oneWayRangeScaler,
                        const std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > >& positionPartialList,
                        const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
                        const std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >& lighTimeCorrectionPartials =
            std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) ):
        ObservationPartial< 1 >( parameterIdentifier ), oneWayRangeScaler_( oneWayRangeScaler ), positionPartialList_( positionPartialList )
    {
        std::pair< boost::function< SingleOneWayRangePartialReturnType(
                    const std::vector< basic_mathematics::Vector6d >&, const std::vector< double >& ) >, bool > lightTimeCorrectionPartial;

        for( unsigned int i = 0; i < lighTimeCorrectionPartials.size( ); i++ )
        {
            lightTimeCorrectionPartial = getLightTimeParameterPartialFunction(
                        parameterIdentifier, lighTimeCorrectionPartials.at( i ) );
            if( lightTimeCorrectionPartial.second != 0 )
            {
                lighTimeCorrectionPartialsFunctions_.push_back( lightTimeCorrectionPartial.first );
            }
        }
    }

    ~OneWayRangePartial( ) { }

    virtual OneWayRangePartialReturnType calculatePartial(
            const std::vector< basic_mathematics::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime );

    boost::shared_ptr< OneWayRangeScaling > getOneWayRangeScaler( )
    {
        return oneWayRangeScaler_;
    }

    int getNumberOfLighTimeCorrectionPartialsFunctions( )
    {
        return lighTimeCorrectionPartialsFunctions_.size( );
    }

protected:

    boost::shared_ptr< OneWayRangeScaling > oneWayRangeScaler_;

    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > > positionPartialList_;

    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > >::iterator positionPartialIterator_;

    std::vector< boost::function< SingleOneWayRangePartialReturnType(
            const std::vector< basic_mathematics::Vector6d >&, const std::vector< double >& ) > >
    lighTimeCorrectionPartialsFunctions_;

    std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lighTimeCorrectionPartials_;
};

}

}



#endif // RANGEPARTIAL_H
