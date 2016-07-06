
#ifndef CREATELIGHTTIMECALCULATOR_H
#define CREATELIGHTTIMECALCULATOR_H

#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/createLightTimeCorrection.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{
namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > > createLightTimeCalculator(
        const boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) >& transmitterCompleteEphemeris,
        const boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) >& receiverCompleteEphemeris,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections,
        const LinkEndId& transmittingLinkEnd,
        const LinkEndId& receivingLinkEnd )
{
    std::vector< LightTimeCorrectionFunction > lightTimeCorrection;

    for( unsigned int i = 0; i < lightTimeCorrections.size( ); i++ )
    {
        LightTimeCorrectionFunction correctionFunction =
                boost::bind( &LightTimeCorrection::calculateLightTimeCorrection, createLightTimeCorrections(
                                        lightTimeCorrections[ i ], bodyMap, transmittingLinkEnd, receivingLinkEnd ), _1, _2, _3, _4 );
        lightTimeCorrection.push_back(  correctionFunction );
    }

    return boost::make_shared< LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > >
            ( transmitterCompleteEphemeris, receiverCompleteEphemeris, lightTimeCorrection );
}

template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
boost::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType, StateScalarType > > createLightTimeCalculator(
        const LinkEndId& transmittingLinkEnd,
        const LinkEndId& receivingLinkEnd,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections )
{
   if( transmittingLinkEnd.second != "" || receivingLinkEnd.second != "" )
   {
       throw std::runtime_error( "Error, body reference points not yet supported when creating light time calculator" );
   }

   return createLightTimeCalculator< ObservationScalarType, TimeType, StateScalarType >(
               boost::bind( &simulation_setup::Body::getTemplatedStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
                            bodyMap.at( transmittingLinkEnd.first ), _1 ),
               boost::bind( &simulation_setup::Body::getTemplatedStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
                            bodyMap.at( receivingLinkEnd.first ), _1 ),
               bodyMap, lightTimeCorrections, transmittingLinkEnd, receivingLinkEnd );
}

}

}
#endif // CREATELIGHTTIMECALCULATOR_H
