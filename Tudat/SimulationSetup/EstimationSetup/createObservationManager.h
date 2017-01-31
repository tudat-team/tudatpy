#ifndef CREATEOBSERVATIONMANAGER_H
#define CREATEOBSERVATIONMANAGER_H

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/ObservationModels/observationManager.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationPartials.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"

namespace tudat
{

namespace observation_models
{

template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > createObservationSimulator(
        const ObservableType observableType,
        const std::vector< LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap &bodyMap,
        const LightTimeCorrectionSettingsMap& singleObservableCorrections =
        LightTimeCorrectionSettingsMap( ) )
{
    std::map< LinkEnds, boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > > observationModels;

    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
//        boost::shared_ptr< ObservationBiasInterface > currentBiasInterface = boost::make_shared< ObservationBiasInterface >( ObservationSize );
//        if( observationBiasInterfaceList.count( observableType ) > 0 )
//        {
//            if( observationBiasInterfaceList.at( observableType ).count( linkEnds[ i ] ) > 0 )
//            {
//                currentBiasInterface = observationBiasInterfaceList.at( observableType ).at( linkEnds[ i ] );
//            }
//        }

        std::vector< boost::shared_ptr< LightTimeCorrectionSettings > > currentLightTimeCorrections;
        if( singleObservableCorrections.count( linkEnds.at( i ) ) > 0 )
        {
            currentLightTimeCorrections = singleObservableCorrections.at( linkEnds.at( i ) );
        }

        observationModels[ linkEnds[ i ] ] = ObservationModelCreator< ObservationSize, ObservationScalarType, TimeType, StateScalarType >::createObservationModel(
                    observableType, linkEnds.at( i ), bodyMap, currentLightTimeCorrections );

    }

    return boost::make_shared< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType, StateScalarType > >(
                observableType, observationModels );


}

template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType > > createObservationManager(
        const ObservableType observableType,
        const std::vector< LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap &bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const LightTimeCorrectionSettingsMap& singleObservableCorrections =
        LightTimeCorrectionSettingsMap( ) )
{
    using namespace observation_models;
    using namespace observation_partials;

    boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType > > observationManager;

    boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType, StateScalarType > > observationSimulator =
            createObservationSimulator< ObservationSize, ObservationScalarType, TimeType, StateScalarType >(
                observableType, linkEnds, bodyMap, singleObservableCorrections );


    PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrectionList =
            getLightTimeCorrectionsList( observationSimulator->getObservationModels( ) );
    boost::shared_ptr< ObservationPartialCreator< ObservationSize, StateScalarType > > observationPartialCreator;
    std::map< LinkEnds, std::pair< std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< ObservationSize > > >,
            boost::shared_ptr< PositionPartialScaling > > > observationPartialsAndScaler;

    if( parametersToEstimate != NULL )
    {
        observationPartialsAndScaler =
            observationPartialCreator->createObservationPartials(
                observableType, linkEnds, bodyMap, parametersToEstimate,
                lightTimeCorrectionList );
    }

    std::map< LinkEnds, std::map< std::pair< int, int >, boost::shared_ptr< observation_partials::ObservationPartial< ObservationSize > > > > observationPartials;
    std::map< LinkEnds, boost::shared_ptr< observation_partials::PositionPartialScaling  > > observationPartialScalers;
    splitObservationPartialsAndScalers( observationPartialsAndScaler, observationPartials, observationPartialScalers );


    observationManager = boost::make_shared< ObservationManager< ObservationSize, ObservationScalarType, TimeType, StateScalarType > >(
                observableType, observationSimulator, observationPartials, observationPartialScalers, stateTransitionMatrixInterface );

    return observationManager;
}

template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType > > createObservationManagerBase(
        const ObservableType observableType,
        const std::vector< LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap &bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const LightTimeCorrectionSettingsMap& singleObservableCorrections =
        LightTimeCorrectionSettingsMap( ) )
{
    boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType > > observationManager;
    switch( observableType )
    {
    case oneWayRange:
        observationManager = createObservationManager< 1, ObservationScalarType, TimeType, StateScalarType >(
                    observableType, linkEnds, bodyMap, parametersToEstimate, stateTransitionMatrixInterface, singleObservableCorrections );
        break;
    case angular_position:
        observationManager = createObservationManager< 2, ObservationScalarType, TimeType, StateScalarType >(
                    observableType, linkEnds, bodyMap, parametersToEstimate, stateTransitionMatrixInterface, singleObservableCorrections );
        break;
    case position_observable:
        observationManager = createObservationManager< 3, ObservationScalarType, TimeType, StateScalarType >(
                    observableType, linkEnds, bodyMap, parametersToEstimate, stateTransitionMatrixInterface, singleObservableCorrections );
        break;
    default:
        std::cerr<<"Error when making observation manager, could not identify observable type "<<observableType<<std::endl;
    }
    return observationManager;
}

template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType > > createObservationManagerBaseWithInterfaceCreator(
        const ObservableType observableType,
        const std::vector< LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap &bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const LightTimeCorrectionSettingsMap& singleObservableCorrections =
        LightTimeCorrectionSettingsMap( ) )
{
    return createObservationManagerBase< ObservationScalarType, TimeType, StateScalarType, StateScalarType >(
                observableType, linkEnds, bodyMap, parametersToEstimate, stateTransitionMatrixInterface, singleObservableCorrections );
}


}


}


#endif // CREATEOBSERVATIONMANAGER_H
