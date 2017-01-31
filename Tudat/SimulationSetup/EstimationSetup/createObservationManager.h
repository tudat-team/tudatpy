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

template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > createObservationSimulator(
        const ObservableType observableType,
        const std::vector< LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap &bodyMap,
        const LightTimeCorrectionSettingsMap& singleObservableCorrections =
        LightTimeCorrectionSettingsMap( ) )
{
    std::map< LinkEnds, boost::shared_ptr< ObservationModel< ObservationSize, ObservationScalarType, TimeType > > > observationModels;

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

        observationModels[ linkEnds[ i ] ] = ObservationModelCreator< ObservationSize, ObservationScalarType, TimeType >::createObservationModel(
                    observableType, linkEnds.at( i ), bodyMap, currentLightTimeCorrections );

    }

    return boost::make_shared< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > >(
                observableType, observationModels );


}

template< int ObservationSize = 1, typename ObservationScalarType = double, typename TimeType = double >
boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > createObservationManager(
        const ObservableType observableType,
        const std::vector< LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap &bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
        const boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const LightTimeCorrectionSettingsMap& singleObservableCorrections =
        LightTimeCorrectionSettingsMap( ) )
{
    using namespace observation_models;
    using namespace observation_partials;

    boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > observationManager;

    boost::shared_ptr< ObservationSimulator< ObservationSize, ObservationScalarType, TimeType > > observationSimulator =
            createObservationSimulator< ObservationSize, ObservationScalarType, TimeType >(
                observableType, linkEnds, bodyMap, singleObservableCorrections );


    PerLinkEndPerLightTimeSolutionCorrections lightTimeCorrectionList =
            getLightTimeCorrectionsList( observationSimulator->getObservationModels( ) );
    boost::shared_ptr< ObservationPartialCreator< ObservationSize, ObservationScalarType > > observationPartialCreator;
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


    observationManager = boost::make_shared< ObservationManager< ObservationSize, ObservationScalarType, TimeType > >(
                observableType, observationSimulator, observationPartials, observationPartialScalers, stateTransitionMatrixInterface );

    return observationManager;
}

template< typename ObservationScalarType = double, typename TimeType = double >
boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > createObservationManagerBase(
        const ObservableType observableType,
        const std::vector< LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap &bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
        const boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const LightTimeCorrectionSettingsMap& singleObservableCorrections =
        LightTimeCorrectionSettingsMap( ) )
{
    boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > observationManager;
    switch( observableType )
    {
    case oneWayRange:
        observationManager = createObservationManager< 1, ObservationScalarType, TimeType >(
                    observableType, linkEnds, bodyMap, parametersToEstimate, stateTransitionMatrixInterface, singleObservableCorrections );
        break;
    case angular_position:
        observationManager = createObservationManager< 2, ObservationScalarType, TimeType >(
                    observableType, linkEnds, bodyMap, parametersToEstimate, stateTransitionMatrixInterface, singleObservableCorrections );
        break;
    case position_observable:
        observationManager = createObservationManager< 3, ObservationScalarType, TimeType >(
                    observableType, linkEnds, bodyMap, parametersToEstimate, stateTransitionMatrixInterface, singleObservableCorrections );
        break;
    default:
        std::cerr<<"Error when making observation manager, could not identify observable type "<<observableType<<std::endl;
    }
    return observationManager;
}

template< typename ObservationScalarType = double, typename TimeType = double >
boost::shared_ptr< ObservationManagerBase< ObservationScalarType, TimeType > > createObservationManagerBaseWithInterfaceCreator(
        const ObservableType observableType,
        const std::vector< LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap &bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
        const boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const LightTimeCorrectionSettingsMap& singleObservableCorrections =
        LightTimeCorrectionSettingsMap( ) )
{
    return createObservationManagerBase< ObservationScalarType, TimeType >(
                observableType, linkEnds, bodyMap, parametersToEstimate, stateTransitionMatrixInterface, singleObservableCorrections );
}


}


}


#endif // CREATEOBSERVATIONMANAGER_H
