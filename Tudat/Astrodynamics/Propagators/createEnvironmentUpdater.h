#ifndef CREATEENVIRONMENTUPDATER_H
#define CREATEENVIRONMENTUPDATER_H

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/Propagators/environmentUpdater.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

namespace tudat
{

namespace propagators
{


void checkValidityOfRequiredEnvironmentUpdates(
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >&
        requestedUpdates,
        const simulation_setup::NamedBodyMap& bodyMap );

void addEnvironmentUpdates(
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >&
        environmentUpdateList,
        const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
        updatesToAdd );

std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
createTranslationalEquationsOfMotionEnvironmentUpdaterSettings(
        const basic_astrodynamics::AccelerationMap& translationalAccelerationModels,
        const simulation_setup::NamedBodyMap& bodyMap );

template< typename StateScalarType >
std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createEnvironmentUpdaterSettings(
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings )
{
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate;

    switch( propagatorSettings->stateType_ )
    {
    case transational_state:
    {
        environmentModelsToUpdate = createTranslationalEquationsOfMotionEnvironmentUpdaterSettings(
                    boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings )->accelerationsMap_,
                    bodyMap );
        break;
    }
    default:
    {
        throw std::runtime_error( "Error, cannot create environment updates for type " +
                                  boost::lexical_cast< std::string >( propagatorSettings->stateType_ ) );
    }
    }
    return environmentModelsToUpdate;

}

std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createFullEnvironmentUpdaterSettings(
        const simulation_setup::NamedBodyMap& bodyMap );


template< typename StateScalarType, typename TimeType >
boost::shared_ptr< propagators::EnvironmentUpdater< StateScalarType, TimeType > > createEnvironmentUpdaterForDynamicalEquations(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap )
{    
    std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > >integratedTypeAndBodyList =
            getIntegratedTypeAndBodyList< StateScalarType >( propagatorSettings );
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
            createEnvironmentUpdaterSettings< StateScalarType >( bodyMap, propagatorSettings );

    return boost::make_shared< EnvironmentUpdater< StateScalarType, TimeType > >(
                bodyMap, environmentModelsToUpdate, integratedTypeAndBodyList );
}

}

}

#endif // CREATEENVIRONMENTUPDATER_H
