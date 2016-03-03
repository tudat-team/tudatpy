#ifndef CREATEENVIRONMENTUPDATER_H
#define CREATEENVIRONMENTUPDATER_H

#include "Tudat/Astrodynamics/Propagators/environmentUpdater.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

namespace tudat
{

namespace propagators
{


void checkValidityOfRequiredEnvironmentUpdates(
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& requestedUpdates,
        const NamedBodyMap& bodyMap );

void addEnvironmentUpdates(
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& environmentUpdateList,
        const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > updatesToAdd );

std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createTranslationalEquationsOfMotionEnvironmentUpdaterSettings(
        const AccelerationMap& translationalAccelerationModels, const NamedBodyMap& bodyMap );

std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createRotationalEquationsOfMotionEnvironmentUpdaterSettings(
        const basic_astrodynamics::TorqueModelMap& torqueModels, const NamedBodyMap& bodyMap );

std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createProperTimeEquationEnvironmentUpdaterSettings(
        const boost::shared_ptr< RelativisticTimeStatePropagatorSettings > stateDerivativeModel, const NamedBodyMap& bodyMap );

template< typename StateScalarType >
std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createEnvironmentUpdaterSettings(
        const NamedBodyMap& bodyMap,
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings )
{
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate;

    switch( propagatorSettings->stateType_ )
    {
    case hybrid:
    {
        boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > singleAccelerationUpdateNeeds;

        for( typename std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
             typeIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             typeIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( ); typeIterator++ )
        {
            for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
            {
                if( typeIterator->first != hybrid )
                {
                    singleAccelerationUpdateNeeds = createEnvironmentUpdaterSettings< StateScalarType >( bodyMap, typeIterator->second.at( i ) );

                    checkValidityOfRequiredEnvironmentUpdates( singleAccelerationUpdateNeeds, bodyMap );
                    addEnvironmentUpdates( environmentModelsToUpdate, singleAccelerationUpdateNeeds );
                }
                else
                {
                    std::cerr<<"Error when making environment updater type list, cannot handle hybrid propagator inside hybrid propagator"<<std::endl;
                }
            }
        }

        break;
    }
    case transational_state:
    {
        environmentModelsToUpdate = createTranslationalEquationsOfMotionEnvironmentUpdaterSettings(
                    boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings )->accelerationsMap_,
                    bodyMap );
        break;
    }
    case rotational_state:
    {
        environmentModelsToUpdate = createRotationalEquationsOfMotionEnvironmentUpdaterSettings(
                    boost::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >( propagatorSettings )->torqueModelMap_,
                    bodyMap );
        break;
    }
    case proper_time:
    {
        environmentModelsToUpdate = createProperTimeEquationEnvironmentUpdaterSettings(
                    boost::dynamic_pointer_cast< RelativisticTimeStatePropagatorSettings >( propagatorSettings ), bodyMap );
        break;
    }
    }
    return environmentModelsToUpdate;

}

std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createFullEnvironmentUpdaterSettings(
        const NamedBodyMap& bodyMap );


template< typename StateScalarType, typename TimeType >
boost::shared_ptr< propagators::EnvironmentUpdater< StateScalarType, TimeType > > createEnvironmentUpdaterForDynamicalEquations(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const NamedBodyMap& bodyMap )
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
