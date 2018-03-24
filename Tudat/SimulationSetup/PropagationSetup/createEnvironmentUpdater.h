/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEENVIRONMENTUPDATER_H
#define TUDAT_CREATEENVIRONMENTUPDATER_H

#include <boost/make_shared.hpp>

#include "Tudat/SimulationSetup/PropagationSetup/environmentUpdater.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

namespace tudat
{

namespace propagators
{

//! Function to check validity of existing environment
/*!
 * Function to check whether the requested environment updates are
 * possible with the existing environment, i.e. if all
 * environment models that are to be updated exist in the bodyMap. The
 * function throws an error if environment cannot
 * be updated as requested.,
 * \param requestedUpdates List of environment updates that are required
 * \param bodyMap List of body objects used in the simulations.
 */
void checkValidityOfRequiredEnvironmentUpdates(
        const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >&
        requestedUpdates,
        const simulation_setup::NamedBodyMap& bodyMap );

//! Get list of required environment model update settings from torque models.
/*!
 * Get list of required environment model update settings from torque models.
 * \param torqueModels List of torque models used in simulation.
 * \param bodyMap List of body objects used in the simulations.
 * \return List of required environment model update settings.
 */
std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
createRotationalEquationsOfMotionEnvironmentUpdaterSettings(
        const basic_astrodynamics::TorqueModelMap& torqueModels, const simulation_setup::NamedBodyMap& bodyMap );

//! Get list of required environment model update settings from translational acceleration models.
/*!
 * Get list of required environment model update settings from translational acceleration models.
 * \param translationalAccelerationModels List of acceleration models used in simulation.
 * \param bodyMap List of body objects used in the simulations.
 * \return List of required environment model update settings.
 */
std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
createTranslationalEquationsOfMotionEnvironmentUpdaterSettings(
        const basic_astrodynamics::AccelerationMap& translationalAccelerationModels,
        const simulation_setup::NamedBodyMap& bodyMap );

//! Get list of required environment model update settings from mass rate models.
/*!
 * Get list of required environment model update settings from mass rate models.
 * \param massRateModels List of mass rate models used in simulation.
 * \param bodyMap List of body objects used in the simulations.
 * \return List of required environment model update settings.
 */
std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
createMassPropagationEnvironmentUpdaterSettings(
        const std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::MassRateModel > > > massRateModels,
        const simulation_setup::NamedBodyMap& bodyMap );

//! Get list of required environment model update settings from a list of propagation settings.
/*!
* Get list of required environment model update settings from a list of propagation settings.
* \param propagatorSettings Object providing the full settings for the
* dynamics that are to be propagated.
* \param bodyMap List of body objects used in the simulations.
* \return List of updates required when propagating dynamics defined by propagatorSettings.
*/
template< typename StateScalarType >
std::map< propagators::EnvironmentModelsToUpdate,
    std::vector< std::string > > createEnvironmentUpdaterSettings(
        const boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap )
{
    std::map< propagators::EnvironmentModelsToUpdate,
        std::vector< std::string > > environmentModelsToUpdate;

    // Check dynamics type
    switch( propagatorSettings->getStateType( ) )
    {
    case hybrid:
    {
        // Cast to derived type
        boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

        // Iterate over all propagation settings in hybrid model
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > singleAccelerationUpdateNeeds;

        for( typename std::map< IntegratedStateType,
             std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >::const_iterator
             typeIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             typeIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( ); typeIterator++ )
        {
            for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
            {
                // Create current propagation model from settings list (must not be hybrid).
                if( typeIterator->first != hybrid )
                {
                    singleAccelerationUpdateNeeds = createEnvironmentUpdaterSettings< StateScalarType >(
                                typeIterator->second.at( i ), bodyMap );

                    // Add single model environment model update settings to full list
                    checkValidityOfRequiredEnvironmentUpdates( singleAccelerationUpdateNeeds, bodyMap );
                    addEnvironmentUpdates( environmentModelsToUpdate, singleAccelerationUpdateNeeds );
                }
                else
                {
                    throw std::runtime_error(
                                "Error when making environment updater type list, cannot handle hybrid propagator inside hybrid propagator" );
                }
            }
        }
        break;
    }
    // Retrieve environment model settings for translational dynamics
    case transational_state:
    {
        environmentModelsToUpdate = createTranslationalEquationsOfMotionEnvironmentUpdaterSettings(
                    boost::dynamic_pointer_cast<
                    TranslationalStatePropagatorSettings< StateScalarType > >(
                        propagatorSettings )->getAccelerationsMap( ),
                    bodyMap );
        break;
    }
    case rotational_state:
    {
        environmentModelsToUpdate = createRotationalEquationsOfMotionEnvironmentUpdaterSettings(
                    boost::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >( propagatorSettings )->getTorqueModelsMap( ),
                    bodyMap );
        break;
    }
    // Retrieve environment model settings for mass rate model
    case body_mass_state:
    {
        environmentModelsToUpdate = createMassPropagationEnvironmentUpdaterSettings(
                    boost::dynamic_pointer_cast<
                    MassPropagatorSettings< StateScalarType > >(
                        propagatorSettings )->getMassRateModelsMap( ), bodyMap );
        break;
    }
    case custom_state:
    {
        break;
    }
    default:
    {
        throw std::runtime_error( "Error, cannot create environment updates for type " +
                                  std::to_string( propagatorSettings->getStateType( ) ) );
    }
    }
    return environmentModelsToUpdate;

}

//! Function to create 'brute-force' update settings, in which each environment model is updated.
/*!
 * Function to create 'brute-force' update settings, in which each
 * environment model (i.e. each member of the Body
 * objects that can be updated) is updated.
 * \param bodyMap List of body objects used in the simulations.
 * \return List of environment model updates, so that each updatable model is updated.
 */
std::map< propagators::EnvironmentModelsToUpdate,
    std::vector< std::string > > createFullEnvironmentUpdaterSettings(
        const simulation_setup::NamedBodyMap& bodyMap );

//! Create environment updater from a list of propagation settings.
/*!
* Get environment updater from a list of propagation settings.
* \param propagatorSettings Object providing the full settings for the
* dynamics that are to be propagated.
* \param bodyMap List of body objects used in the simulations.
* \return Object that updates the environment in bodyMap, as per thhe
* requirements set by propagatorSettings.
*/
template< typename StateScalarType, typename TimeType >
boost::shared_ptr< propagators::EnvironmentUpdater< StateScalarType, TimeType > >
createEnvironmentUpdaterForDynamicalEquations(
        const boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap )
{
    // Create environment update settings.
    std::map< IntegratedStateType,
        std::vector< std::pair< std::string, std::string > > >integratedTypeAndBodyList =
            getIntegratedTypeAndBodyList< StateScalarType >( propagatorSettings );

    std::map< propagators::EnvironmentModelsToUpdate,
        std::vector< std::string > > environmentModelsToUpdate =
            createEnvironmentUpdaterSettings< StateScalarType >( propagatorSettings, bodyMap );

    // Create and return environment updater object.
    return boost::make_shared< EnvironmentUpdater< StateScalarType, TimeType > >(
                bodyMap, environmentModelsToUpdate, integratedTypeAndBodyList );
}

} // namespace propagators

} // namespace tudat

#endif // TUDAT_CREATEENVIRONMENTUPDATER_H
