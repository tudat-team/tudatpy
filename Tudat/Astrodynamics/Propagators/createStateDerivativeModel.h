/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATESTATEDERIVATIVEMODEL_H
#define TUDAT_CREATESTATEDERIVATIVEMODEL_H

#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"
#include "Tudat/Astrodynamics/Propagators/nBodyCowellStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyEnckeStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/bodyMassStateDerivative.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace propagators
{

//! Function to create object handling frame origin transformations during numerical integration
/*!
 * Function to create object handling frame origin transformations during numerical integration
 * (needed by NBodyStateDerivative).
 * \param centralBodies Names of central bodies, belonging to the
 * entries in the bodiesToIntegrate vector of same index.
 * \param bodiesToIntegrate Names of bodies that are to be integrated numerically.
 * \param bodyMap List of body objects used in simulation.
 * \return Object handling frame origin transformations during numerical integration
 */
template< typename StateScalarType, typename TimeType >
boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > createCentralBodyData(
        const std::vector< std::string >& centralBodies,
        const std::vector< std::string >& bodiesToIntegrate,
        const simulation_setup::NamedBodyMap& bodyMap )
{

    // Check whether the bodies that are to e integrated exist in bodyMap
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        if( bodyMap.count( bodiesToIntegrate.at( i ) ) == 0 )
        {
            throw std::runtime_error(
                        "Warning when creating CentralBodyData, body " + bodiesToIntegrate.at( i )
                        + " not present in provided body map." );
        }
    }

    std::vector< std::string > centralBodiesToUse;
    // Generate central body data for requested settings; if no
    // central bodies provided, get inertial SSB as central body
    if( centralBodies.size( ) == 0 )
    {
        for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
        {
            centralBodiesToUse.push_back( "SSB" );
        }
    }
    else
    {
        centralBodiesToUse = centralBodies;
    }

    std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >
        ( const TimeType ) > > bodyStateFunctions;

    // Retrieve frame origin state functions
    for( unsigned int i = 0; i < centralBodiesToUse.size( ); i++ )
    {
        if( centralBodiesToUse.at( i )  != "SSB" )
        {
            bodyStateFunctions[ centralBodiesToUse.at( i ) ] =
                    boost::bind( &simulation_setup::Body::getTemplatedStateInBaseFrameFromEphemeris
                                 < StateScalarType, TimeType >,
                                 bodyMap.at( centralBodiesToUse.at( i ) ), _1 );
        }
        else
        {
            bodyStateFunctions[ centralBodiesToUse.at( i ) ] = boost::lambda::constant(
                        Eigen::Matrix< StateScalarType, 6, 1 >::Zero( ) );
        }
    }
    return boost::make_shared< CentralBodyData< StateScalarType, TimeType > >(
                centralBodiesToUse, bodiesToIntegrate, bodyStateFunctions );
}

//! Function to create a translational state derivative model.
/*!
 *  Function to create a translational state derivative model from
 *  propagation settings and environment.
 *  \param translationPropagatorSettings Settings for the translational dynamics model.
 *  \param bodyMap List of body objects in the environment
 *  \param propagationStartTime Time from which numerical propagation starts.
 *  \return Translational state derivative model (instance of derived class of NBodyStateDerivative)
 */
template< typename StateScalarType = double, typename TimeType = double >
boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > >
createTranslationalStateDerivativeModel(
        const boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
        translationPropagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const TimeType propagationStartTime )
{

    // Create object for frame origin transformations.
    boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData =
            createCentralBodyData< StateScalarType, TimeType >(
                translationPropagatorSettings->centralBodies_,
                translationPropagatorSettings->bodiesToIntegrate_,
                bodyMap );

    boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel;

    // Check propagator type and create corresponding state derivative object.
    switch( translationPropagatorSettings->propagator_ )
    {
    case cowell:
    {
        stateDerivativeModel = boost::make_shared<
                               NBodyCowellStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->accelerationsMap_, centralBodyData,
                  translationPropagatorSettings->bodiesToIntegrate_ );
        break;
    }
    case encke:
    {
        // Calculate initial Kepler elements for Encke propagator
        std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > initialKeplerElements;
        initialKeplerElements.resize( translationPropagatorSettings->bodiesToIntegrate_.size( ) );
        std::vector< std::string > centralBodies = translationPropagatorSettings->centralBodies_;

        for( unsigned int i = 0; i < translationPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
        {
            if( bodyMap.count( centralBodies[ i ] ) == 0 )
            {
                std::cerr<<"Error when creating Encke propagator, did not find central body "<<centralBodies[ i ]<<std::endl;
            }
            initialKeplerElements[ i ] = orbital_element_conversions::convertCartesianToKeplerianElements< StateScalarType >(
                        translationPropagatorSettings->getInitialStates( ).segment( i * 6, 6 ), static_cast< StateScalarType >(
                            bodyMap.at( centralBodies[ i ] )->getGravityFieldModel( )->getGravitationalParameter( ) ) );
        }

        // Create Encke state derivative object.
        stateDerivativeModel = boost::make_shared< NBodyEnckeStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->accelerationsMap_, centralBodyData, translationPropagatorSettings->bodiesToIntegrate_,
                  initialKeplerElements, propagationStartTime );

        break;
    }
    default:
        throw std::runtime_error(
            "Error, did not recognize translational state propagation type: " +
            boost::lexical_cast< std::string >( translationPropagatorSettings->propagator_ ) );
    }
    return stateDerivativeModel;
}

//! Function to create a mass state derivative model.
/*!
 *  Function to create a mass state derivative model from propagation settings and environment.
 *  \param massPropagatorSettings Settings for the mass dynamics model.
 *  \param bodyMap List of body objects in the environment
 *  \return Mass state derivative model.
 */
template< typename StateScalarType = double, typename TimeType = double >
boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > createBodyMassStateDerivativeModel(
        const boost::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings,
        const  simulation_setup::NamedBodyMap& bodyMap )
{
    return boost::make_shared< propagators::BodyMassStateDerivative< StateScalarType, TimeType > >(
                massPropagatorSettings->massRateModels_,
                massPropagatorSettings->bodiesWithMassToPropagate_ );
}

//! Function to create a state derivative model.
/*!
 *  Function to create a state derivative model from propagation settings and the environment.
 *  \param propagatorSettings Settings for the dynamical model.
 *  \param bodyMap List of body objects in the environment
 *  \param propagationStartTime Time from which numerical propagation starts.
 *  \return State derivative model (instance of required derived class of SingleStateTypeDerivative)
 */
template< typename StateScalarType = double, typename TimeType = double >
boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > >
createStateDerivativeModel(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const TimeType propagationStartTime )
{
    boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel;

    // Check dynamics type and call associated function to create
    // specific type of state derivative model.
    switch( propagatorSettings->stateType_ )
    {
    case transational_state:
    {
        // Check input consistency.
        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
                translationPropagatorSettings =
                boost::dynamic_pointer_cast<
                TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        if( translationPropagatorSettings == NULL )
        {
            throw std::runtime_error(
                "Error, expected translational state propagation settings when making state derivative model" );
        }
        else
        {
            stateDerivativeModel = createTranslationalStateDerivativeModel< StateScalarType, TimeType >(
                        translationPropagatorSettings, bodyMap, propagationStartTime );
        }
        break;
    }
    case body_mass_state:
    {
        // Check input consistency.
        boost::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings =
                boost::dynamic_pointer_cast< MassPropagatorSettings< StateScalarType > >( propagatorSettings );
        if( massPropagatorSettings == NULL )
        {
            throw std::runtime_error(
                "Error, expected mass propagation settings when making state derivative model" );
        }
        else
        {
            stateDerivativeModel = createBodyMassStateDerivativeModel< StateScalarType, TimeType >(
                        massPropagatorSettings, bodyMap );
        }
        break;
    }
    default:
        throw std::runtime_error(
                    "Error, could not process state type "
                    + boost::lexical_cast< std::string >( propagatorSettings->stateType_ )
                    + " when making state derivative model" );
    }
    return stateDerivativeModel;
}

//! Function to create a list of state derivative models.
/*!
 *  Function to create a list of state derivative models from
 *  propagation settings and the environment.
 *  \param propagatorSettings Settings for the dynamical model.
 *  \param bodyMap List of body objects in the environment.
 *  \param propagationStartTime Time from which numerical propagation starts.
 *  \return List of state derivative models (instances of required
 *  derived class of SingleStateTypeDerivative)
 */
template< typename StateScalarType = double, typename TimeType = double >
std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > >
createStateDerivativeModels(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const TimeType propagationStartTime )
{
    std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > >
    stateDerivativeModels;

    // Check type of state derivative model and call associated create function.
    switch( propagatorSettings->stateType_ )
    {
    // If hybrid, call create function separately for each entry.
    case hybrid:
    {
        boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

        // Iterate over all propagation settings
        for( typename std::map< IntegratedStateType,
             std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::iterator
             propagatorIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             propagatorIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( ); propagatorIterator++ )
        {
            for( unsigned int i = 0; i < propagatorIterator->second.size( ); i++ )
            {
                // Call create function for current propagation settings.
                if( propagatorIterator->first != hybrid )
                {
                    stateDerivativeModels.push_back( createStateDerivativeModel< StateScalarType, TimeType >(
                                                         propagatorIterator->second.at( i ), bodyMap, propagationStartTime ) );
                }
                else
                {
                    throw std::runtime_error(
                                "Error when making state derivative model, cannot process nested hybrid propagators" );
                }
            }
        }
        break;
    }
    // If not hybrid, call create function for single object directly.
    default:
        stateDerivativeModels.push_back( createStateDerivativeModel< StateScalarType, TimeType >(
                                             propagatorSettings, bodyMap, propagationStartTime ) );
    }

    return stateDerivativeModels;
}

} // namespace propagators
} // namespace tudat

#endif // TUDAT_CREATESTATEDERIVATIVEMODEL_H
