/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEESTIMATABLEPARAMETERS_H
#define TUDAT_CREATEESTIMATABLEPARAMETERS_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/SimulationSetup/EstimationSetup/estimatableParameterSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create interface object for estimating parameters representing an initial dynamical state.
/*!
 *  Function to create interface object for estimating parameters representing an initial dynamical state.
 *  \param bodyMap Map of body objects containing the fll simulation environment.
 *  \param parameterSettings Object defining the parameter interface that is to be created.
 *  \return Interface object for estimating an initial state.
 */
template< typename InitialStateParameterType = double >
boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix
< InitialStateParameterType, Eigen::Dynamic, 1 > > > createInitialDynamicalStateParameterToEstimate(
        const NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSettings >& parameterSettings )
{
    using namespace tudat::estimatable_parameters;

    boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > >
            initialStateParameterToEstimate;

    // Check consistency of input.
    if( !isParameterDynamicalPropertyInitialState( parameterSettings->parameterType_.first ) )
    {
        std::string errorMessage = "Error when requesting to make initial state parameter " +
                std::to_string( parameterSettings->parameterType_.first ) + " of " +
                parameterSettings->parameterType_.second.first +
                ", parameter is not an initial state parameter ";
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Identify state that is to be estimation
        switch( parameterSettings->parameterType_.first )
        {
        case initial_body_state:

            // Check consistency of input.
            if( boost::dynamic_pointer_cast<
                    InitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                        parameterSettings ) == NULL )
            {
                throw std::runtime_error( "Error when making body initial state parameter, settings type is incompatible" );
            }
            else
            {
                boost::shared_ptr< InitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >
                        initialStateSettings =
                        boost::dynamic_pointer_cast<
                        InitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                            parameterSettings );

                Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialTranslationalState;

                // If initial time is not defined, use preset initial state
                if( ! ( initialStateSettings->initialTime_ == initialStateSettings->initialTime_  ) )
                {
                    initialTranslationalState = initialStateSettings->initialStateValue_;


                }
                // Compute initial state from environment
                else
                {
                    initialTranslationalState = propagators::getInitialStateOfBody
                            < double, InitialStateParameterType >(
                                initialStateSettings->parameterType_.second.first, initialStateSettings->centralBody_,
                                bodyMap, initialStateSettings->initialTime_ );

                }

                // Create translational state estimation interface object
                initialStateParameterToEstimate =
                        boost::make_shared< InitialTranslationalStateParameter< InitialStateParameterType > >(
                            initialStateSettings->parameterType_.second.first, initialTranslationalState,
                            initialStateSettings->centralBody_,
                            initialStateSettings->frameOrientation_ );
            }
            break;
        case arc_wise_initial_body_state:
            if( boost::dynamic_pointer_cast< ArcWiseInitialTranslationalStateEstimatableParameterSettings<
                    InitialStateParameterType > >( parameterSettings ) == NULL )
            {
                throw std::runtime_error(
                            "Error when making body initial state parameter, settings type is incompatible" );
            }
            else
            {
                boost::shared_ptr< ArcWiseInitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >
                        initialStateSettings =  boost::dynamic_pointer_cast<
                        ArcWiseInitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                            parameterSettings );

                if( initialStateSettings->isStateSet_ )
                {
                    initialStateParameterToEstimate = boost::make_shared< ArcWiseInitialTranslationalStateParameter<
                            InitialStateParameterType > >(
                                initialStateSettings->parameterType_.second.first,
                                initialStateSettings->arcStartTimes_,
                                initialStateSettings->initialStateValue_,
                                initialStateSettings->centralBody_,
                                initialStateSettings->frameOrientation_ );
                }
                else
                {
                    initialStateParameterToEstimate = boost::make_shared< ArcWiseInitialTranslationalStateParameter<
                            InitialStateParameterType > >(
                                initialStateSettings->parameterType_.second.first, initialStateSettings->arcStartTimes_,
                                propagators::getInitialArcWiseStateOfBody< double, InitialStateParameterType >(
                                    initialStateSettings->parameterType_.second.first,
                                    initialStateSettings->centralBody_, bodyMap,
                                    initialStateSettings->arcStartTimes_ ),
                                initialStateSettings->centralBody_, initialStateSettings->frameOrientation_ );
                }
            }
            break;
        default:
            std::string errorMessage = "Error, could not create parameter for initial state of type " +
                    std::to_string( parameterSettings->parameterType_.first );
            throw std::runtime_error( errorMessage );
        }
    }

    return initialStateParameterToEstimate;
}

//! Function to create an interface object for estimating a parameter defined by a single double value
/*!
 * Function to create an interface object for estimating a parameter defined by a single double value
 * \param doubleParameterName Object defining the parameter interface that is to be created.
 * \param bodyMap Map of body objects containing the fll simulation environment.
 * \param accelerationModelMap List of acceleration models used in simulations; empty by default (only required for
 * selected parameters).
 * \return Interface object for estimating parameter.
 */
boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > createDoubleParameterToEstimate(
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSettings >& doubleParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap =
        basic_astrodynamics::AccelerationMap( ) );

//! Function to create an interface object for estimating a parameter defined by a list of double values
/*!
 * Function to create an interface object for estimating a parameter defined by a list of single double values
 * \param vectorParameterName Object defining the parameter interface that is to be created.
 * \param bodyMap Map of body objects containing the fll simulation environment.
 * \param accelerationModelMap List of acceleration models used in simulations; empty by default (only required for
 * selected parameters).
 * \return Interface object for estimating parameter.
 */
boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > createVectorParameterToEstimate(
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSettings >& vectorParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap =
        basic_astrodynamics::AccelerationMap( ) );

//! Function to create the interface object for estimating any number/type of parameters.
/*!
 *  Function to create the interface object for estimating any number/type of parameters. This can include both
 *  environmental parameters and initial dynamical states. The types of parameters are defined by the parameterNames m
 *  input variables
 *  \param parameterNames List of objects defining the parameters that are to be estimated.
 *  \param bodyMap Map of body objects containing the fll simulation environment.
 *  \param accelerationModelMap List of acceleration models used in simulations; empty by default (only required for
 *  selected parameters).
 *  \return Interface object for estimating a set of parameters.
 */
template< typename InitialStateParameterType = double >
boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > createParametersToEstimate(
        const std::vector< boost::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >& parameterNames,
        const NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap = basic_astrodynamics::AccelerationMap( ) )

{
    using namespace tudat::estimatable_parameters;

    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
            initialDynamicalParametersToEstimate;
    std::vector< boost::shared_ptr< EstimatableParameter< double > > > doubleParametersToEstimate;
    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParametersToEstimate;

    // Iterate over all parameters.
    for( unsigned int i = 0; i < parameterNames.size( ); i++ )
    {
        // Create initial dynamical parameters.
        if( isParameterDynamicalPropertyInitialState( parameterNames.at( i )->parameterType_.first ) )
        {
            initialDynamicalParametersToEstimate.push_back(
                        createInitialDynamicalStateParameterToEstimate< InitialStateParameterType >(
                            bodyMap, parameterNames.at( i ) ) );
        }
        // Create parameters defined by single double value
        else if( isDoubleParameter( parameterNames[ i ]->parameterType_.first ) == true )
        {
            doubleParametersToEstimate.push_back( createDoubleParameterToEstimate(
                                                      parameterNames[ i ], bodyMap, accelerationModelMap ) );
        }
        // Create parameters defined by list of double values
        else if( isDoubleParameter( parameterNames[ i ]->parameterType_.first ) == false )
        {
            vectorParametersToEstimate.push_back( createVectorParameterToEstimate(
                                                      parameterNames[ i ], bodyMap, accelerationModelMap ) );
        }
        else
        {
            std::string errorMessage = "Error, parameter type of  " +
                    std::string( parameterNames[ i ]->parameterType_.second.first ) + "of " +
                    std::to_string( parameterNames[ i ]->parameterType_.first ) +
                    "not recognized when making estimatable parameter set.";

            throw std::runtime_error( errorMessage );
        }
    }

    return boost::make_shared< EstimatableParameterSet< InitialStateParameterType > >(
                doubleParametersToEstimate, vectorParametersToEstimate, initialDynamicalParametersToEstimate );
}


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEESTIMATABLEPARAMETERS_H
