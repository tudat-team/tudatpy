/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DETERMINEPARAMETERPOSTFITINFLUENCE_H
#define TUDAT_DETERMINEPARAMETERPOSTFITINFLUENCE_H

#include <algorithm>

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/ObservationModels/simulateObservations.h"
#include "Tudat/SimulationSetup/EstimationSetup/orbitDeterminationManager.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"

namespace tudat
{

namespace simulation_setup
{

//! Function that determines how well the translational dynamics of N bodies can absorb  the influence of a change in an
//! environmental parameter value
/*!
 *  Function that determines how well the translational dynamics of N bodies can absorb  the influence of a change in an
 *  environmental parameter value. The methods outlined by e.g. Dirkx et al. (2016); Planetary and Space Science 134:82-95.
 *  This function uses a nominal dynamical model to simulate ideal observations of a set of N bodies' 3-dimensional positions.
 *  These observations are then used as input to an orbit determination routing in which a set of parameters have their values
 *  adjusted w.r.t. the nominal case. In the orbit determination, only the initial states of the N bodies are estimated. As such,
 *  this function provides the degree to which a change (e.g. uncertainty) in the environment can be mimicked by a change in the
 *  bodies initial conditions.
 *  \param bodyMap List of body objects that comprises the environment
 *  \param integratorSettings Settings for numerical integrator.
 *  \param propagatorSettings Settings for propagator.
 *  \param perturbedParameterSettings Type of parameter that is to be adjusted in analysis.
 *  \param simulatedObservationInterval Time interval between consecutive simulated 3-dimensional position observations
 *  \param parameterPerturbations Perturbations in the parameter vector that are to be used
 *  \param parameterIndices Indices in the vector of perturbed parameter at which to apply the perturbations in
 *  parameterPerturbations, e.g. index parameterIndices( i ) of the parameter vector gets perturbation parameterPerturbations( i )
 *  \param numberOfIterations Number of iterations to use in the orbit determination loop
 *  \return Pair of estimation output (first) and adjustment to initial state vectors (second)
 */
template< typename TimeType = double, typename StateScalarType = double >
std::pair< std::shared_ptr< PodOutput< StateScalarType > >, Eigen::VectorXd > determinePostfitParameterInfluence(
        const NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > perturbedParameterSettings,
        const double simulatedObservationInterval,
        const std::vector< double > parameterPerturbations,
        const std::vector< int > parameterIndices,
        const int numberOfIterations = 2 )
{
    using namespace observation_models;
    using namespace estimatable_parameters;

    // Check input consistency
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< StateScalarType > > translationalPropagatorSettings =
            std::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings ) ;
    if( translationalPropagatorSettings == nullptr )
    {
        throw std::runtime_error(
                    "Error in determinePostfitParameterInfluence, only single-arc translational dynamics currently supported" );
    }

    // Getlist of bodies for which the dynamics is to be fit
    std::vector< std::string > observedBodies = translationalPropagatorSettings->bodiesToIntegrate_;

    // Create list of ideal observation settings and initial states to estimate
    std::vector< LinkEnds > linkEndsList;
    ObservationSettingsMap observationSettingsMap;
    std::vector< std::shared_ptr< EstimatableParameterSettings > > initialStateParameterNames;

    for( unsigned int i = 0; i < observedBodies.size( ); i++ )
    {
        // Add current body to list of observed bodies
        LinkEnds observationLinkEnds;
        observationLinkEnds[ observed_body ] = std::make_pair( observedBodies.at( i ), "" );
        linkEndsList.push_back( observationLinkEnds );
        observationSettingsMap.insert(
                    std::make_pair( observationLinkEnds, std::make_shared< ObservationSettings >(
                                        position_observable ) ) );

        // Add current body to list of estimated bodies
        initialStateParameterNames.push_back(
                    std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                        observedBodies.at( i ), translationalPropagatorSettings->getInitialStates( ).segment( i * 6, 6 ),
                        translationalPropagatorSettings->centralBodies_.at( i ) ) );
    }

    // Create initial state estimation objects
    std::shared_ptr< EstimatableParameterSet< StateScalarType > > initialStateParametersToEstimate =
            createParametersToEstimate< StateScalarType >( initialStateParameterNames, bodyMap );

    // Create orbit determination object.
    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                bodyMap, initialStateParametersToEstimate, observationSettingsMap,
                integratorSettings, propagatorSettings );

    // Retrieve nominal (e.g. pre-fit) body states
    Eigen::VectorXd nominalBodyStates = initialStateParametersToEstimate->template getFullParameterValues< double >( );

    // Get range over which observations are to be simulated.
    std::pair< double, double > dataTimeInterval =
            getTabulatedEphemerisSafeInterval( bodyMap.at( observedBodies.at( 0 ) )->getEphemeris( ) );
    double startTime = dataTimeInterval.first;
    double endTime = dataTimeInterval.second;

    // Define list of times at which observations are to be simulated
    std::vector< TimeType > baseTimeList;
    double currentTime = startTime;
    while( currentTime < endTime )
    {
        baseTimeList.push_back( currentTime );
        currentTime += simulatedObservationInterval;
    }
    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > measurementSimulationInput;
    for( unsigned int i = 0; i < linkEndsList.size( ); i++ )
    {
        measurementSimulationInput[ position_observable ][ linkEndsList.at( i ) ] =
                std::make_pair( baseTimeList, observed_body );
    }

    // Simulate ideal observations
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, LinkEndType > > > SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;
    PodInputDataType observationsAndTimes = simulateObservations< StateScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ) );
    //input_output::writeMatrixToFile( observationsAndTimes.begin( )->second.begin( )->second.first, "preFitObservations.dat" );

    // Define estimation input
    std::shared_ptr< PodInput< StateScalarType, TimeType > > podInput =
            std::make_shared< PodInput< StateScalarType, TimeType > >(
                observationsAndTimes, initialStateParametersToEstimate->getParameterSetSize( ) );
    podInput->defineEstimationSettings( true, true, false, true, true );

    // Create parameters that are to be perturbed
    std::vector< std::shared_ptr< EstimatableParameterSettings > > perturbedParameterSettingsList;
    perturbedParameterSettingsList.push_back( perturbedParameterSettings );
    std::shared_ptr< EstimatableParameterSet< StateScalarType > > perturbedParameters =
            createParametersToEstimate< StateScalarType >( perturbedParameterSettingsList, bodyMap );

    // Perturb parameters by required amount
    Eigen::VectorXd parameterVectorToPerturb = perturbedParameters->template getFullParameterValues< double >( );
    Eigen::VectorXd unperturbedParameterVector = parameterVectorToPerturb;
    for( unsigned int i = 0; i < parameterPerturbations.size( ); i++ )
    {
        if( parameterIndices.at( i ) >= parameterVectorToPerturb.rows( ) )
        {
            throw std::runtime_error( "Error, inconsistent input to determinePostfitParameterInfluence function, perturbation index is too high" );
        }
        else
        {
            parameterVectorToPerturb( parameterIndices.at( i ) ) += parameterPerturbations.at( i );
        }
    }
    perturbedParameters->resetParameterValues( parameterVectorToPerturb );

    // Fit nominal dynamics to pertrubed dynamical model
    std::shared_ptr< PodOutput< StateScalarType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( numberOfIterations ) );

    // Reset parameter to nominal value
    perturbedParameters->resetParameterValues( unperturbedParameterVector );

    return std::make_pair(
                podOutput, initialStateParametersToEstimate->template getFullParameterValues< double >( ) - nominalBodyStates );



}

}

}
#endif // TUDAT_DETERMINEPARAMETERPOSTFITINFLUENCE_H
