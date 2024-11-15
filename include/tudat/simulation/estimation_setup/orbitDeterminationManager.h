/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ORBITDETERMINATIONMANAGER_H
#define TUDAT_ORBITDETERMINATIONMANAGER_H

#include <algorithm>



#include "tudat/io/basicInputOutput.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/astro/observation_models/observationManager.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/estimation_setup/createObservationManager.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType >
void checkObservationResidualDiscontinuities(
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
        const std::pair< int, int > observableStartAndSize,
        const observation_models::ObservableType observableType )
{
    if( observableType == observation_models::angular_position || observableType == observation_models::euler_angle_313_observable || observableType == observation_models::relative_angular_position )
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualsBlock = residuals.block( observableStartAndSize.first, 0, observableStartAndSize.second, 1 );
        for( int i = 1; i < residualsBlock.rows( ); i++ )
        {
            if( std::fabs( residualsBlock( i, 0 ) - residualsBlock( i - 1, 0 ) ) > 6.0 )
            {
                if( residualsBlock( i, 0 ) > 0 )
                {
                    residualsBlock( i, 0 ) = residualsBlock( i, 0 ) - 2.0 * mathematical_constants::PI;
                }
                else
                {
                    residualsBlock( i, 0 ) = residualsBlock( i, 0 ) + 2.0 * mathematical_constants::PI;
                }
            }
            else if( std::fabs( residualsBlock( i, 0 ) - residualsBlock( i - 1, 0 ) ) > 3.0 )
            {
                std::cerr<<"Warning, detected jump in observation residual of size "<<std::fabs( residualsBlock( i, 0 ) - residualsBlock( i - 1, 0 ) )<<
                         " for observable type "<<observableType<<std::endl;
            }
        }
        residuals.block( observableStartAndSize.first, 0, observableStartAndSize.second, 1 ) = residualsBlock;
    }
}

template< typename ObservationScalarType = double, typename TimeType = double,
    typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateResiduals(
    const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationsCollection,
    const std::map< observation_models::ObservableType,
        std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >& observationSimulator,
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals )
{
    residuals = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( observationsCollection->getTotalObservableSize( ), 1 );

    typename observation_models::ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets
        sortedObservations = observationsCollection->getObservationsSets( );

    // Iterate over all observable types in observationsAndTimes
    for( auto observablesIterator : sortedObservations )
    {
        observation_models::ObservableType currentObservableType = observablesIterator.first;

        // Iterate over all link ends for current observable type in observationsAndTimes
        for( auto dataIterator : observablesIterator.second )
        {
            observation_models::LinkEnds currentLinkEnds = dataIterator.first;
            for( unsigned int i = 0; i < dataIterator.second.size( ); i++ )
            {
                std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > currentObservations =
                    dataIterator.second.at( i );
                std::pair< int, int > observationIndices = observationsCollection->getObservationSetStartAndSize( ).at(
                    currentObservableType ).at( currentLinkEnds ).at( i );

                // Compute estimated ranges and range partials from current parameter estimate.
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector;
                observationSimulator.at( currentObservableType )->
                    computeObservations(
                    currentObservations->getObservationTimes( ), currentLinkEnds,
                    currentObservations->getReferenceLinkEnd( ),
                    currentObservations->getAncilliarySettings( ),
                    observationsVector );

                residuals.block( observationIndices.first, 0, observationIndices.second, 1 ) =
                        ( currentObservations->getObservationsVector( ) - observationsVector );

            }
        }

        std::pair< int, int > observableStartAndSize = observationsCollection->getObservationTypeStartAndSize( ).at( currentObservableType );

        checkObservationResidualDiscontinuities< ObservationScalarType >(
            residuals, observableStartAndSize, currentObservableType );
    }
}


//! Function to calculate the observation partials matrix and residuals
/*!
 *  This function calculates the observation partials matrix and residuals, based on the state transition matrix,
 *  sensitivity matrix and body states resulting from the previous numerical integration iteration.
 *  Partials and observations are calculated by the observationManagers_.
 *  \param observationsAndTimes Observable values and associated time tags, per observable type and set of link ends.
 *  \param parameterVectorSize Length of the vector of estimated parameters
 *  \param totalObservationSize Total number of observations in observationsAndTimes map.
 *  \param residualsAndPartials Pair of residuals of computed w.r.t. input observable values and partials of
 *  observables w.r.t. parameter vector (return by reference).
 */
template< typename ObservationScalarType = double, typename TimeType = double,
    typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateDesignMatrixAndResiduals(
    std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationsCollection,
    const std::map< observation_models::ObservableType,
        std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >& observationManagers,
    const int totalNumberParameters,
    const int totalObservationSize,
    Eigen::MatrixXd& designMatrix,
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
    const bool calculateResiduals = true,
    const bool calculatePartials = true )
{
    if( calculatePartials && totalNumberParameters <= 0 )
    {
        throw std::runtime_error( "Error when computing observation partials; number of parameters is 0 or smaller: " + std::to_string( totalNumberParameters ) );
    }

    // Initialize return data.
    if( calculatePartials )
    {
        designMatrix = Eigen::MatrixXd::Zero( totalObservationSize, totalNumberParameters );
    }

    if( calculateResiduals )
    {
        residuals = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObservationSize, 1 );
    }

    typename observation_models::ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets
        sortedObservations = observationsCollection->getObservationsSets( );

    // Iterate over all observable types in observationsAndTimes
    for( auto observableIt : sortedObservations )
    {
        observation_models::ObservableType currentObservableType = observableIt.first;

        // Iterate over all link ends for current observable type in observationsAndTimes
        for( auto linkEndIt : observableIt.second )
        {
            observation_models::LinkEnds currentLinkEnds = linkEndIt.first;
            for( unsigned int i = 0; i < linkEndIt.second.size( ); i++ )
            {
                std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > currentObservations = linkEndIt.second.at( i );
                std::pair< int, int > observationIndices = observationsCollection->getObservationSetStartAndSize( ).at(
                    currentObservableType ).at( currentLinkEnds ).at( i );

//                std::cout<<"Current size "<<currentObservations->getObservationTimes( ).size( )<<
//                " "<<currentObservations->getObservations( ).size( )<<" "<<observationIndices.first<<" "<<observationIndices.second<<std::endl;
                if( observationIndices.second > 0 )
                {
                    // Compute estimated ranges and range partials from current parameter estimate.
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector;
                    Eigen::MatrixXd partialsMatrix;
                    observationManagers.at( currentObservableType )->
                            computeObservationsWithPartials(
                            currentObservations->getObservationTimes( ),
                            currentLinkEnds,
                            currentObservations->getReferenceLinkEnd( ),
                            currentObservations->getAncilliarySettings( ),
                            observationsVector,
                            partialsMatrix,
                            calculateResiduals,
                            calculatePartials );

                    if( calculatePartials )
                    {
//                        std::cout<<designMatrix.rows( )<<" "<<designMatrix.cols( )<<std::endl;
//                        std::cout<<observationIndices.first<<" "<<0<<" "<<observationIndices.second<<" "<<totalNumberParameters<<std::endl;
//                        std::cout<<partialsMatrix.rows( )<<" "<<partialsMatrix.cols( )<<std::endl<<std::endl;

                        // Set current observation partials in matrix of all partials
                        designMatrix.block( observationIndices.first, 0, observationIndices.second,
                                            totalNumberParameters ) = partialsMatrix;
                    }

                    // Compute residuals for current link ends and observable type.
                    if( calculateResiduals )
                    {
                        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualsVector =
                                currentObservations->getObservationsVector( ) - observationsVector;
                        residuals.block( observationIndices.first, 0, observationIndices.second, 1 ) = ( residualsVector );
                    }
                }
            }
        }

        if( calculateResiduals )
        {
            std::pair< int, int > observableStartAndSize = observationsCollection->getObservationTypeStartAndSize( ).at( currentObservableType );
            checkObservationResidualDiscontinuities< ObservationScalarType >( residuals, observableStartAndSize, currentObservableType );
        }
    }
}

template< typename ObservationScalarType = double, typename TimeType = double,
    typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateDesignMatrix(
    const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationsCollection,
    const std::map< observation_models::ObservableType,
        std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >& observationManagers,
    const int totalNumberParameters,
    const int totalObservationSize,
    Eigen::MatrixXd& designMatrix )
{
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > dummyVector;
    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >(
        observationsCollection, observationManagers, totalNumberParameters, totalObservationSize, designMatrix, dummyVector, false, true );

}


template< typename ObservationScalarType = double, typename TimeType = double,
    typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateResiduals(
    const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationsCollection,
    const std::map< observation_models::ObservableType,
        std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >& observationManagers,
    const int totalObservationSize,
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals )
{
    Eigen::VectorXd dummyMatrix;
    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >(
        observationsCollection, observationManagers, 0, totalObservationSize, dummyMatrix, residuals, true, false );

}

//! Top-level class for performing orbit determination.
/*!
 *  Top-level class for performing orbit determination. All required propagation/estimation settings are provided to
 *  this class, which then creates all objects needed for the propagation and estimation process. The parameter
 *  estimation itself is performed by providing measurement data and related metadata (as EstimationInput) to the estimateParameters
 *  function.
 */
template< typename ObservationScalarType = double, typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class OrbitDeterminationManager
{
public:

    //! Typedef for vector of observations.
    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;

    //! Typedef for vector of parameters.
    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ParameterVectorType;

    //    //! Typedef for observations per link ends, with associated times and reference link end.
    //    typedef std::map< observation_models::LinkEnds, std::pair< ObservationVectorType,
    //    std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > SingleObservableEstimationInputType;

    //    //! Typedef for complete set of observations data, as used in orbit determination
    //    typedef std::map< observation_models::ObservableType, SingleObservableEstimationInputType > EstimationInputType;

    //    //! Typedef for complete set of observations data in alternative form, as used in orbit determination, convertible to EstimationInputType
    //    //! by convertEstimationInput function.
    //    typedef std::map< observation_models::ObservableType,
    //    std::map< observation_models::LinkEnds, std::pair< std::map< TimeType, ObservationScalarType >,
    //    observation_models::LinkEndType > > > AlternativeEstimationInputType;

    std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > preprocessDeprecatedIntegratorSettings(
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
            const std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettings,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const int integratorIndexOffset = 0 )
    {
        std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > independentIntegratorSettingsList =
                utilities::cloneDuplicatePointers( integratorSettings );
        if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< ObservationScalarType, TimeType > >( propagatorSettings ) != nullptr )
        {
            std::shared_ptr< propagators::MultiArcPropagatorSettings< ObservationScalarType, TimeType > > multiArcPropagatorSettings =
                    std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< ObservationScalarType, TimeType > >( propagatorSettings );

            std::vector< double > arcStartTimes = estimatable_parameters::getMultiArcStateEstimationArcStartTimes(
                            parametersToEstimate, ( integratorIndexOffset == 0 ) );
            if( multiArcPropagatorSettings->getSingleArcSettings( ).size( ) != arcStartTimes.size( ) )
            {
                throw std::runtime_error( "Error when processing deprecated integrator/propagator settings in estimation; inconsistent number of arcs" );
            }
            for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
            {
                multiArcPropagatorSettings->getSingleArcSettings( ).at( i )->resetInitialTime( arcStartTimes.at( i ) );
            }
        }
        else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< ObservationScalarType, TimeType > >( propagatorSettings ) != nullptr )
        {
            independentIntegratorSettingsList = preprocessDeprecatedIntegratorSettings(
                        parametersToEstimate, integratorSettings,
                        std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< ObservationScalarType, TimeType > >( propagatorSettings )->getMultiArcPropagatorSettings( ),
                        1 );
        }
        return independentIntegratorSettingsList;
    }

    OrbitDeterminationManager(
            const SystemOfBodies &bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
            parametersToEstimate,
            const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
            const std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettings,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const bool propagateOnCreation = true ):
        parametersToEstimate_( parametersToEstimate )
    {

        std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > processedIntegratorSettings =
                preprocessDeprecatedIntegratorSettings( parametersToEstimate, integratorSettings, propagatorSettings );
        initializeOrbitDeterminationManager(
                    bodies, observationSettingsList, propagators::validateDeprecatePropagatorSettings( processedIntegratorSettings, propagatorSettings ),
                    propagateOnCreation );
    }

    //! Constructor
    /*!
     *  Constructor
     *  \param bodies Map of body objects with names of bodies, storing all environment models used in simulation.
     *  \param parametersToEstimate Container object for all parameters that are to be estimated
     *  \param observationSettingsMap Sets of observation model settings per link ends (i.e. transmitter, receiver, etc.)
     *  per observable type for which measurement data is to be provided in orbit determination process
     *  (through estimateParameters function)
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param propagateOnCreation Boolean denoting whether initial propagatoon is to be performed upon object creation (default
     *  true)
     */
    OrbitDeterminationManager(
            const SystemOfBodies &bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
            parametersToEstimate,
            const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const bool propagateOnCreation = true ):
        parametersToEstimate_( parametersToEstimate ),
        bodies_( bodies )
    {
        std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > processedIntegratorSettings;
        if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< ObservationScalarType, TimeType > >( propagatorSettings ) != nullptr )
        {
            processedIntegratorSettings = { integratorSettings };
        }
        else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< ObservationScalarType, TimeType > >( propagatorSettings ) != nullptr )
        {
            int numberOfArcs = estimatable_parameters::getMultiArcStateEstimationArcStartTimes(
                        parametersToEstimate, true ).size( );
            std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > unprocessedIntegratorSettings =
                    std::vector<std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > >(
                        numberOfArcs, integratorSettings );
            processedIntegratorSettings =
                    preprocessDeprecatedIntegratorSettings( parametersToEstimate, unprocessedIntegratorSettings, propagatorSettings );
        }
        else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< ObservationScalarType, TimeType > >( propagatorSettings ) != nullptr )
        {
            int numberOfArcs = estimatable_parameters::getMultiArcStateEstimationArcStartTimes(
                        parametersToEstimate, false ).size( );
            std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > unprocessedIntegratorSettings =
                    std::vector<std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > >(
                        numberOfArcs + 1, integratorSettings );
            processedIntegratorSettings =
                    preprocessDeprecatedIntegratorSettings( parametersToEstimate, unprocessedIntegratorSettings, propagatorSettings );
        }

        initializeOrbitDeterminationManager(
                    bodies, observationSettingsList, propagators::validateDeprecatePropagatorSettings(
                        processedIntegratorSettings, propagatorSettings ),
                    propagateOnCreation );
    }
//        parametersToEstimate_( parametersToEstimate )
//    {
//        std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > independentIntegratorSettingsList =
//                { integratorSettings };
//        if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< ObservationScalarType, TimeType > >( propagatorSettings ) != nullptr )
//        {
//            std::vector< double > arcStartTimes = estimatable_parameters::getMultiArcStateEstimationArcStartTimes(
//                        parametersToEstimate, true );
//            std::vector<std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettingsList(
//                        arcStartTimes.size( ), integratorSettings);
//            independentIntegratorSettingsList = utilities::deepcopyDuplicatePointers( integratorSettingsList );

//            for( unsigned int i = 0; i < independentIntegratorSettingsList.size( ); i++ )
//            {
//                independentIntegratorSettingsList.at( i )->initialTime_ = arcStartTimes.at( i );
//            }
//        }
//        propagatorSettings =

//        initializeOrbitDeterminationManager( bodies, observationSettingsList, propagators::validateDeprecatePropagatorSettings( integratorSettings, propagatorSettings ),
//                                             propagateOnCreation );
//    }




    OrbitDeterminationManager(
            const SystemOfBodies &bodies,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >
            parametersToEstimate,
            const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const bool propagateOnCreation = true ):
        parametersToEstimate_( parametersToEstimate ),
        considerParameters_( parametersToEstimate_->getConsiderParameters( ) ),
        bodies_( bodies )
    {
        initializeOrbitDeterminationManager( bodies, observationSettingsList, propagatorSettings, propagateOnCreation );
    }

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > getParametersToEstimate( )
    {
        return parametersToEstimate_;
    }

    SystemOfBodies getBodies( )
    {
        return bodies_;
    }


    //! Function to retrieve map of all observation managers
    /*!
     *  Function to retrieve map of all observation managers. A single observation manager can simulate observations and
     *  calculate observation partials for all link ends involved in the given observable type.
     *  \return Map of observation managers for all observable types involved in current orbit determination.
     */
    std::map< observation_models::ObservableType,
    std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >
    getObservationManagers( ) const
    {
        return observationManagers_;
    }

    //! Function to retrieve map of all observation simulators
    /*!
     *  Function to retrieve map of all observation simulators. A single observation simulators can simulate observations all
     *  link ends involved in the given observable type. The observation simulators are retrieved from the observation manager
     *  objects (that are stored in the observationManagers_ map).
     *  \return Map of observation simulators for all observable types involved in current orbit determination.
     */
    std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >
    getObservationSimulators( ) const
    {
        std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators;

        for( typename std::map< observation_models::ObservableType,
             std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >::const_iterator
             managerIterator = observationManagers_.begin( ); managerIterator != observationManagers_.end( );
             managerIterator++ )
        {
            observationSimulators.push_back( managerIterator->second->getObservationSimulator( ) );
        }

        return observationSimulators;
    }

    Eigen::MatrixXd normalizeAprioriCovariance(
            const Eigen::MatrixXd& inverseAPrioriCovariance,
            const Eigen::VectorXd& normalizationValues )
    {
        int numberOfEstimatedParameters = inverseAPrioriCovariance.rows( );
        Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix = Eigen::MatrixXd::Zero(
                    numberOfEstimatedParameters, numberOfEstimatedParameters );

        for( int j = 0; j < numberOfEstimatedParameters; j++ )
        {
            for( int k = 0; k < numberOfEstimatedParameters; k++ )
            {
                normalizedInverseAprioriCovarianceMatrix( j, k ) = inverseAPrioriCovariance( j, k ) /
                        ( normalizationValues( j ) * normalizationValues( k ) );
            }
        }
        return normalizedInverseAprioriCovarianceMatrix;
    }

    Eigen::MatrixXd normalizeCovariance(
            const Eigen::MatrixXd& covariance,
            const Eigen::VectorXd& normalizationFactors )
    {
        int numberParameters = covariance.rows( );
        Eigen::MatrixXd normalizedCovariance = Eigen::MatrixXd::Zero( numberParameters, numberParameters );
        for( int j = 0; j < numberParameters; j++ )
        {
            for( int k = 0; k < numberParameters; k++ )
            {
                normalizedCovariance( j, k ) = covariance( j, k ) * ( normalizationFactors( j ) * normalizationFactors( k ) );
            }
        }
        return normalizedCovariance;
    }

    //! Function to normalize the matrix of partial derivatives so that each column is in the range [-1,1]
    /*!
     * Function to normalize the matrix of partial derivatives so that each column is in the range [-1,1]
     * \param observationMatrix Matrix of partial derivatives. Matrix modified by this function, and normalized matrix is
     * returned by reference
     * \return Vector with scaling values used for normalization
     */
    Eigen::VectorXd normalizeDesignMatrix( Eigen::MatrixXd& observationMatrix )
    {
        Eigen::VectorXd normalizationTerms = Eigen::VectorXd( observationMatrix.cols( ) );

        for( int i = 0; i < observationMatrix.cols( ); i++ )
        {
            Eigen::VectorXd currentVector = observationMatrix.block( 0, i, observationMatrix.rows( ), 1 );
            double minimum = currentVector.minCoeff( );
            double maximum = currentVector.maxCoeff( );
            if( std::fabs( minimum ) > maximum )
            {
                normalizationTerms( i ) = minimum;
            }
            else
            {
                normalizationTerms( i ) = maximum;
            }
            if( normalizationTerms( i ) == 0.0 )
            {
                normalizationTerms( i ) = 1.0;
            }
            currentVector = currentVector / normalizationTerms( i );

            observationMatrix.block( 0, i, observationMatrix.rows( ), 1 ) = currentVector;
        }

        //        for( unsigned int i = 0; i < observationLinkParameterIndices_.size( ); i++ )
        //        {
        //            int currentColumn = observationLinkParameterIndices_.at( i );
        //            int startIndex = -1;
        //            int endIndex = -1;
        //            std::vector< double > partialMaximum;
        //            bool isInRange = 0;

        //            for( int j = 0; j < observationMatrix.rows( ); j++ )
        //            {
        //                if( observationMatrix( j, currentColumn ) != 0.0 )
        //                {
        //                   if( isInRange == 0 )
        //                   {
        //                       isInRange = 1;
        //                       startIndex = j;
        //                   }
        //                }
        //                else if( ( startIndex != -1 ) && isInRange && ( observationMatrix( j, currentColumn ) == 0.0 ) )
        //                {
        //                    isInRange = 0;
        //                    endIndex = j;
        //                }

        //            }
        //        }
        return normalizationTerms;
    }

//    void saveResultsFromCurrentIteration(
//            std::vector< std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > >& dynamicsHistoryPerIteration )
//    {
//        if( std::dynamic_pointer_cast< propagators::HybridArcVariationalEquationsSolver< ObservationScalarType, TimeType > >(  variationalEquationsSolver_) != nullptr )
//        {

//            std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< ObservationScalarType, TimeType > > hybridArcSolver =
//                    std::dynamic_pointer_cast< propagators::HybridArcVariationalEquationsSolver< ObservationScalarType, TimeType > >(  variationalEquationsSolver_);

//            std::vector< std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > currentDynamicsSolution;
//            std::vector< std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > currentDependentVariableSolution;

//            currentDynamicsSolution = hybridArcSolver->getSingleArcSolver( )->getDynamicsSimulator( )->getEquationsOfMotionNumericalSolutionBase( );
//            currentDependentVariableSolution = hybridArcSolver->getSingleArcSolver( )->getDynamicsSimulator( )->getDependentVariableNumericalSolutionBase( );

//            auto multiArcDynamicsSolution =
//                    hybridArcSolver->getMultiArcSolver( )->getDynamicsSimulatorBase( )->getEquationsOfMotionNumericalSolutionBase( );
//            auto multiArcDependentVariableSolution =
//                    hybridArcSolver->getMultiArcSolver( )->getDynamicsSimulatorBase( )->getDependentVariableNumericalSolutionBase( );
//            currentDynamicsSolution.insert( currentDynamicsSolution.end( ), multiArcDynamicsSolution.begin( ), multiArcDynamicsSolution.end( ) );
//            currentDependentVariableSolution.insert( currentDependentVariableSolution.end( ), multiArcDependentVariableSolution.begin( ), multiArcDependentVariableSolution.end( ) );

//            dynamicsHistoryPerIteration.push_back( currentDynamicsSolution );
//            dependentVariableHistoryPerIteration.push_back( currentDependentVariableSolution );
//        }
//        else
//        {
//            dynamicsHistoryPerIteration.push_back(
//                        variationalEquationsSolver_->getDynamicsSimulatorBase( )->getEquationsOfMotionNumericalSolutionBase( ));
//            dependentVariableHistoryPerIteration.push_back(
//                        variationalEquationsSolver_->getDynamicsSimulatorBase( )->getDependentVariableNumericalSolutionBase( ));

//        }
//    }

    std::shared_ptr< CovarianceAnalysisOutput< ObservationScalarType, TimeType > > computeCovariance(
            const std::shared_ptr< CovarianceAnalysisInput< ObservationScalarType, TimeType > > estimationInput )
    {
        // Get total number of observations
        int totalNumberOfObservations = estimationInput->getObservationCollection( )->getTotalObservableSize( );

        // Define full parameters values
        ParameterVectorType parameterValues = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );
        ParameterVectorType fullParameterEstimate;
        fullParameterEstimate.resize( totalNumberParameters_ );
        fullParameterEstimate.segment( 0, numberEstimatedParameters_ ) = parameterValues;
        if ( considerParametersIncluded_ )
        {
            fullParameterEstimate.segment( numberEstimatedParameters_, numberConsiderParameters_ ) = considerParametersValues_;
        }

        // Compute design matrices (estimated and consider), and residuals (empty for covariance analysis)
        bool exceptionDuringPropagation = false;
        std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > simulationResults;
        std::pair< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > designMatricesAndResiduals =
                performPreEstimationSteps( estimationInput, fullParameterEstimate, false, 0, exceptionDuringPropagation, simulationResults );
        Eigen::MatrixXd designMatrixEstimatedParameters = designMatricesAndResiduals.first.first;
        Eigen::MatrixXd designMatrixConsiderParameters;
        if ( considerParametersIncluded_ )
        {
            designMatrixConsiderParameters = designMatricesAndResiduals.first.second;
        }
        else
        {
            designMatrixConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
        }

        // Normalise partials and inverse a priori covariance
        Eigen::VectorXd normalizationTerms = normalizeDesignMatrix( designMatrixEstimatedParameters );
        Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix = normalizeAprioriCovariance(
                estimationInput->getInverseOfAprioriCovariance( numberEstimatedParameters_ ), normalizationTerms );

        // Normalise partials w.r.t. consider parameters and consider covariance
        Eigen::VectorXd considerNormalizationTerms;
        Eigen::MatrixXd normalizedConsiderCovariance;
        if ( considerParametersIncluded_ )
        {
            considerNormalizationTerms = normalizeDesignMatrix( designMatrixConsiderParameters );
            normalizedConsiderCovariance = normalizeCovariance( estimationInput->getConsiderCovariance( ), considerNormalizationTerms );
        }
        else
        {
            considerNormalizationTerms = Eigen::VectorXd::Zero( 0 );
            normalizedConsiderCovariance = Eigen::MatrixXd::Zero( 0, 0 );
        }


        // Retrieve constraints
        Eigen::MatrixXd constraintStateMultiplier;
        Eigen::VectorXd constraintRightHandSide;
        parametersToEstimate_->getConstraints( constraintStateMultiplier, constraintRightHandSide );

        // Compute inverse of updated covariance
        Eigen::MatrixXd inverseNormalizedCovariance = linear_algebra::calculateInverseOfUpdatedCovarianceMatrix(
                designMatrixEstimatedParameters.block( 0, 0, designMatrixEstimatedParameters.rows( ), numberEstimatedParameters_ ),
                estimationInput->getWeightsMatrixDiagonals( ),
                normalizedInverseAprioriCovarianceMatrix, constraintStateMultiplier, constraintRightHandSide, estimationInput->getLimitConditionNumberForWarning( ) );

        // Compute contribution consider parameters
        Eigen::MatrixXd covarianceContributionConsiderParameters;
        if ( considerParametersIncluded_ )
        {
            covarianceContributionConsiderParameters = linear_algebra::calculateConsiderParametersCovarianceContribution(
                    inverseNormalizedCovariance.inverse( ), designMatrixEstimatedParameters, estimationInput->getWeightsMatrixDiagonals( ),
                    designMatrixConsiderParameters, normalizedConsiderCovariance );
        }
        else
        {
            covarianceContributionConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
        }

        // Create covariance output object
        std::shared_ptr< CovarianceAnalysisOutput< ObservationScalarType, TimeType > > estimationOutput =
                std::make_shared< CovarianceAnalysisOutput< ObservationScalarType, TimeType > >(
                     designMatrixEstimatedParameters, estimationInput->getWeightsMatrixDiagonals( ), normalizationTerms,
                    inverseNormalizedCovariance, designMatrixConsiderParameters, considerNormalizationTerms, covarianceContributionConsiderParameters,
                    exceptionDuringPropagation );

        return estimationOutput;
    }

    //! Function to perform parameter estimation from measurement data.
    /*!
     *  Function to perform parameter estimation, including orbit determination, i.e. body initial states, from measurement data.
     *  All observable types and link ends per obsevable types that are included in the measurement data input must have been
     *  provided to the constructor by the observationSettingsMap parameter.
     *  \param estimationInput Object containing all measurement data, associated metadata, including measurement weight, and a priori
     *  estimate for covariance matrix and parameter adjustment.
     *  \param convergenceChecker Object used to check convergence/termination of algorithm
     *  \return Object containing estimated parameter value and associateed data, such as residuals and observation partials.
     */
    std::shared_ptr< EstimationOutput< ObservationScalarType, TimeType > > estimateParameters(
            std::shared_ptr< EstimationInput< ObservationScalarType, TimeType > > estimationInput )

    {
        currentParameterEstimate_ = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );

        // Get number of observations
        int totalNumberOfObservations = estimationInput->getObservationCollection( )->getTotalObservableSize( );

        if( estimationInput->getWeightsMatrixDiagonals( ).rows( ) != totalNumberOfObservations )
        {
            throw std::runtime_error( "Error when estimating parameters, size of weights diagonal (" +
                std::to_string( estimationInput->getWeightsMatrixDiagonals( ).rows( ) ) + ") is not compatible with number of observations (" +
                std::to_string( totalNumberOfObservations ) + ")" );
        }
        // Declare variables to be returned (i.e. results from best iteration)
        double bestResidual = TUDAT_NAN;
        ParameterVectorType bestParameterEstimate = ParameterVectorType::Constant( numberEstimatedParameters_, TUDAT_NAN );
        Eigen::VectorXd bestTransformationData = Eigen::VectorXd::Constant( numberEstimatedParameters_, TUDAT_NAN );
        Eigen::VectorXd bestResiduals = Eigen::VectorXd::Constant( totalNumberOfObservations, TUDAT_NAN );
        Eigen::MatrixXd bestDesignMatrixEstimatedParameters = Eigen::MatrixXd::Constant( totalNumberOfObservations, totalNumberParameters_, TUDAT_NAN );
        Eigen::VectorXd bestWeightsMatrixDiagonal = Eigen::VectorXd::Constant( totalNumberOfObservations, TUDAT_NAN );
        Eigen::MatrixXd bestInverseNormalizedCovarianceMatrix = Eigen::MatrixXd::Constant( numberEstimatedParameters_, numberEstimatedParameters_, TUDAT_NAN );


        Eigen::VectorXd bestConsiderTransformationData;
        Eigen::MatrixXd bestDesignMatrixConsiderParameters, bestConsiderCovarianceContribution;
        if ( considerParametersIncluded_ )
        {
            bestConsiderTransformationData = Eigen::VectorXd::Constant( numberConsiderParameters_, TUDAT_NAN );
            bestDesignMatrixConsiderParameters = Eigen::MatrixXd::Constant( totalNumberOfObservations, numberConsiderParameters_, TUDAT_NAN );
            bestConsiderCovarianceContribution = Eigen::MatrixXd::Constant( numberEstimatedParameters_, numberEstimatedParameters_, TUDAT_NAN );
        }
        else
        {
            bestConsiderTransformationData = Eigen::VectorXd::Zero( 0 );
            bestDesignMatrixConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
            bestConsiderCovarianceContribution = Eigen::MatrixXd::Zero( 0, 0 );
        }

        std::vector< Eigen::VectorXd > residualHistory;
        std::vector< ParameterVectorType > parameterHistory;
        std::vector< std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > > simulationResultsPerIteration;

        // Declare residual bookkeeping variables
        std::vector< double > rmsResidualHistory;
        double residualRms;

        // Set current parameter estimate as both previous and current estimate
        ParameterVectorType newParameterEstimate = currentParameterEstimate_;
        ParameterVectorType oldParameterEstimate = currentParameterEstimate_;
        ParameterVectorType newFullParameterEstimate;
        newFullParameterEstimate.resize( totalNumberParameters_ );

        bool exceptionDuringPropagation = false, exceptionDuringInversion = false;

        // Iterate until convergence (at least once)
        int bestIteration = -1;
        int numberOfIterations = 0;
        while( true )
        {
            oldParameterEstimate = newParameterEstimate;
            newFullParameterEstimate.segment( 0, numberEstimatedParameters_ ) = newParameterEstimate;
            if ( considerParametersIncluded_ )
            {
                newFullParameterEstimate.segment( numberEstimatedParameters_, numberConsiderParameters_ ) = considerParametersValues_;
            }

            // Compute design matrices (for estimated and consider parameters) and residuals.
            std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > > simulationResults;
            std::pair< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > designMatricesAndResiduals =
                    performPreEstimationSteps( estimationInput, newFullParameterEstimate, true, numberOfIterations, exceptionDuringPropagation, simulationResults );
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residuals = designMatricesAndResiduals.second;
            Eigen::MatrixXd designMatrixEstimatedParameters = designMatricesAndResiduals.first.first;
            Eigen::MatrixXd designMatrixConsiderParameters;
            if ( considerParametersIncluded_ )
            {
                designMatrixConsiderParameters = designMatricesAndResiduals.first.second;
            }
            else
            {
                designMatrixConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
            }

            // Set simulation results
            if( estimationInput->getSaveStateHistoryForEachIteration( ) )
            {
                simulationResultsPerIteration.push_back( simulationResults );
            }

            // Normalise estimated parameters partials and inverse apriori covariance
            Eigen::VectorXd normalizationTerms = normalizeDesignMatrix( designMatrixEstimatedParameters );
            Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix = normalizeAprioriCovariance(
                    estimationInput->getInverseOfAprioriCovariance( numberEstimatedParameters_ ), normalizationTerms );

            // Normalise partials w.r.t. consider parameters, consider covariance and parameters deviations
            Eigen::VectorXd normalizationTermsConsider, normalizedConsiderParametersDeviation;
            Eigen::MatrixXd normalizedConsiderCovariance;
            if ( considerParametersIncluded_ )
            {
                normalizationTermsConsider = normalizeDesignMatrix( designMatrixConsiderParameters );
                normalizedConsiderCovariance = normalizeCovariance( estimationInput->getConsiderCovariance( ), normalizationTermsConsider );
                normalizedConsiderParametersDeviation = estimationInput->considerParametersDeviations_.cwiseProduct( normalizationTermsConsider );
            }
            else
            {
                normalizationTermsConsider = Eigen::VectorXd::Zero( 0 );
                normalizedConsiderCovariance = Eigen::MatrixXd::Zero( 0, 0 );
                normalizedConsiderParametersDeviation = Eigen::VectorXd::Zero( 0 );
            }

            // Perform least squares calculation for correction to parameter vector.
            std::pair< Eigen::VectorXd, Eigen::MatrixXd > leastSquaresOutput;
            try
            {
                // Get constraints
                Eigen::MatrixXd constraintStateMultiplier;
                Eigen::VectorXd constraintRightHandSide;
                parametersToEstimate_->getConstraints( constraintStateMultiplier, constraintRightHandSide );

                double conditionNumberCheck = estimationInput->getLimitConditionNumberForWarning( );
                if( numberOfIterations > 0 && estimationInput->conditionNumberWarningEachIteration_ == false )
                {
                    conditionNumberCheck = TUDAT_NAN;
                }
                // Perform LSQ inversion
                leastSquaresOutput = std::move( linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix(
                        designMatrixEstimatedParameters, residuals.template cast< double >( ), estimationInput->getWeightsMatrixDiagonals( ),
                        normalizedInverseAprioriCovarianceMatrix, conditionNumberCheck, constraintStateMultiplier, constraintRightHandSide,
                        designMatrixConsiderParameters, normalizedConsiderParametersDeviation ) );

                if( constraintStateMultiplier.rows( ) > 0 )
                {
                    leastSquaresOutput.first.conservativeResize( numberEstimatedParameters_ );
                }
            }
            catch( std::runtime_error& error )
            {
                std::cerr<<"Error when solving normal equations during parameter estimation: "<<std::endl<<error.what( )<<
                           std::endl<<"Terminating estimation"<<std::endl;
                exceptionDuringInversion = true;
                break;
            }

            ParameterVectorType parameterAddition =
                    ( leastSquaresOutput.first.cwiseQuotient( normalizationTerms.segment( 0, numberEstimatedParameters_ ) ) ).template cast< ObservationScalarType >( );

            // Compute contribution consider parameters
            Eigen::MatrixXd covarianceContributionConsiderParameters;
            if ( considerParametersIncluded_ )
            {
                covarianceContributionConsiderParameters = linear_algebra::calculateConsiderParametersCovarianceContribution(
                        ( leastSquaresOutput.second ).inverse( ), designMatrixEstimatedParameters, estimationInput->getWeightsMatrixDiagonals( ),
                        designMatrixConsiderParameters, normalizedConsiderCovariance );
            }
            else
            {
                covarianceContributionConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
            }

            // Calculate mean residual for current iteration.
            residualRms = linear_algebra::getVectorEntryRootMeanSquare( residuals.template cast< double >( ) );
            rmsResidualHistory.push_back( residualRms );

            if( estimationInput->getSaveResidualsAndParametersFromEachIteration( ) )
            {
                residualHistory.push_back( residuals.template cast< double >( ) );
                if ( numberOfIterations == 0 )
                {
                    parameterHistory.push_back( oldParameterEstimate );
                }
            }

            if( estimationInput->getPrintOutput( ) )
            {
                std::cout << "Current residual: " << residualRms << std::endl;
            }

            // If current iteration is better than previous one, update 'best' data.
            if( residualRms < bestResidual || !( bestResidual == bestResidual ) )
            {
                bestResidual = residualRms;
                bestParameterEstimate = oldParameterEstimate;
                bestResiduals = std::move( residuals.template cast< double >( ) );
                estimationInput->getObservationCollection( )->setResiduals( residuals );
                if( estimationInput->getSaveDesignMatrix( ) )
                {
                    bestDesignMatrixEstimatedParameters = std::move( designMatrixEstimatedParameters );
                    if ( considerParametersIncluded_ )
                    {
                        bestDesignMatrixConsiderParameters = std::move( designMatrixConsiderParameters );
                    }
                }
                bestWeightsMatrixDiagonal = std::move( estimationInput->getWeightsMatrixDiagonals( ) );
                bestTransformationData = std::move( normalizationTerms );
                bestInverseNormalizedCovarianceMatrix = std::move( leastSquaresOutput.second );
                bestIteration = numberOfIterations;

                if ( considerParametersIncluded_ )
                {
                    bestConsiderTransformationData = std::move( normalizationTermsConsider );
                    bestConsiderCovarianceContribution = covarianceContributionConsiderParameters;
                }
            }


            // Increment number of iterations
            numberOfIterations++;

            // Check for convergence
            bool applyParameterCorrection = true;
            bool terminateLoop = false;
            if(  estimationInput->getConvergenceChecker( )->isEstimationConverged( numberOfIterations, rmsResidualHistory ) )
            {
                terminateLoop = true;
                applyParameterCorrection = estimationInput->applyFinalParameterCorrection_;
            }

            if( applyParameterCorrection )
            {
                // Update value of parameter vector
                newParameterEstimate = oldParameterEstimate + parameterAddition;
                parametersToEstimate_->template resetParameterValues<ObservationScalarType>( newParameterEstimate );
                newParameterEstimate = parametersToEstimate_->template getFullParameterValues<ObservationScalarType>( );

                if ( estimationInput->getSaveResidualsAndParametersFromEachIteration( ) )
                {
                    parameterHistory.push_back( newParameterEstimate );
                }

                if ( estimationInput->getPrintOutput( ) )
                {
                    std::cout << "Parameter update" << parameterAddition.transpose( ) << std::endl;
                }
            }

            if( terminateLoop )
            {
                break;
            }
        }

        if( estimationInput->getPrintOutput( ) )
        {
            std::cout << "Final residual: " << bestResidual << std::endl;
        }

        // Create estimation output object
        std::shared_ptr< EstimationOutput< ObservationScalarType, TimeType > > estimationOutput =
                std::make_shared< EstimationOutput< ObservationScalarType, TimeType > >(
                    bestParameterEstimate, bestResiduals, bestDesignMatrixEstimatedParameters, bestWeightsMatrixDiagonal,
                    bestTransformationData, bestInverseNormalizedCovarianceMatrix, bestResidual, bestIteration,
                    residualHistory, parameterHistory, bestDesignMatrixConsiderParameters, bestConsiderTransformationData,
                    bestConsiderCovarianceContribution, exceptionDuringInversion, exceptionDuringPropagation );

        if( estimationInput->getSaveStateHistoryForEachIteration( ) )
        {
            estimationOutput->setSimulationResults( simulationResultsPerIteration );
        }

        return estimationOutput;
    }


    template< int ObservationSize >
    void computePartialsAndObservations(
        const observation_models::LinkEnds& linkEnds,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings,
        const observation_models::ObservableType observableType,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< TimeType >& times,
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationsVector,
        Eigen::MatrixXd& partials )
    {
        std::shared_ptr< observation_models::ObservationManager< ObservationSize, ObservationScalarType, TimeType > > observationManager =
            std::dynamic_pointer_cast< observation_models::ObservationManager< ObservationSize, ObservationScalarType, TimeType > >( observationManagers_.at(
            observableType ) );

        // Extract observation model
        std::shared_ptr< observation_models::ObservationModel< ObservationSize, ObservationScalarType, TimeType > > observationModel =
                observationManager->getObservationModel( linkEnds );

        // Compute analytical partials
        observationManager->computeObservationsWithPartials(
            times, linkEnds, referenceLinkEnd, ancilliarySettings, observationsVector, partials, true, true );

    }
    //! Function to reset the current parameter estimate.
    /*!
     *  Function to reset the current parameter estimate; reintegrates the variational equations and equations of motion with new estimate.
     *  \param newParameterEstimate New estimate of parameter vector.
     *  \param reintegrateVariationalEquations Boolean denoting whether the variational equations are to be reintegrated
     */
    void resetParameterEstimate( const ParameterVectorType& newParameterEstimate, const bool reintegrateVariationalEquations = 1 )
    {
        if( integrateAndEstimateOrbit_ )
        {
            variationalEquationsSolver_->resetParameterEstimate( newParameterEstimate, reintegrateVariationalEquations );
        }
        else
        {
            parametersToEstimate_->template resetParameterValues< ObservationScalarType>( newParameterEstimate );
        }
        currentParameterEstimate_ = newParameterEstimate;
    }

    //    //! Function to convert from one representation of all measurement data to the other
    //    /*!
    //     *  Function to convert from one representation of all measurement data (AlternativeEstimationInputType) to the other (EstimationInputType).
    //     *  In the former, the vector of times and vector of associated observations are stored separately.  In the latter, they are
    //     *  stored as a map with time as key and observation as value.
    //     *  \param alternativeEstimationInput Original representation of measurement data
    //     *  \return Converted measurement data.
    //     */
    //    static EstimationInputType convertEstimationInput( const AlternativeEstimationInputType& alternativeEstimationInput )
    //    {
    //        // Declare return data structure.
    //        EstimationInputType convertedObservationAndTimes;

    //        // Iterate over all observable types.
    //        for( typename AlternativeEstimationInputType::const_iterator inputIterator = alternativeEstimationInput.begin( );
    //             inputIterator != alternativeEstimationInput.end( ); inputIterator++ )
    //        {
    //            // Declare data structure for converted measurement data at single observable type.
    //            std::map< observation_models::LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > singleTypeObservations;

    //            // Iterate over all link ends in current observable types
    //            for( typename std::map< observation_models::LinkEnds, std::pair< std::map< TimeType, ObservationScalarType >, observation_models::LinkEndType > >::const_iterator stationIterator =
    //                 inputIterator->second.begin( ); stationIterator != inputIterator->second.end( ); stationIterator++ )
    //            {
    //                // Initialize vector of observation times.
    //                std::vector< TimeType > times;
    //                times.resize( stationIterator->second.first.size( ) );

    //                // Initialize vector of observation values.
    //                ObservationVectorType observations = ObservationVectorType::Zero( stationIterator->second.first.size( ) );

    //                // Iterate over all observations in input map, and split time and observation.
    //                int counter = 0;
    //                for( typename std::map< TimeType, ObservationScalarType >::const_iterator singleDataSetIterator = stationIterator->second.first.begin( );
    //                     singleDataSetIterator != stationIterator->second.first.end( ); singleDataSetIterator++ )
    //                {
    //                    times[ counter ] = singleDataSetIterator->first;
    //                    observations( counter ) = singleDataSetIterator->second;
    //                    counter++;
    //                }

    //                // Set converted data for current link nends and observable.
    //                singleTypeObservations[ stationIterator->first ] = std::make_pair( observations, std::make_pair( times, stationIterator->second.second ) );
    //            }

    //            // Set converted data for current observable.
    //            convertedObservationAndTimes[ inputIterator->first ] = singleTypeObservations;
    //        }
    //        return convertedObservationAndTimes;
    //    }

    //! Function to retrieve the object to numerical integrate and update the variational equations and equations of motion
    /*!
     *  Function to retrieve the object to numerical integrate and update the variational equations and equations of motion
     *  \return Object to numerical integrate and update the variational equations and equations of motion
     */
    std::shared_ptr< propagators::VariationalEquationsSolver< ObservationScalarType, TimeType > >
    getVariationalEquationsSolver( ) const
    {
        return variationalEquationsSolver_;
    }

    //! Function to retrieve an observation manager
    /*!
     *  Function to retrieve an observation manager for a single observable type. The observation manager can simulate observations and
     *  calculate observation partials for all link ends involved in the given observable type.
     *  \param observableType Type of observable for which manager is to be retrieved.
     *  \return Observation manager for given observable type.
     */
    std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > getObservationManager(
            const observation_models::ObservableType observableType ) const
    {
        // Check if manager exists for requested observable type.
        if( observationManagers_.count( observableType ) == 0 )
        {
            throw std::runtime_error(
                        "Error when retrieving observation manager of type " + std::to_string(
                            observableType ) + ", manager not found" );
        }

        return observationManagers_.at( observableType );
    }

    //! Function to retrieve the current paramater estimate.
    /*!
     *  Function to retrieve the current paramater estimate.
     *  \return Current paramater estimate.
     */
    ParameterVectorType getCurrentParameterEstimate( )
    {
        return currentParameterEstimate_;
    }

    //! Function to retrieve the object used to propagate/process the solution of the variational equations/dynamics.
    /*!
     *  Function to retrieve the object used to propagate/process the numerical solution of the variational.
     *  equations/dynamics.
     *  \return Object used to propagate/process the numerical solution of the variational equations/dynamics.
     */
    std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface >
    getStateTransitionAndSensitivityMatrixInterface( )
    {
        return stateTransitionAndSensitivityMatrixInterface_;
    }

protected:

    //! Function called by either constructor to initialize the object.
    /*!
     *  Function called by either constructor to initialize the object.
     *  \param bodies Map of body objects with names of bodies, storing all environment models used in simulation.
     *  \param observationSettingsMap Sets of observation model settings per link ends (i.e. transmitter, receiver, etc.)
     *  for which measurement data is to be provided in orbit determination process
     *  (through estimateParameters function)
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param propagateOnCreation Boolean denoting whether initial propagatoon is to be performed upon object creation (default
     *  true)
     */
    void initializeOrbitDeterminationManager(
            const SystemOfBodies &bodies,
            const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const bool propagateOnCreation = true )
    {
        propagators::toggleIntegratedResultSettings< ObservationScalarType, TimeType >( propagatorSettings );
        using namespace numerical_integrators;
        using namespace orbit_determination;
        using namespace observation_models;

        // Detect whether consider parameters are included
        considerParametersIncluded_ = false;
        if ( considerParameters_ != nullptr )
        {
            considerParametersIncluded_ = true;
        }

        // Create full set of parameters (estimated + consider parameters combined), and define corresponding indices
        setFullParametersSet( );
        getEstimatedAndConsiderParametersIndices( );

        // Retrieve size of estimated and consider parameters
        totalNumberParameters_ = fullParameters_->getParameterSetSize( );
        numberEstimatedParameters_ = parametersToEstimate_->getParameterSetSize( );
        numberConsiderParameters_ = 0;
        if ( considerParameters_ != nullptr )
        {
            numberConsiderParameters_ = considerParameters_->getParameterSetSize( );
        }

        // Check if any dynamics is to be estimated
        std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > > initialDynamicalStates =
                estimatable_parameters::getListOfInitialDynamicalStateParametersEstimate< ObservationScalarType >( fullParameters_ );
        if( initialDynamicalStates.size( ) > 0 )
        {
            integrateAndEstimateOrbit_ = true;
        }
        else
        {
            integrateAndEstimateOrbit_ = false;
        }

        propagatorSettings->getOutputSettingsBase( )->setUpdateDependentVariableInterpolator( true );
        if( integrateAndEstimateOrbit_ )
        {
            variationalEquationsSolver_ = simulation_setup::createVariationalEquationsSolver< ObservationScalarType, TimeType >(
                    bodies, propagatorSettings, fullParameters_, propagateOnCreation );
        }

        if( integrateAndEstimateOrbit_ )
        {
            stateTransitionAndSensitivityMatrixInterface_ = variationalEquationsSolver_->getStateTransitionMatrixInterface( );
        }
        else if( propagatorSettings == nullptr )
        {
            stateTransitionAndSensitivityMatrixInterface_ = createStateTransitionAndSensitivityMatrixInterface< ObservationScalarType, TimeType >(
                        propagatorSettings, fullParameters_, 0, totalNumberParameters_ );
        }
        else
        {
            throw std::runtime_error( "Error, cannot parse propagator settings without estimating dynamics in OrbitDeterminationManager" );
        }

        // TODO correct this when moving dependent variable interface into results object
        if( std::dynamic_pointer_cast< propagators::HybridArcVariationalEquationsSolver< ObservationScalarType, TimeType > >( variationalEquationsSolver_ ) == nullptr )
        {
            dependentVariablesInterface_ = variationalEquationsSolver_->getDynamicsSimulatorBase( )->getDependentVariablesInterface( );
        }


        // Iterate over all observables and create observation managers.
        observationManagers_ = createObservationManagersBase(
            observationSettingsList, bodies, fullParameters_,
            stateTransitionAndSensitivityMatrixInterface_, dependentVariablesInterface_ );

        // Set current parameter estimate from body initial states and parameter set.
        currentParameterEstimate_ = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );
        currentFullParameterValues_ = fullParameters_->template getFullParameterValues< ObservationScalarType >( );
        if ( considerParametersIncluded_ )
        {
            considerParametersValues_ = considerParameters_->template getFullParameterValues<ObservationScalarType>( );
        }
        else
        {
            considerParametersValues_ = ParameterVectorType::Zero( 0 );
        }

    }

    //! Function to create full parameters set with estimated and consider parameters.
    void setFullParametersSet( )
    {
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > fullDoubleParameters = parametersToEstimate_->getEstimatedDoubleParameters( );
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > fullVectorParameters = parametersToEstimate_->getEstimatedVectorParameters( );

        // Check if consider parameters are included in full set of parameters
        if ( considerParametersIncluded_ )
        {
            std::vector< std::string > parametersDescriptions = parametersToEstimate_->getParametersDescriptions( );
            std::vector< std::string > considerParametersDescriptions = considerParameters_->getParametersDescriptions( );
            for ( unsigned int i = 0 ; i < considerParametersDescriptions.size( ) ; i++ )
            {
                if (std::find(parametersDescriptions.begin(), parametersDescriptions.end(), considerParametersDescriptions[i]) != parametersDescriptions.end()) {
                    throw std::runtime_error("Error when initialising orbit determination manager, the following consider parameter is already included as estimated parameter: "
                                             + considerParametersDescriptions[i]);
                }
            }

            for ( unsigned int i = 0 ; i < considerParameters_->getEstimatedDoubleParameters( ).size( ) ; i++ )
            {
                fullDoubleParameters.push_back( considerParameters_->getEstimatedDoubleParameters( )[ i ] );
            }
            for ( unsigned int i = 0 ; i < considerParameters_->getEstimatedVectorParameters( ).size( ) ; i++ )
            {
                fullVectorParameters.push_back( considerParameters_->getEstimatedVectorParameters( )[ i ] );
            }
            if ( considerParameters_->getEstimatedInitialStateParameters( ).size( ) != 0 )
            {
                throw std::runtime_error( "Error when initialising orbit determination manager, consider parameters cannot include initial states parameters." );
            }
        }

        fullParameters_ = std::make_shared< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > >( fullDoubleParameters, fullVectorParameters,
                                                                                                                        parametersToEstimate_->getEstimatedInitialStateParameters( ) );
    }

    void getEstimatedAndConsiderParametersIndices( )
    {
        indicesAndSizeConsiderParameters_.clear( );
        indicesAndSizeEstimatedParameters_.clear( );
        if ( considerParametersIncluded_ )
        {
            std::vector< std::string > considerParametersDescriptions = considerParameters_->getParametersDescriptions();
            for ( unsigned int i = 0; i < considerParametersDescriptions.size( ); i++ )
            {
                std::pair< int, int > indicesInFullParametersSet = fullParameters_->getIndicesForParameterDescription( considerParametersDescriptions[ i ] );
                std::pair< int, int > indicesInConsiderParametersSet = considerParameters_->getIndicesForParameterDescription( considerParametersDescriptions[ i ] );
                indicesAndSizeConsiderParameters_.push_back( std::make_pair( std::make_pair( indicesInConsiderParametersSet.first, indicesInFullParametersSet.first ),
                                                                             indicesInFullParametersSet.second ) );
            }
        }

        std::vector< std::string > estimatedParametersDescriptions = parametersToEstimate_->getParametersDescriptions( );
        for ( unsigned int i = 0 ; i < estimatedParametersDescriptions.size( ) ; i++ )
        {
            std::pair< int, int > indicesInFullParametersSet = fullParameters_->getIndicesForParameterDescription( estimatedParametersDescriptions[ i ] );
            std::pair< int, int > indicesInEstimatedParametersSet = parametersToEstimate_->getIndicesForParameterDescription( estimatedParametersDescriptions[ i ] );
            indicesAndSizeEstimatedParameters_.push_back( std::make_pair( std::make_pair( indicesInEstimatedParametersSet.first, indicesInFullParametersSet.first ),
                                                                          indicesInFullParametersSet.second ) );
        }
    }


    std::pair< std::pair< Eigen::MatrixXd, Eigen::MatrixXd >, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > performPreEstimationSteps(
            std::shared_ptr< CovarianceAnalysisInput< ObservationScalarType, TimeType > > estimationInput,
            ParameterVectorType& newParameterEstimate,
            const bool calculateResiduals,
            const int numberOfIterations,
            bool& exceptionDuringPropagation,
            std::shared_ptr< propagators::SimulationResults< ObservationScalarType, TimeType > >& simulationResults )
    {
        // Get number of observations
        int totalNumberOfObservations = estimationInput->getObservationCollection( )->getTotalObservableSize( );

        // Re-integrate equations of motion and variational equations with new parameter estimate.
        try
        {
            if( ( numberOfIterations > 0 ) || ( estimationInput->getReintegrateEquationsOnFirstIteration( ) ) )
            {
                resetParameterEstimate( newParameterEstimate, estimationInput->getReintegrateVariationalEquations( ) );
            }

            if( std::dynamic_pointer_cast< EstimationInput< ObservationScalarType, TimeType > >( estimationInput ) != nullptr )
            {
                if ( std::dynamic_pointer_cast< EstimationInput< ObservationScalarType, TimeType > >( estimationInput )->getSaveStateHistoryForEachIteration( ) )
                {
                    simulationResults = variationalEquationsSolver_->getVariationalPropagationResults( );
                }
            }
        }
        catch( std::runtime_error& error )
        {
            std::cerr<<"Error when resetting parameters during parameter estimation: "<<std::endl<<
                     error.what( )<<std::endl<<"Terminating estimation"<<std::endl;
            exceptionDuringPropagation = true;
        }

        if( estimationInput->getPrintOutput( ) )
        {
            std::cout << "Calculating residuals and partials " << totalNumberOfObservations << std::endl;
        }

        // Calculate residuals and observation matrix for current parameter estimate.
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residuals;
        Eigen::MatrixXd designMatrix;
        if ( calculateResiduals )
        {
            calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >(
                    estimationInput->getObservationCollection( ), observationManagers_, totalNumberParameters_, totalNumberOfObservations, designMatrix, residuals, true );
        }
        else
        {
            calculateDesignMatrix< ObservationScalarType, TimeType >(
                    estimationInput->getObservationCollection( ), observationManagers_, totalNumberParameters_, totalNumberOfObservations, designMatrix );
        }

        // Divide partials matrix between estimated and consider parameters
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > designMatrices = separateEstimatedAndConsiderDesignMatrices( designMatrix, totalNumberOfObservations );

        return std::make_pair( designMatrices, residuals );
    }

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > separateEstimatedAndConsiderDesignMatrices(
            const Eigen::MatrixXd& designMatrix,
            const int numberObservations )
    {
        Eigen::MatrixXd designMatrixEstimatedParameters = Eigen::MatrixXd::Zero( numberObservations, numberEstimatedParameters_ );
        for ( unsigned int i = 0 ; i < indicesAndSizeEstimatedParameters_.size( ) ; i++ )
        {
            designMatrixEstimatedParameters.block( 0, indicesAndSizeEstimatedParameters_[ i ].first.first, numberObservations,
                                                   indicesAndSizeEstimatedParameters_[ i ].second )
                    = designMatrix.block( 0, indicesAndSizeEstimatedParameters_[ i ].first.second, numberObservations, indicesAndSizeEstimatedParameters_[ i ].second );
        }
        Eigen::MatrixXd designMatrixConsiderParameters = Eigen::MatrixXd::Zero( numberObservations, numberConsiderParameters_ );
        for ( unsigned int i = 0 ; i < indicesAndSizeConsiderParameters_.size( ) ; i++ )
        {
            designMatrixConsiderParameters.block( 0, indicesAndSizeConsiderParameters_[ i ].first.first, numberObservations,
                                                  indicesAndSizeConsiderParameters_[ i ].second )
                    = designMatrix.block( 0, indicesAndSizeConsiderParameters_[ i ].first.second, numberObservations, indicesAndSizeConsiderParameters_[ i ].second );
        }
        return std::make_pair( designMatrixEstimatedParameters, designMatrixConsiderParameters );
    }


    //! Boolean to denote whether any dynamical parameters are estimated
    bool integrateAndEstimateOrbit_;

    //! Object used to propagate/process the numerical solution of the variational equations/dynamics
    std::shared_ptr< propagators::VariationalEquationsSolver< ObservationScalarType, TimeType > >
    variationalEquationsSolver_;

    //! List of object that compute the values/partials of the observables
    std::map< observation_models::ObservableType,
    std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > > observationManagers_;

    //! Container object for all parameters that are to be estimated
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate_;

    //! Container object for consider parameters (if any)
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > considerParameters_;

    //! Container object for estimated and consider parameters (combined)
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > fullParameters_;

    SystemOfBodies bodies_;

    //! Current values of the vector of estimated parameters
    ParameterVectorType currentParameterEstimate_;

    //! Current values of the full vector of estimated and consider parameters
    ParameterVectorType currentFullParameterValues_;

    //! Consider parameters values
    ParameterVectorType considerParametersValues_;

    //! Total number of parameters (estimated and consider parameters together)
    unsigned int totalNumberParameters_;

    //! Number of estimated parameters
    unsigned int numberEstimatedParameters_;

    //! Number of consider parameters
    unsigned int numberConsiderParameters_;

    //std::vector< int > observationLinkParameterIndices_;

    //! Object used to interpolate the numerically integrated result of the state transition/sensitivity matrices.
    std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface >
    stateTransitionAndSensitivityMatrixInterface_;

    //! Object used to interpolate the numerically integrated result of the dependent variables.
    std::shared_ptr< propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface_;

    //! Container object for indices and sizes of consider parameters in the full estimated parameters set.
    std::vector< std::pair< std::pair< int, int >, int > > indicesAndSizeConsiderParameters_;

    //! Container object for indices and sizes of estimated parameters in the full estimated parameters set.
    std::vector< std::pair< std::pair< int, int >, int > > indicesAndSizeEstimatedParameters_;

    //! Boolean denoting whether consider parameters are included in the orbit determination
    bool considerParametersIncluded_;

};


//extern template class OrbitDeterminationManager< double, double >;



}

}
#endif // TUDAT_ORBITDETERMINATIONMANAGER_H
