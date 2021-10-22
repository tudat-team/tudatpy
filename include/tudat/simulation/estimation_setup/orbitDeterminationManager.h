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

#include <boost/make_shared.hpp>

#include "tudat/io/basicInputOutput.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/astro/observation_models/observationManager.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/estimation_setup/createObservationManager.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"

namespace tudat
{

namespace simulation_setup
{

////! Function to create a single vector of observation weights from weights sorted by link ends/observation type
///*!
// *  Function to create a single vector of observation weights from weights sorted by link ends/observation type,
// *  the ruesulting vector is sorted according to the iteration order of the inner and outer maps in the weightsData
// *  input
// *  \param weightsData Weights sorted by link ends and observation type
// *  \return Concatenated vector of weights
// */
//template< typename ObservationScalarType = double >
//Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedWeightsVector(
//        const typename std::map< observation_models::ObservableType, std::map<
//        observation_models::LinkEnds, Eigen::VectorXd > >& weightsData )
//{
//    typedef std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > >
//            WeightsDataStructure;

//    // Get total required size of weights vector
//    int totalNumberOfObservations = 0;
//    for( typename WeightsDataStructure::const_iterator observablesIterator =
//         weightsData.begin( ); observablesIterator != weightsData.end( ); observablesIterator++ )
//    {
//        for( typename std::map< observation_models::LinkEnds, Eigen::VectorXd >::const_iterator dataIterator =
//             observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
//        {
//            totalNumberOfObservations += dataIterator->second.rows( );
//        }
//    }
//    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > concatenatedWeights = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalNumberOfObservations, 1 );

//    // Iterate over all observations and concatenate the weight vectors.
//    int currentIndex = 0;
//    for( typename WeightsDataStructure::const_iterator observablesIterator =
//         weightsData.begin( ); observablesIterator != weightsData.end( ); observablesIterator++ )
//    {
//        for( typename std::map< observation_models::LinkEnds, Eigen::VectorXd >::const_iterator dataIterator =
//             observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
//        {
//            concatenatedWeights.segment( currentIndex, dataIterator->second.rows( ) ) =
//                    dataIterator->second.template cast< ObservationScalarType >( );
//            currentIndex += dataIterator->second.rows( );
//        }
//    }

//    return concatenatedWeights;
//}

//! Top-level class for performing orbit determination.
/*!
 *  Top-level class for performing orbit determination. All required propagation/estimation settings are provided to
 *  this class, which then creates all objects needed for the propagation and estimation process. The parameter
 *  estimation itself is performed by providing measurement data and related metadata (as PodInput) to the estimateParameters
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
//    std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > SingleObservablePodInputType;

//    //! Typedef for complete set of observations data, as used in orbit determination
//    typedef std::map< observation_models::ObservableType, SingleObservablePodInputType > PodInputType;

//    //! Typedef for complete set of observations data in alternative form, as used in orbit determination, convertible to PodInputType
//    //! by convertPodInput function.
//    typedef std::map< observation_models::ObservableType,
//    std::map< observation_models::LinkEnds, std::pair< std::map< TimeType, ObservationScalarType >,
//    observation_models::LinkEndType > > > AlternativePodInputType;

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
        parametersToEstimate_( parametersToEstimate )
    {
        initializeOrbitDeterminationManager( bodies, observationSettingsList, { integratorSettings }, propagatorSettings,
                                             propagateOnCreation );
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
        initializeOrbitDeterminationManager( bodies, observationSettingsList, integratorSettings, propagatorSettings,
                                             propagateOnCreation );
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

//    //! Function to determine the number of observations per link end.
//    /*!
//     *  Function to determine the number of observations per link end from a map of observations for each link ends.
//     *  The input type is directly related to the data stored for a single observable in PodInput::PodInputDataType.
//     *  \param dataPerLinkEnd Map of observations and times for a set of link ends.
//     *  \return Vector of size of number of observations in input map (in order of forward iterator over input map).
//     */
//    static std::vector< int > getNumberOfObservationsPerLinkEnd(
//            const SingleObservablePodInputType& dataPerLinkEnd )
//    {
//        // Declare output vector.
//        std::vector< int > numberOfObservations;

//        // Iterate over all link ends.
//        for( typename SingleObservablePodInputType::const_iterator dataIterator =
//             dataPerLinkEnd.begin( ); dataIterator != dataPerLinkEnd.end( ); dataIterator++  )
//        {
//            // Add number of observations for current link ends.
//            numberOfObservations.push_back( dataIterator->second.first.rows( ) );
//        }

//        return numberOfObservations;
//    }

//    //! Function to determine total number of observation and number of observations per observable
//    /*!
//     *  Function to determine total number of observation and number of observations per observable from the complete set
//     *  of measurement data.
//     *  \param observationsAndTimes Set of measurement data per obsevable type and link ends
//     *  \return Pair first: map with number of observations per observable type, second: total number of observations
//     *  (i.e. sum of valus of first)
//     */
//    static std::pair< std::map< observation_models::ObservableType, int >, int > getNumberOfObservationsPerObservable(
//            const PodInputType& observationsAndTimes )
//    {
//        // Initialize counters.
//        std::map< observation_models::ObservableType, int > numberOfObservations;
//        int totalNumberOfObservations = 0;

//        // Iterate over all observabel types.
//        for( typename PodInputType::const_iterator observablesIterator = observationsAndTimes.begin( );
//             observablesIterator != observationsAndTimes.end( ); observablesIterator++ )
//        {
//            // Initialize number of observations for current observable
//            numberOfObservations[ observablesIterator->first ] = 0;

//            // Iterate over all link ends.
//            for( typename SingleObservablePodInputType::const_iterator dataIterator = observablesIterator->second.begin( );
//                 dataIterator != observablesIterator->second.end( ); dataIterator++  )
//            {
//                // Add number of observations with given link ends.
//                numberOfObservations[ observablesIterator->first ] += dataIterator->second.first.size( );
//            }

//            // Add to total number of observations.
//            totalNumberOfObservations += numberOfObservations[ observablesIterator->first ];
//        }

//        return std::make_pair( numberOfObservations, totalNumberOfObservations );
//    }

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
    void calculateObservationMatrixAndResiduals(
            const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationsCollection,
            const int parameterVectorSize, const int totalObservationSize,
            std::pair< Eigen::VectorXd, Eigen::MatrixXd >& residualsAndPartials  )
    {
        // Initialize return data.
        residualsAndPartials.second = Eigen::MatrixXd::Zero( totalObservationSize, parameterVectorSize );
        residualsAndPartials.first = Eigen::VectorXd::Zero( totalObservationSize );

        typename observation_models::ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets
                sortedObservations = observationsCollection->getObservations( );

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
                    std::pair< ObservationVectorType, Eigen::MatrixXd > observationsWithPartials;
                    observationsWithPartials = observationManagers_[ currentObservableType ]->
                            computeObservationsWithPartials(
                                currentObservations->getObservationTimes( ), currentLinkEnds,
                                currentObservations->getReferenceLinkEnd( ) );

                    // Compute residuals for current link ends and observabel type.
                    residualsAndPartials.first.segment( observationIndices.first, observationIndices.second ) =
                            ( currentObservations->getObservationsVector( ) - observationsWithPartials.first ).template cast< double >( );

                    // Set current observation partials in matrix of all partials
                    residualsAndPartials.second.block( observationIndices.first, 0, observationIndices.second, parameterVectorSize ) =
                            observationsWithPartials.second;

                }

            }

            std::pair< int, int > observableStartAndSize = observationsCollection->getObservationTypeStartAndSize( ).at( currentObservableType );

            observation_models::checkObservationResidualDiscontinuities(
                        residualsAndPartials.first.block( observableStartAndSize.first, 0, observableStartAndSize.second, 1 ),
                        currentObservableType );
        }

    }



    //! Function to normalize the matrix of partial derivatives so that each column is in the range [-1,1]
    /*!
     * Function to normalize the matrix of partial derivatives so that each column is in the range [-1,1]
     * \param observationMatrix Matrix of partial derivatives. Matrix modified by this function, and normalized matrix is
     * returned by reference
     * \return Vector with scaling values used for normalization
     */
    Eigen::VectorXd normalizeObservationMatrix( Eigen::MatrixXd& observationMatrix )
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

    //! Function to perform parameter estimation from measurement data.
    /*!
     *  Function to perform parameter estimation, including orbit determination, i.e. body initial states, from measurement data.
     *  All observable types and link ends per obsevable types that are included in the measurement data input must have been
     *  provided to the constructor by the observationSettingsMap parameter.
     *  \param podInput Object containing all measurement data, associated metadata, including measurement weight, and a priori
     *  estimate for covariance matrix and parameter adjustment.
     *  \param convergenceChecker Object used to check convergence/termination of algorithm
     *  \return Object containing estimated parameter value and associateed data, such as residuals and observation partials.
     */
    std::shared_ptr< PodOutput< ObservationScalarType, TimeType > > estimateParameters(
            const std::shared_ptr< PodInput< ObservationScalarType, TimeType > > podInput,
            std::shared_ptr< EstimationConvergenceChecker > convergenceChecker =
            std::make_shared< EstimationConvergenceChecker >( ) )

    {
        currentParameterEstimate_ = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );

        // Get size of parameter vector and number of observations (total and per type)
        int parameterVectorSize = currentParameterEstimate_.size( );
        int totalNumberOfObservations = podInput->getObservationCollection( )->getTotalObservableSize( );

        // Declare variables to be returned (i.e. results from best iteration)
        double bestResidual = TUDAT_NAN;
        ParameterVectorType bestParameterEstimate = ParameterVectorType::Constant( parameterVectorSize, TUDAT_NAN );
        Eigen::VectorXd bestTransformationData = Eigen::VectorXd::Constant( parameterVectorSize, TUDAT_NAN );
        Eigen::VectorXd bestResiduals = Eigen::VectorXd::Constant( totalNumberOfObservations, TUDAT_NAN );
        Eigen::MatrixXd bestInformationMatrix = Eigen::MatrixXd::Constant( totalNumberOfObservations, parameterVectorSize, TUDAT_NAN );
        Eigen::VectorXd bestWeightsMatrixDiagonal = Eigen::VectorXd::Constant( totalNumberOfObservations, TUDAT_NAN );
        Eigen::MatrixXd bestInverseNormalizedCovarianceMatrix = Eigen::MatrixXd::Constant( parameterVectorSize, parameterVectorSize, TUDAT_NAN );

        std::vector< Eigen::VectorXd > residualHistory;
        std::vector< Eigen::VectorXd > parameterHistory;
        std::vector< std::vector< std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > > dynamicsHistoryPerIteration;
        std::vector< std::vector< std::map< TimeType, Eigen::VectorXd > > > dependentVariableHistoryPerIteration;

        // Declare residual bookkeeping variables
        std::vector< double > rmsResidualHistory;
        double residualRms;

        // Declare variables to be used in loop.

        // Set current parameter estimate as both previous and current estimate
        ParameterVectorType newParameterEstimate = currentParameterEstimate_ +
                podInput->getInitialParameterDeviationEstimate( );
        ParameterVectorType oldParameterEstimate = currentParameterEstimate_;

        int numberOfEstimatedParameters = parameterVectorSize;

        bool exceptionDuringPropagation = false, exceptionDuringInversion = false;
        // Iterate until convergence (at least once)
        int numberOfIterations = 0;
        do
        {
//            try
            {
                // Re-integrate equations of motion and variational equations with new parameter estimate.
                if( ( numberOfIterations > 0 ) ||( podInput->getReintegrateEquationsOnFirstIteration( ) ) )
                {
                    resetParameterEstimate( newParameterEstimate, podInput->getReintegrateVariationalEquations( ) );
                }

/*
                if( podInput->getSaveStateHistoryForEachIteration( ) )
                {
                    dynamicsHistoryPerIteration.push_back(
                                variationalEquationsSolver_->getDynamicsSimulatorBase( )->getEquationsOfMotionNumericalSolutionBase( ) );
                    dependentVariableHistoryPerIteration.push_back(
                                variationalEquationsSolver_->getDynamicsSimulatorBase( )->getDependentVariableNumericalSolutionBase( ) );
                }*/

                // Michael
                if( podInput->getSaveStateHistoryForEachIteration( ) )
                {
                    if( std::dynamic_pointer_cast< propagators::HybridArcVariationalEquationsSolver< ObservationScalarType, TimeType > >(  variationalEquationsSolver_) != nullptr )
                    {

                        std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< ObservationScalarType, TimeType > > hybridArcSolver =
                                std::dynamic_pointer_cast< propagators::HybridArcVariationalEquationsSolver< ObservationScalarType, TimeType > >(  variationalEquationsSolver_);

                    dynamicsHistoryPerIteration.push_back(
                            hybridArcSolver->getMultiArcSolver( )->getDynamicsSimulatorBase( )->getEquationsOfMotionNumericalSolutionBase( ));
                    dependentVariableHistoryPerIteration.push_back(
                        hybridArcSolver->getMultiArcSolver( )->getDynamicsSimulatorBase( )->getDependentVariableNumericalSolutionBase( ));

                       /* dynamicsHistoryPerIteration.push_back(
                                hybridArcSolver->getSingleArcSolver( )->getDynamicsSimulatorBase( )->getEquationsOfMotionNumericalSolutionBase( ));
                        dependentVariableHistoryPerIteration.push_back(
                                hybridArcSolver->getSingleArcSolver( )->getDynamicsSimulatorBase( )->getDependentVariableNumericalSolutionBase( ));
*/
                    }
                }

            }
//            catch( std::runtime_error& error )
//            {
//                std::cerr<<"Error when resetting parameters during parameter estimation: "<<std::endl<<
//                           error.what( )<<std::endl<<"Terminating estimation"<<std::endl;
//                exceptionDuringPropagation = true;
//                break;
//            }

            oldParameterEstimate = newParameterEstimate;

            if( podInput->getPrintOutput( ) )
            {
                std::cout << "Calculating residuals and partials " << totalNumberOfObservations << std::endl;
            }
            // Calculate residuals and observation matrix for current parameter estimate.
            std::pair< Eigen::VectorXd, Eigen::MatrixXd > residualsAndPartials;
            calculateObservationMatrixAndResiduals(
                        podInput->getObservationCollection( ), parameterVectorSize, totalNumberOfObservations, residualsAndPartials );

            Eigen::VectorXd transformationData = normalizeObservationMatrix( residualsAndPartials.second );

            Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix = Eigen::MatrixXd::Zero(
                        numberOfEstimatedParameters, numberOfEstimatedParameters );

            Eigen::MatrixXd inverseAPrioriCovariance = podInput->getInverseOfAprioriCovariance( );
            for( int j = 0; j < numberOfEstimatedParameters; j++ )
            {
                for( int k = 0; k < numberOfEstimatedParameters; k++ )
                {
                    normalizedInverseAprioriCovarianceMatrix( j, k ) = inverseAPrioriCovariance( j, k ) /
                            ( transformationData( j ) * transformationData( k ) );
                }
            }

            // Perform least squares calculation for correction to parameter vector.
            std::pair< Eigen::VectorXd, Eigen::MatrixXd > leastSquaresOutput;
            try
            {
                Eigen::MatrixXd constraintStateMultiplier;
                Eigen::VectorXd constraintRightHandSide;
                parametersToEstimate_->getConstraints( constraintStateMultiplier, constraintRightHandSide );
                leastSquaresOutput =
                        std::move( linear_algebra::performLeastSquaresAdjustmentFromInformationMatrix(
                                       residualsAndPartials.second.block( 0, 0, residualsAndPartials.second.rows( ), numberOfEstimatedParameters ),
                                       residualsAndPartials.first, podInput->getWeightsMatrixDiagonals( ),
                                       normalizedInverseAprioriCovarianceMatrix, 1, 1.0E8, constraintStateMultiplier, constraintRightHandSide ) );

                if( constraintStateMultiplier.rows( ) > 0 )
                {
                    leastSquaresOutput.first.conservativeResize( parameterVectorSize );
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
                    ( leastSquaresOutput.first.cwiseQuotient( transformationData.segment( 0, numberOfEstimatedParameters ) ) ).
                    template cast< ObservationScalarType >( );

            //            std::cout<<"LSQ: "<<leastSquaresOutput.first<<std::endl<<
            //                       transformationData.segment( 0, numberOfEstimatedParameters ).transpose( )<<std::endl;



            // Update value of parameter vector
            newParameterEstimate = oldParameterEstimate + parameterAddition;
            parametersToEstimate_->template resetParameterValues< ObservationScalarType >( newParameterEstimate );
            newParameterEstimate = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );

            if( podInput->getSaveResidualsAndParametersFromEachIteration( ) )
            {
                residualHistory.push_back( residualsAndPartials.first );
                if( numberOfIterations == 0 )
                {
                    parameterHistory.push_back( oldParameterEstimate.template cast< double >( ) );
                }
                parameterHistory.push_back( newParameterEstimate.template cast< double >( ) );
            }

            oldParameterEstimate = newParameterEstimate;

            if( podInput->getPrintOutput( ) )
            {
                std::cout << "Parameter update" << parameterAddition.transpose( ) << std::endl;
            }

            // Calculate mean residual for current iteration.
            residualRms = linear_algebra::getVectorEntryRootMeanSquare( residualsAndPartials.first );

            rmsResidualHistory.push_back( residualRms );
            if( podInput->getPrintOutput( ) )
            {
                std::cout << "Current residual: " << residualRms << std::endl;
            }

            // If current iteration is better than previous one, update 'best' data.
            if( residualRms < bestResidual || !( bestResidual == bestResidual ) )
            {
                bestResidual = residualRms;
                bestParameterEstimate = std::move( oldParameterEstimate );
                bestResiduals = std::move( residualsAndPartials.first );
                if( podInput->getSaveInformationMatrix( ) )
                {
                    bestInformationMatrix = std::move( residualsAndPartials.second );
                }
                bestWeightsMatrixDiagonal = std::move( podInput->getWeightsMatrixDiagonals( ) );
                bestTransformationData = std::move( transformationData );
                bestInverseNormalizedCovarianceMatrix = std::move( leastSquaresOutput.second );
            }


            // Increment number of iterations
            numberOfIterations++;

            // Check for convergence
        } while( convergenceChecker->isEstimationConverged( numberOfIterations, rmsResidualHistory ) == false );

        if( podInput->getPrintOutput( ) )
        {
            std::cout << "Final residual: " << bestResidual << std::endl;
        }


        std::shared_ptr< PodOutput< ObservationScalarType, TimeType > > podOutput =
                std::make_shared< PodOutput< ObservationScalarType, TimeType > >(
                    bestParameterEstimate, bestResiduals, bestInformationMatrix, bestWeightsMatrixDiagonal, bestTransformationData,
                    bestInverseNormalizedCovarianceMatrix, bestResidual, residualHistory, parameterHistory, exceptionDuringInversion,
                    exceptionDuringPropagation );

        if( podInput->getSaveStateHistoryForEachIteration( ) )
        {
            podOutput->setStateHistories(
                        dynamicsHistoryPerIteration, dependentVariableHistoryPerIteration );
        }

        return podOutput;
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
//     *  Function to convert from one representation of all measurement data (AlternativePodInputType) to the other (PodInputType).
//     *  In the former, the vector of times and vector of associated observations are stored separately.  In the latter, they are
//     *  stored as a map with time as key and obsevation as value.
//     *  \param alternativePodInput Original representation of measurement data
//     *  \return Converted measurement data.
//     */
//    static PodInputType convertPodInput( const AlternativePodInputType& alternativePodInput )
//    {
//        // Declare return data structure.
//        PodInputType convertedObservationAndTimes;

//        // Iterate over all observable types.
//        for( typename AlternativePodInputType::const_iterator inputIterator = alternativePodInput.begin( );
//             inputIterator != alternativePodInput.end( ); inputIterator++ )
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
            const std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettings,
            const std::shared_ptr< propagators::PropagatorSettings< ObservationScalarType > > propagatorSettings,
            const bool propagateOnCreation = true )
    {
        using namespace numerical_integrators;
        using namespace orbit_determination;
        using namespace observation_models;

        // Check if any dynamics is to be estimated
        std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > >
                initialDynamicalStates =
                estimatable_parameters::getListOfInitialDynamicalStateParametersEstimate< ObservationScalarType >(
                    parametersToEstimate_ );
        if( initialDynamicalStates.size( ) > 0 )
        {
            integrateAndEstimateOrbit_ = true;
        }
        else
        {
            integrateAndEstimateOrbit_ = false;
        }

        if( integrateAndEstimateOrbit_ )
        {
            variationalEquationsSolver_ =
                    simulation_setup::createVariationalEquationsSolver(
                        bodies, integratorSettings, propagatorSettings, parametersToEstimate_, 1,
                        std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), 0, propagateOnCreation );
        }

        if( integrateAndEstimateOrbit_ )
        {
            stateTransitionAndSensitivityMatrixInterface_ = variationalEquationsSolver_->getStateTransitionMatrixInterface( );
        }
        else if( propagatorSettings == nullptr )
        {
            stateTransitionAndSensitivityMatrixInterface_ = createStateTransitionAndSensitivityMatrixInterface(
                        propagatorSettings, 0, parametersToEstimate_->getParameterSetSize( ) );
        }
        else
        {
            throw std::runtime_error( "Error, cannot parse propagator settings without estimating dynamics in OrbitDeterminationManager" );
        }

        // Iterate over all observables and create observation managers.
        std::map< ObservableType, std::vector< std::shared_ptr< ObservationModelSettings > > > sortedObservationSettingsList =
                sortObservationModelSettingsByType( observationSettingsList );
        for( auto it : sortedObservationSettingsList )
        {
            // Call createObservationSimulator of required observation size
            ObservableType observableType = it.first;

            // Create observation manager for current observable.
            observationManagers_[ observableType ] =
                    createObservationManagerBase< ObservationScalarType, TimeType >(
                        observableType,
                        it.second,
                        bodies, parametersToEstimate_,
                        stateTransitionAndSensitivityMatrixInterface_ );
        }

        // Set current parameter estimate from body initial states and parameter set.
        currentParameterEstimate_ = parametersToEstimate_->template getFullParameterValues< ObservationScalarType >( );

        //        std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParameters =
        //                parametersToEstimate_->getDoubleParameters( );
        //        for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator
        //             parameterIterator = doubleParameters.begin( ); parameterIterator != doubleParameters.end( ); parameterIterator++ )
        //        {
        //            if( estimatable_parameters::isParameterObservationLinkProperty(
        //                        parameterIterator->second->getParameterName( ).first ) )
        //            {
        //                observationLinkParameterIndices_.push_back( parameterIterator->first );
        //            }
        //        }

        //        std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParameters =
        //                parametersToEstimate_->getVectorParameters( );
        //        for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >::iterator
        //             parameterIterator = vectorParameters.begin( ); parameterIterator != vectorParameters.end( ); parameterIterator++ )
        //        {
        //            if( estimatable_parameters::isParameterObservationLinkProperty(
        //                        parameterIterator->second->getParameterName( ).first ) )
        //            {
        //                for( int i = 0; i < parameterIterator->second->getParameterSize( ); i++ )
        //                {
        //                    observationLinkParameterIndices_.push_back( parameterIterator->first + i );
        //                }
        //            }
        //        }

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

    //! Current values of the vector of estimated parameters
    ParameterVectorType currentParameterEstimate_;

    //std::vector< int > observationLinkParameterIndices_;

    //! Object used to interpolate the numerically integrated result of the state transition/sensitivity matrices.
    std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface >
    stateTransitionAndSensitivityMatrixInterface_;

};

extern template class OrbitDeterminationManager< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class OrbitDeterminationManager< double, Time >;
extern template class OrbitDeterminationManager< long double, double >;
extern template class OrbitDeterminationManager< long double, Time >;
#endif



}

}
#endif // TUDAT_ORBITDETERMINATIONMANAGER_H
