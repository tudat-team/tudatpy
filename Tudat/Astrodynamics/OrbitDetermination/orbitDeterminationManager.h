#ifndef ORBITDETERMINATIONMANAGER_H
#define ORBITDETERMINATIONMANAGER_H

#include <algorithm>

#include <boost/make_shared.hpp>

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Mathematics/Interpolators/lookupScheme.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/ObservationModels/observationManager.h"
#include "Tudat/Astrodynamics/OrbitDetermination/podInputOutputTypes.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrectionPartials.h"
#include "Tudat/SimulationSetup/EstimationSetup/createAccelerationPartials.h"
#include "Tudat/SimulationSetup/PropagationSetup/variationalEquationsSolver.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationManager.h"

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType = double, typename TimeType = double >
Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getConcatenatedWeightsVector(
        const typename std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > >& weightsData )
{
    typedef std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > > WeightsDataStructure;

    int totalNumberOfObservations = 0;
    for( typename WeightsDataStructure::const_iterator observablesIterator =
         weightsData.begin( ); observablesIterator != weightsData.end( ); observablesIterator++ )
    {
        for( typename std::map< observation_models::LinkEnds, Eigen::VectorXd >::const_iterator dataIterator =
             observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
        {
            totalNumberOfObservations += dataIterator->second.rows( );
        }
    }

    // Iterate over all observations and concatenate the time vectors.
    Eigen::VectorXd concatenatedWeights = Eigen::VectorXd::Zero( totalNumberOfObservations, 1 );

    int currentIndex = 0;
    for( typename WeightsDataStructure::const_iterator observablesIterator =
         weightsData.begin( ); observablesIterator != weightsData.end( ); observablesIterator++ )
    {
        for( typename std::map< observation_models::LinkEnds, Eigen::VectorXd >::const_iterator dataIterator =
             observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
        {
            concatenatedWeights.segment( currentIndex, dataIterator->second.rows( ) ) = dataIterator->second;
            currentIndex += dataIterator->second.rows( );
        }
    }

    return concatenatedWeights;
}


class EstimationConvergenceChecker
{
public:
    EstimationConvergenceChecker(
            const unsigned int maximumNumberOfIterations = 5,
            const double minimumResidualChange = 0.0,
            const double minimumResidual = 1.0E-20,
            const int numberOfIterationsWithoutImprovement = 2 ):
        maximumNumberOfIterations_( maximumNumberOfIterations ), minimumResidualChange_( minimumResidualChange ),
        minimumResidual_( minimumResidual ),
        numberOfIterationsWithoutImprovement_( numberOfIterationsWithoutImprovement ),
        previousResidual_( 1.0E300 )
    {

    }

    bool isEstimationConverged( const int numberOfIterations, const std::vector< double > rmsResidualHistory )
    {
        bool isConverged = 0;
        if( numberOfIterations >= maximumNumberOfIterations_ )
        {
            std::cout<<"Maximum number of iterations reached"<<std::endl;
            isConverged = 1;
        }
        if( rmsResidualHistory[ rmsResidualHistory.size( ) - 1 ] < minimumResidual_ )
        {
            std::cout<<"Required residual level achieved"<<std::endl;
            isConverged = 1;
        }
        if( ( findIndexOfMaximumEntry( rmsResidualHistory ) - rmsResidualHistory.size( ) ) < numberOfIterationsWithoutImprovement_ )
        {
            std::cout<<"Too many iterations without parameter improvement"<<std::endl;
            isConverged = 1;
        }
        if( rmsResidualHistory.size( ) > 1 )
        {
            if( std::fabs( rmsResidualHistory.at( rmsResidualHistory.size( )  - 1 ) -
                    rmsResidualHistory.at( rmsResidualHistory.size( )  - 2 ) ) < minimumResidualChange_ )
            {
                isConverged = 1;
            }
        }
        return isConverged;
    }
protected:

    unsigned int findIndexOfMaximumEntry( std::vector< double > rmsResidualHistory )
    {
        int maximumEntry = 0;
        double maximumValue = rmsResidualHistory[ 0 ];
        for( unsigned int i = 0; i < rmsResidualHistory.size( ); i++ )
        {
            if( rmsResidualHistory[ i ] >= maximumValue )
            {
                maximumEntry = i;
                maximumValue = rmsResidualHistory[ i ];
            }
        }
        return maximumEntry;
    }

    int maximumNumberOfIterations_;
    double minimumResidualChange_;

    double minimumResidual_;
    unsigned int numberOfIterationsWithoutImprovement_;

    double previousResidual_;
};

//! Top-level class for performing orbit determination.
/*!
 *  Top-level class for performing orbit determination. Acceleration models, parameters to estimate, observable types and link ends
 *  as well as various settings for observation and dynamical simulation are provided to the class (see constructor). The parameter
 *  estimation process is then performed by providing measurement data and related metadata (as PodInput) to the estimateParameters function.
 *  \tparam ObservationScalarType Scalar Data type for observations (double, long double)
 *  \tparam ObservationScalarType TimeType Data type for time representation (double, tudat::Time)
 *  \tparam StateScalarType Data type for time entries of state vector of integrated bodies (double, tudat::Time)
 *  (double, long double)
 */
template< typename ObservationScalarType = double, typename TimeType = double, typename StateScalarType = ObservationScalarType >
class OrbitDeterminationManager
{
public:

    //! Typedef for vector of observations.
    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;

    //! Typedef for vector of parameters.
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > ParameterVectorType;


    typedef std::map< observation_models::LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > SingleObservablePodInputType;

    //! Typedef for complete set of observations data, as used in orbit determination
    typedef std::map< observation_models::ObservableType, SingleObservablePodInputType > PodInputType;

    //! Typedef for complete set of observations data in alternative form, as used in orbit determination, convertible to PodInputType
    //! by convertPodInput function.
    typedef std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, std::pair< std::map< TimeType, ObservationScalarType >, observation_models::LinkEndType > > > AlternativePodInputType;

    //! Constructor
    /*!
     *
     *  \param bodyMap Map of body objects with names of bodies, storing all environment models used in simulation.
     *  \param parametersToEstimate Container object for all parameters that are to be estimated (excluding initial states of bodies)
     *  \param bodiesToIntegrate List of bodies that are to be numerically integrated in simulation
     *  \param bodiesToEstimate List of bodies for which initial states are to be estimated (note that it is advised that this
     *  vector be equal to the bodiesToIntegrate, to prevent effects of parameter adjustments not covered by the partial derivatives)
     *  \param linkEndsPerObservable Sets of link ends (i.e. transmitter, receiver, etc.) per observable type for which measurement data is
     *  to be provided in orbit determination process (through estimateParameters function)
     *  \param initialBodyStates Initial states of bodies to integrate, in frames centered on their respective central bodies
     *  (defined in in propagatorSettings type)
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param observableCorrections List of correction function to for each observable type and set of link ends (relativistic, troposphere, etc. )
     *  \param observationViabilitySettings List of observation viability functions to for each observable type and set of link ends
     *  (minimum elevation angle, body occultation, etc. )
     */
    OrbitDeterminationManager(
            const NamedBodyMap &bodyMap,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const std::map< observation_models::ObservableType, std::vector< observation_models::LinkEnds > >& linkEndsPerObservable,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::map< observation_models::ObservableType, observation_models::LightTimeCorrectionSettingsMap >& observableCorrections =
            ( std::map< observation_models::ObservableType, observation_models::LightTimeCorrectionSettingsMap >( ) ),
            const int integrateOnCreation = 1 ):
        parametersToEstimate_( parametersToEstimate ), linkEndsPerObservable_( linkEndsPerObservable )
    {
        using namespace numerical_integrators;
        using namespace orbit_determination;
        using namespace observation_models;

        std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > > initialDynamicalStates =
                estimatable_parameters::getListOfInitialDynamicalStateParametersEstimate< StateScalarType >( parametersToEstimate );

        if( initialDynamicalStates.size( ) > 0 )
        {
            if( integrateOnCreation )
            {
                integrateAndEstimateOrbit_ = true;

                    variationalEquationsSolver_ = boost::make_shared< propagators::SingleArcVariationalEquationsSolver
                            < StateScalarType, TimeType, StateScalarType > >(
                                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, 1,
                                boost::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), 0 );
            }
            else
            {
                integrateAndEstimateOrbit_ = true;

                    variationalEquationsSolver_ = boost::make_shared< propagators::SingleArcVariationalEquationsSolver
                            < StateScalarType, TimeType, StateScalarType > >(
                                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, 1,
                                boost::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), 1, 0 );
            }
        }
        else
        {
            integrateAndEstimateOrbit_ = false;
        }

        std::cout<<"Dynamics simulator set and intialized "<<std::endl;

        if( integrateAndEstimateOrbit_ && integrateOnCreation )
        {
            stateTransitionAndSensitivityMatrixInterface_ =
                    variationalEquationsSolver_->getStateTransitionMatrixInterface( );
        }
        else if( propagatorSettings == NULL )
        {
            std::cout<<"Making st interface"<<std::endl;
            stateTransitionAndSensitivityMatrixInterface_ = boost::make_shared<
                    propagators::SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                        0, parametersToEstimate_->getParameterSetSize( ) );
        }
        else if( propagatorSettings != NULL )
        {
            std::cout<<"Making st interface"<<std::endl;
            stateTransitionAndSensitivityMatrixInterface_ = boost::make_shared<
                    propagators::SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                        parametersToEstimate_->getInitialDynamicalStateParameterSize( ), parametersToEstimate_->getParameterSetSize( ) );
        }

        // Iterate over all observables and create observation managers.
        for( std::map< observation_models::ObservableType, std::vector< observation_models::LinkEnds > >::const_iterator observablesIterator = linkEndsPerObservable.begin( );
             observablesIterator != linkEndsPerObservable.end( ); observablesIterator++ )
        {
            // Get observable corrections for current observable.
            std::map< observation_models::LinkEnds, std::vector< boost::shared_ptr< LightTimeCorrectionSettings > > > singleObservableCorrections;
            if( observableCorrections.count( observablesIterator->first ) > 0 )
            {
                singleObservableCorrections = observableCorrections.at( observablesIterator->first );
            }

            // Create observation manager for current observable.
            observationManagers_[ observablesIterator->first ] =
                    createObservationManagerBase< ObservationScalarType, TimeType, StateScalarType >(
                        observablesIterator->first, observablesIterator->second, bodyMap, parametersToEstimate,
                        stateTransitionAndSensitivityMatrixInterface_,
                        singleObservableCorrections );
        }

        std::cout<<"Observation managers created "<<std::endl;

        // Set current parameter estimate from body initial states and parameter set.
        currentParameterEstimate_ = parametersToEstimate_->template getFullParameterValues< StateScalarType >( );
    }

    //! Function to retrieve map of all observation managers
    /*!
     *  Function to retrieve map of all observation managers. A single observation manager can simulate observations and
     *  calculate observation partials for all link ends involved in the given observable type.
     *  \return Map of observation managers for all observable types involved in current orbit determination.
     */
    std::map< observation_models::ObservableType, boost::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType > > >
    getObservationManagers( ) const
    {
        return observationManagers_;
    }

    //! Function to determine the number of observations per link end.
    /*!
     *  Function to determine the number of observations per link end from a map of observations for each link ends. The input type is
     *  directly related to the data stored for a single observable in PodInput::PodInputDataType.
     *  \param Map of observations and times for a set of link ends.
     *  \return Vector of size of number of observations in input map (in order of forward iterator over input map).
     */
    static std::vector< int > getNumberOfObservationsPerLinkEnd(
            const SingleObservablePodInputType& dataPerLinkEnd )
    {
        // Declare output vector.
        std::vector< int > numberOfObservations;

        // Iterate over all link ends.
        for( typename SingleObservablePodInputType::const_iterator dataIterator =
             dataPerLinkEnd.begin( ); dataIterator != dataPerLinkEnd.end( ); dataIterator++  )
        {
            // Add number of observations for current link ends.
            numberOfObservations.push_back( dataIterator->second.first.rows( ) );
        }

        return numberOfObservations;
    }

    //! Function to determine total number of observation and number of observations per observable
    /*!
     *  Function to determine total number of observation and number of observations per observable from the complete set of measurement data.
     *  \param observationsAndTimes Set of measurement data per obsevable type and link ends
     *  \return Pair first: map with number of observations per observable type, second: total number of observations (i.e. sum of valus of first)
     */
    static std::pair< std::map< observation_models::ObservableType, int >, int > getNumberOfObservationsPerObservable(
            const PodInputType& observationsAndTimes )
    {
        //int currentObservationSize;

        // Initialize counters.
        std::map< observation_models::ObservableType, int > numberOfObservations;
        int totalNumberOfObservations = 0;

        // Iterate over all observabel types.
        for( typename PodInputType::const_iterator observablesIterator = observationsAndTimes.begin( );
             observablesIterator != observationsAndTimes.end( ); observablesIterator++ )
        {
            // Initialize number of observations for current observable
            numberOfObservations[ observablesIterator->first ] = 0;

            // Iterate over all link ends.
            for( typename SingleObservablePodInputType::const_iterator dataIterator = observablesIterator->second.begin( );
                 dataIterator != observablesIterator->second.end( ); dataIterator++  )
            {
                // Add number of observations with given link ends.
                numberOfObservations[ observablesIterator->first ] += dataIterator->second.first.size( );
            }

            // Add to total number of observations.
            totalNumberOfObservations += numberOfObservations[ observablesIterator->first ];
        }

        return std::make_pair( numberOfObservations, totalNumberOfObservations );
    }

    //! Function to calculate the observation partials matrix and residuals
    /*!
     *  This function calculates the observation partials matrix and residuals, based on the state transition matrix, sensitivity matrix and
     *  body states resulting from the previous numerical integration iteration. Partials and observations are calculated by the
     *  observationManagers_.
     *  \param observationsAndTime Observable values and associated time tags, per observable type and set of link ends.
     *  \param parameterVectorSize Length of the vector of estimated parameters
     *  \param totalNumberOfObservations Total number of observations in observationsAndTimes map.
     *  \return Pair of residuals of computed w.r.t. input observable values and partials of observables w.r.t. parameter vector.
     */
    void calculateObservationMatrixAndResiduals(
            const PodInputType& observationsAndTimes, const int parameterVectorSize, const int totalObservationSize,
            std::pair< Eigen::VectorXd, Eigen::MatrixXd >& residualsAndPartials  )
    {
        // Initialize return data.
        residualsAndPartials.second = Eigen::MatrixXd::Zero( totalObservationSize, parameterVectorSize );
        residualsAndPartials.first = Eigen::VectorXd::Zero( totalObservationSize );

        // Declare variable denoting current index in vector of all observations.
        int startIndex = 0;

        //boost::function< TimeType( const TimeType& ) > inputTimeConversionFunction;
        //bool isTimeConverterFound = 0;
        std::vector< TimeType > simulationInputTime;

        // Iterate over all observable types in observationsAndTimes
        for( typename PodInputType::const_iterator observablesIterator = observationsAndTimes.begin( );
             observablesIterator != observationsAndTimes.end( ); observablesIterator++ )
        {
            // Iterate over all link ends for current observable type in observationsAndTimes
            for( typename SingleObservablePodInputType::const_iterator dataIterator = observablesIterator->second.begin( );
                 dataIterator != observablesIterator->second.end( ); dataIterator++  )
            {

                simulationInputTime.clear( );

                    simulationInputTime = dataIterator->second.second.first;


                // Compute estimated ranges and range partials from current parameter estimate.
                std::pair< ObservationVectorType, Eigen::MatrixXd > observationsWithPartials;

                observationsWithPartials = observationManagers_[ observablesIterator->first ]->computeObservationsWithPartials(
                            simulationInputTime, dataIterator->first, dataIterator->second.second.second );

                // Compute residuals for current link ends and observabel type.
                residualsAndPartials.first.segment( startIndex, dataIterator->second.first.size( ) ) =
                        ( dataIterator->second.first - observationsWithPartials.first ).template cast< double >( );


                // Set current observation partials in matrix of all partials
                residualsAndPartials.second.block( startIndex, 0, dataIterator->second.first.size( ), parameterVectorSize ) =
                        observationsWithPartials.second;
                // Increment current index of observation.
                startIndex += dataIterator->second.first.size( );

            }
        }
    }

    Eigen::VectorXd normalizeObservationMatrix( Eigen::MatrixXd& observationMatrix )
    {
        Eigen::VectorXd range = Eigen::VectorXd( observationMatrix.cols( ) );

        for( int i = 0; i < observationMatrix.cols( ); i++ )
        {
            Eigen::VectorXd currentVector = observationMatrix.block( 0, i, observationMatrix.rows( ), 1 );
            double minimum = currentVector.minCoeff( );
            double maximum = currentVector.maxCoeff( );
            if( std::fabs( minimum ) > maximum )
            {
                range( i ) = minimum;
            }
            else
            {
                range( i ) = maximum;
            }
            currentVector = currentVector / range( i );
            //for( int j = 0; j < currentVector.rows( ); j++ )
            //{
            //    currentVector( j ) = ( currentVector( j ) - minimum ) / range( i );
            //}
            //currentVector = currentVector - minimum;
            //currentVector = currentVector / range;
            observationMatrix.block( 0, i, observationMatrix.rows( ), 1 ) = currentVector;
        }
        return range;
    }

    //! Function to perform parameter estimation from measurement data.
    /*!
     *  Function to perform parameter estimation, including orbit determination, i.e. body initial states, from measurement data.
     *  All observable types and link ends per obsevable types that are included in the measurement data input must have been
     *  provided to the constructor by the linkEndsPerObservable parameter.
     *  \param podInput Object containing all measurement data, associated metadata, including measurement weight, and a priori
     *  estimate for covariance matrix and parameter adjustment.
     *  \return Object containing estimated parameter value and associated data, such as residuals and observation partials.
     */
    boost::shared_ptr< PodOutput< StateScalarType > > estimateParameters(
            const boost::shared_ptr< PodInput< ObservationScalarType, TimeType, StateScalarType > >& podInput,
            const int reintegrateVariationalEquations = 1,
            const boost::shared_ptr< EstimationConvergenceChecker > convergenceChecker = ( boost::make_shared< EstimationConvergenceChecker >( ) ),
            const bool saveInformationmatrix = 1,
            const bool printOutput = 1 )
    {

        // Get size of parameter vector and number of observations (total and per type)
        int parameterVectorSize = currentParameterEstimate_.size( );
        std::pair< std::map< observation_models::ObservableType, int >, int > observationNumberPair =
                getNumberOfObservationsPerObservable( podInput->observationsAndTimes_ );
        int totalNumberOfObservations = observationNumberPair.second;

        // Declare variables to be returned (i.e. results from best iteration)
        double bestResidual = 1.0E100;
        ParameterVectorType bestParameterEstimate = ParameterVectorType::Zero( parameterVectorSize );
        Eigen::VectorXd bestTransformationData = Eigen::VectorXd::Zero( parameterVectorSize );
        Eigen::VectorXd bestResiduals = Eigen::VectorXd::Zero( totalNumberOfObservations );
        Eigen::MatrixXd bestInformationMatrix = Eigen::MatrixXd::Zero( totalNumberOfObservations, parameterVectorSize );
        Eigen::VectorXd bestWeightsMatrixDiagonal = Eigen::VectorXd::Zero( totalNumberOfObservations );
        Eigen::MatrixXd bestInverseNormalizedCovarianceMatrix = Eigen::MatrixXd::Zero( parameterVectorSize, parameterVectorSize );

        // Declare residual bookkeeping variables
        std::vector< double > rmsResidualHistory;
        double residualRms;

        // Declare variables to be used in loop.

        // Set current parameter estimate as both previous and current estimate
        ParameterVectorType newParameterEstimate = currentParameterEstimate_ + podInput->initialParameterDeviationEstimate_;
        ParameterVectorType oldParameterEstimate = currentParameterEstimate_;

        int numberOfEstimatedParameters = parameterVectorSize;//parametersToEstimate_->getEstimatedParameterSetSize( );

        // Iterate until convergence (at least once)
        int numberOfIterations = 0;
        do
        {
            // Re-integrate equations of motion and variational equations with new parameter estimate.
            if( ( numberOfIterations > 0 ) ||( podInput->reintegrateEquationsOnFirstIteration_ ) )
            {
                resetParameterEstimate( newParameterEstimate, reintegrateVariationalEquations );
            }
            oldParameterEstimate = newParameterEstimate;

            if( printOutput )
            {
                std::cout<<"Calculating residuals and partials "<<totalNumberOfObservations<<std::endl;
            }
            // Calculate residuals and observation matrix for current parameter estimate.
            std::pair< Eigen::VectorXd, Eigen::MatrixXd > residualsAndPartials;
            calculateObservationMatrixAndResiduals(
                        podInput->observationsAndTimes_, parameterVectorSize, totalNumberOfObservations, residualsAndPartials );

            Eigen::VectorXd transformationData = normalizeObservationMatrix( residualsAndPartials.second );

            Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix = Eigen::MatrixXd::Zero(
                        numberOfEstimatedParameters, numberOfEstimatedParameters );

            for( int j = 0; j < numberOfEstimatedParameters; j++ )
            {
                normalizedInverseAprioriCovarianceMatrix( j, j ) = podInput->inverseOfAprioriCovariance_( j, j ) /
                        ( transformationData( j ) * transformationData( j ) );
            }

            input_output::writeMatrixToFile(  residualsAndPartials.second,  "obsPartials.dat" );

            // Perform least squares calculation for correction to parameter vector.
            std::pair< Eigen::VectorXd, Eigen::MatrixXd > leastSquaresOutput =
                    linear_algebra::performLeastSquaresAdjustmentFromInformationMatrix(
                        residualsAndPartials.second.block( 0, 0, residualsAndPartials.second.rows( ), numberOfEstimatedParameters ),
                        residualsAndPartials.first, getConcatenatedWeightsVector( podInput->weightsMatrixDiagonals_ ),
                        normalizedInverseAprioriCovarianceMatrix, Eigen::VectorXd::Zero( numberOfEstimatedParameters, 0.0 ) );
            ParameterVectorType parameterAddition =
                    ( leastSquaresOutput.first.cwiseQuotient( transformationData.segment( 0, numberOfEstimatedParameters ) ) ).
                    template cast< StateScalarType >( );

            // Update value of parameter vector
            newParameterEstimate = oldParameterEstimate + parameterAddition;
            oldParameterEstimate = newParameterEstimate;

            if( printOutput )
            {
                std::cout<<"Parameter update"<<parameterAddition.transpose( )<<std::endl;
            }
            // Calculate mean residual for current iteration.
            residualRms = 0.0;
            double residualSum = 1.0;
            for( int i = 0; i < residualsAndPartials.first.size( ); i++ )
            {
                residualRms += std::fabs( residualsAndPartials.first[ i ] );// * weightMatrix( i, i );
                residualSum += 1.0; //weightMatrix( i, i );
            }
            residualRms = residualRms / ( residualSum );
            rmsResidualHistory.push_back( residualRms );
            if( printOutput )
            {
                std::cout<<"Current residual: "<<residualRms<<std::endl;
            }

            // If current iteration is better than previous one, update 'best' data.
            if( residualRms < bestResidual )
            {
                bestResidual = residualRms;
                bestParameterEstimate = oldParameterEstimate;
                bestResiduals = residualsAndPartials.first;
                if( saveInformationmatrix )
                {
                    bestInformationMatrix = residualsAndPartials.second;
                }
                bestWeightsMatrixDiagonal = getConcatenatedWeightsVector( podInput->weightsMatrixDiagonals_ );
                bestTransformationData = transformationData;
                bestInverseNormalizedCovarianceMatrix = leastSquaresOutput.second;
            }


            // Increment number of iterations
            numberOfIterations++;

            // Check for convergence
        } while( convergenceChecker->isEstimationConverged( numberOfIterations, rmsResidualHistory ) == false );


        std::cout<<"Final residual: "<<bestResidual<<std::endl;

        return boost::make_shared< PodOutput< StateScalarType > >(
                    bestParameterEstimate, bestResiduals, bestInformationMatrix, bestWeightsMatrixDiagonal, bestTransformationData,
                    bestInverseNormalizedCovarianceMatrix, bestResidual );
    }

    //! Function to reset the current parameter estimate.
    /*!
     *  Function to reset the current parameter estimate; reintegrates the variational equations and equations of motion with new estimate.
     *  \param newParameterEstimate New estimate of parameter vector.
     */
    void resetParameterEstimate( const ParameterVectorType& newParameterEstimate, const bool reintegrateVariationalEquations = 1 )
    {
        if( integrateAndEstimateOrbit_ )
        {
            variationalEquationsSolver_->resetParameterEstimate( newParameterEstimate, reintegrateVariationalEquations );
        }
        else
        {
            parametersToEstimate_->template resetParameterValues< StateScalarType>( newParameterEstimate );
        }
        currentParameterEstimate_ = newParameterEstimate;
    }

    //! Function to convert from one representation of all measurement data to the other
    /*!
     *  Function to convert from one representation of all measurement data (AlternativePodInputType) to the other (PodInputType).
     *  In the former, the vector of times and vector of associated observations are stored separately.  In the latter, they are
     *  stored as a map with time as key and obsevation as value.
     *  \param alternativePodInput Original representation of measurement data
     *  \return Converted measurement data.
     */
    static PodInputType convertPodInput( const AlternativePodInputType& alternativePodInput )
    {
        // Declare return data structure.
        PodInputType convertedObservationAndTimes;

        // Iterate over all observable types.
        for( typename AlternativePodInputType::const_iterator inputIterator = alternativePodInput.begin( );
             inputIterator != alternativePodInput.end( ); inputIterator++ )
        {
            // Declare data structure for converted measurement data at single observable type.
            std::map< observation_models::LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > singleTypeObservations;

            // Iterate over all link ends in current observable types
            for( typename std::map< observation_models::LinkEnds, std::pair< std::map< TimeType, ObservationScalarType >, observation_models::LinkEndType > >::const_iterator stationIterator =
                 inputIterator->second.begin( ); stationIterator != inputIterator->second.end( ); stationIterator++ )
            {
                // Initialize vector of observation times.
                std::vector< TimeType > times;
                times.resize( stationIterator->second.first.size( ) );

                // Initialize vector of observation values.
                ObservationVectorType observations = ObservationVectorType::Zero( stationIterator->second.first.size( ) );

                std::cerr<<"Warning, conversion between POD input types does not support non-double observables"<<std::endl;

                // Iterate over all observations in input map, and split time and observation.
                int counter = 0;
                for( typename std::map< TimeType, ObservationScalarType >::const_iterator singleDataSetIterator = stationIterator->second.first.begin( );
                     singleDataSetIterator != stationIterator->second.first.end( ); singleDataSetIterator++ )
                {
                    times[ counter ] = singleDataSetIterator->first;
                    observations( counter ) = singleDataSetIterator->second;
                    counter++;
                }

                // Set converted data for current link nends and observable.
                singleTypeObservations[ stationIterator->first ] = std::make_pair( observations, std::make_pair( times, stationIterator->second.second ) );
            }

            // Set converted data for current observable.
            convertedObservationAndTimes[ inputIterator->first ] = singleTypeObservations;
        }
        return convertedObservationAndTimes;
    }

    //! Function to retrieve the object to numerical integrate and update the variational equations and equations of motion
    /*!
     *  Function to retrieve the object to numerical integrate and update the variational equations and equations of motion
     *  \return Object to numerical integrate and update the variational equations and equations of motion
     */
    boost::shared_ptr< propagators::VariationalEquationsSolver< StateScalarType, TimeType, StateScalarType > >
    getVariationalEquationsSolver( ) const
    {
        return variationalEquationsSolver_;
    }

    //! Function to retrieve an observation manager
    /*!
     *  Function to retrieve an observation manager for a single observable type. The observation manager can simulate observations and
     *  calculate observation partials for all link ends involved in the given observable type.
     *  \param observation_models::ObservableType Type of observable for which manager is to be retrieved.
     *  \return Observation manager for given observable type.
     */
    boost::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType > > getObservationManager(
            const observation_models::ObservableType observableType ) const
    {
        // Check if manager exists for requested observable type.
        if( observationManagers_.count( observableType ) == 0 )
        {
            std::cerr<<"Error when retrieving observation manager of type "<<observableType<<", manager not found"<<std::endl;
        }

        return observationManagers_.at( observableType );
    }

    //! Function to retrieve the complete set of link ends for all observables.
    /*!
     *  Function to retrieve the complete set of link ends for all observables.
     *  \return Complete set of link ends for all observables.
     */
    std::map< observation_models::ObservableType, std::vector< observation_models::LinkEnds > > getLinkEndsPerObservable( )
    {
        return linkEndsPerObservable_;
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

    boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > getStateTransitionAndSensitivityMatrixInterface( )
    {
        return stateTransitionAndSensitivityMatrixInterface_;
    }


protected:

    bool integrateAndEstimateOrbit_;

    boost::shared_ptr< propagators::VariationalEquationsSolver< StateScalarType, TimeType, StateScalarType > > variationalEquationsSolver_;

    std::map< observation_models::ObservableType, boost::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType, StateScalarType > > >
    observationManagers_;

    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate_;

    std::map< observation_models::ObservableType, std::vector< observation_models::LinkEnds > > linkEndsPerObservable_;

    ParameterVectorType currentParameterEstimate_;

    boost::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionAndSensitivityMatrixInterface_;

};



}

}
#endif // ORBITDETERMINATIONMANAGER_H
