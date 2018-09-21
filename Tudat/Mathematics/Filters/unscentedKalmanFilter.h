/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Wan, E. and Van Der Merwe, R., “The Unscented Kalman Filter for Nonlinear Estimation,” in Adaptive Systems
 *          for Signal Processing, Communications, and Control Symposium. Institute of Electrical and Electronics
 *          Engineers, 2000, pp. 153–158.
 *      Jah, M., Lisano, M., Born, G., and Axelrad, P., “Mars Aerobraking Spacecraft State Estimation By Processing
 *          Inertial Measurement Unit Data,” Journal of Guidance, Control, and Dynamics, vol. 31, no. 6, pp. 1802–1812,
 *          November–December 2008.
 *      Challa, M., Moore, J., and Rogers, D., “A Simple Attitude Unscented Kalman Filter: Theory and Evaluation in
 *          a Magnetometer-Only Spacecraft Scenario,” IEEE Access, vol. 4, pp. 1845–1858, 2016.
 *      Vittaldev, V. (2010). The unified state model: Derivation and application in astrodynamics
 *          and navigation. Master's thesis, Delft University of Technology.
 */

#ifndef TUDAT_UNSCENTED_KALMAN_FILTER_H
#define TUDAT_UNSCENTED_KALMAN_FILTER_H

#include "Tudat/Mathematics/Filters/kalmanFilter.h"

namespace tudat
{

namespace filters
{

//! Enumeration for value of contant parameters.
enum ConstantParameterIndices
{
    alpha_index = 0,
    beta_index = 1,
    gamma_index = 2,
    kappa_index = 3,
    lambda_index = 4
};

//! Enumeration for value of contant parameters.
enum ConstantParameterReferences
{
    reference_Wan_and_Van_der_Merwe = 0,        // reference [Wan, E., et al.]
    reference_Lisano_and_Born_and_Axelrad = 1,  // reference [Jah, M., et al.]
    reference_Challa_and_Moore_and_Rogers = 2,  // reference [Challa, M., et al.]
    custom_parameters = 3
};

//! Unscented Kalman filter class.
/*!
 *  Class for the set up and use of the unscented Kalman filter.
 *  \tparam IndependentVariableType Type of independent variable. Default is double.
 *  \tparam DependentVariableType Type of dependent variable. Default is double.
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class UnscentedKalmanFilter: public KalmanFilterBase< IndependentVariableType, DependentVariableType >
{
public:

    //! Inherit typedefs from base class.
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::DependentVector DependentVector;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::DependentMatrix DependentMatrix;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::Function Function;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::IntegratorSettings IntegratorSettings;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::Integrator Integrator;

    //! Default constructor.
    /*!
     *  Default constructor. This constructor takes the system and measurement functions as models for the simulation.
     *  These functions can be a function of time and state vector.
     *  \param systemFunction Function returning the state as a function of time and state vector. Can be a differential
     *      equation if the integratorSettings is set (i.e., if it is not a nullptr).
     *  \param measurementFunction Function returning the measurement as a function of time and state.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param filteringStepSize Scalar representing the value of the constant filtering time step.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param integratorSettings Pointer to integration settings defining the integrator to be used to propagate the state.
     *  \param constantValueReference Reference to be used for the values of the \f$ \alpha \f$ and \f$ \kappa \f$ parameters. This
     *      variable has to be part of the ConstantParameterReferences enumeration (custom parameters are supported).
     *  \param customConstantParameters Values of the constant parameters \f$ \alpha \f$ and \f$ \kappa \f$, in case the custom_parameters
     *      enumeration is used in the previous field.
     */
    UnscentedKalmanFilter( const Function& systemFunction,
                           const Function& measurementFunction,
                           const DependentMatrix& systemUncertainty,
                           const DependentMatrix& measurementUncertainty,
                           const IndependentVariableType filteringStepSize,
                           const IndependentVariableType initialTime,
                           const DependentVector& initialStateVector,
                           const DependentMatrix& initialCovarianceMatrix,
                           const std::shared_ptr< IntegratorSettings > integratorSettings = nullptr,
                           const ConstantParameterReferences constantValueReference = reference_Wan_and_Van_der_Merwe,
                           const std::pair< DependentVariableType, DependentVariableType > customConstantParameters =
            std::make_pair( static_cast< DependentVariableType >( TUDAT_NAN ),
                            static_cast< DependentVariableType >( TUDAT_NAN ) ) ) :
        KalmanFilterBase< IndependentVariableType, DependentVariableType >( systemUncertainty, measurementUncertainty,
                                                                            filteringStepSize, initialTime, initialStateVector,
                                                                            initialCovarianceMatrix, integratorSettings ),
        inputSystemFunction_( systemFunction ), inputMeasurementFunction_( measurementFunction )
    {
        // Set dimensions
        stateDimension_ = systemUncertainty.rows( );
        measurementDimension_ = measurementUncertainty.rows( );

        // Set contant parameter values
        setConstantParameterValues( constantValueReference, customConstantParameters );

        // Generate weights for state and covariance estimations
        generateEstimationWeights( );

        // Create bases for augmented state vector and covariance matrix
        augmentedStateVector_ = DependentVector::Zero( augmentedStateDimension_ );
        augmentedCovarianceMatrix_ = DependentMatrix::Zero( augmentedStateDimension_, augmentedStateDimension_ );
        augmentedCovarianceMatrix_.block( stateDimension_, stateDimension_, stateDimension_, stateDimension_ ) = systemUncertainty;
        augmentedCovarianceMatrix_.block( 2 * stateDimension_, 2 * stateDimension_,
                                          measurementDimension_, measurementDimension_ ) = measurementUncertainty;
    }

    //! Destructor.
    ~UnscentedKalmanFilter( ){ }

    //! Function to update the filter with the new step data.
    /*!
     *  Function to update the filter with the new step data.
     *  \param currentMeasurementVector Vector representing current measurement.
     */
    void updateFilter( const DependentVector& currentMeasurementVector )
    {
        // Compute sigma points
        computeSigmaPoints( this->aPosterioriStateEstimate_, this->aPosterioriCovarianceEstimate_ );
        historyOfSigmaPoints_[ this->currentTime_ ] = mapOfSigmaPoints_; // store points

        // Prediction step
        // Compute series of state estimates based on sigma points
        std::map< unsigned int, DependentVector > sigmaPointsStateEstimates;
        for ( sigmaPointConstantIterator_ = mapOfSigmaPoints_.begin( );
              sigmaPointConstantIterator_ != mapOfSigmaPoints_.end( ); sigmaPointConstantIterator_++ )
        {
            currentSigmaPoint_ = sigmaPointConstantIterator_->first;
            sigmaPointsStateEstimates[ currentSigmaPoint_ ] = this->predictState(
                        sigmaPointConstantIterator_->second.segment( 0, stateDimension_ ) );
        }

        // Compute the weighted average to find the a-priori state vector
        DependentVector aPrioriStateEstimate = DependentVector::Zero( stateDimension_ );
        computeWeightedAverageFromSigmaPointEstimates( aPrioriStateEstimate, sigmaPointsStateEstimates );

        // Compute the weighted average to find the a-priori covariance matrix
        DependentMatrix aPrioriCovarianceEstimate = DependentMatrix::Zero( stateDimension_, stateDimension_ );
        computeWeightedAverageFromSigmaPointEstimates( aPrioriCovarianceEstimate, aPrioriStateEstimate, sigmaPointsStateEstimates );

        // Re-compute sigma points
        computeSigmaPoints( aPrioriStateEstimate, aPrioriCovarianceEstimate );

        // Compute series of measurement estimates based on sigma points
        std::map< unsigned int, DependentVector > sigmaPointsMeasurementEstimates;
        for ( sigmaPointConstantIterator_ = mapOfSigmaPoints_.begin( );
              sigmaPointConstantIterator_ != mapOfSigmaPoints_.end( ); sigmaPointConstantIterator_++ )
        {
            currentSigmaPoint_ = sigmaPointConstantIterator_->first;
            sigmaPointsMeasurementEstimates[ currentSigmaPoint_ ] = this->measurementFunction_(
                        this->currentTime_, sigmaPointConstantIterator_->second.segment( 0, stateDimension_ ) );
        }

        // Compute the weighted average to find the expected measurement vector
        DependentVector measurementEstimate = DependentVector::Zero( measurementDimension_ );
        computeWeightedAverageFromSigmaPointEstimates( measurementEstimate, sigmaPointsMeasurementEstimates );

        // Compute innovation and cross-correlation matrices
        DependentMatrix innovationMatrix = DependentMatrix::Zero( measurementDimension_, measurementDimension_ );
        computeWeightedAverageFromSigmaPointEstimates( innovationMatrix, measurementEstimate, sigmaPointsMeasurementEstimates );
        DependentMatrix crossCorrelationMatrix = DependentMatrix::Zero( stateDimension_, measurementDimension_ );
        for ( sigmaPointConstantIterator_ = sigmaPointsStateEstimates.begin( );
              sigmaPointConstantIterator_ != sigmaPointsStateEstimates.end( ); sigmaPointConstantIterator_++ )
        {
            crossCorrelationMatrix += covarianceEstimationWeights_.at( sigmaPointConstantIterator_->first ) *
                    ( sigmaPointConstantIterator_->second - aPrioriStateEstimate ) *
                    ( sigmaPointsMeasurementEstimates[ sigmaPointConstantIterator_->first ] - measurementEstimate ).transpose( );
        }

        // Compute Kalman gain
        DependentMatrix kalmanGain = crossCorrelationMatrix * innovationMatrix.inverse( );

        // Correction step
        this->currentTime_ += this->filteringStepSize_;
        this->correctState( aPrioriStateEstimate, currentMeasurementVector, measurementEstimate, kalmanGain );
        correctCovariance( aPrioriCovarianceEstimate, innovationMatrix, kalmanGain );
    }

    //! Function to return the history of sigma points.
    /*!
     *  Function to return the history of sigma points.
     *  \return History of map of sigma points for each time step.
     */
    std::map< IndependentVariableType, DependentMatrix > getHistoryOfSigmaPoints( )
    {
        // Convert map of maps into map of Eigen::Matrix
        DependentVector vectorOfSigmaPoints;
        std::map< IndependentVariableType, DependentMatrix > mapOfSigmaPointsHistory;
        for ( typename std::map< IndependentVariableType, std::map< unsigned int, DependentVector > >::const_iterator
              mapIterator = historyOfSigmaPoints_.begin( ); mapIterator != historyOfSigmaPoints_.end( ); mapIterator++ )
        {
            // Extract current map of sigma points and turn it into a vector
            vectorOfSigmaPoints = utilities::createConcatenatedEigenMatrixFromMapValues( mapIterator->second );

            // Reshape the vector into a matrix and store it into the output map
            Eigen::Map< Eigen::MatrixXd > matrixOfSigmaPoints( vectorOfSigmaPoints.data( ),
                                                               augmentedStateDimension_, numberOfSigmaPoints_ );
            mapOfSigmaPointsHistory[ mapIterator->first ] = matrixOfSigmaPoints;
        }

        // Give output
        return mapOfSigmaPointsHistory;
    }

private:

    //! Function to create the function that defines the system model.
    /*!
     *  Function to create the function that defines the system model. The output of this function is then bound
     *  to the systemFunction_ variable, via the std::bind command.
     *  \param currentTime Scalar representing the current time.
     *  \param currentStateVector Vector representing the current state.
     *  \return Vector representing the estimated state.
     */
    DependentVector createSystemFunction( const IndependentVariableType currentTime,
                                          const DependentVector& currentStateVector )
    {
        return inputSystemFunction_( currentTime, currentStateVector ) +
                mapOfSigmaPoints_[ currentSigmaPoint_ ].segment( stateDimension_, stateDimension_ ); // add system noise
    }

    //! Function to create the function that defines the system model.
    /*!
     *  Function to create the function that defines the system model. The output of this function is then bound
     *  to the measurementFunction_ variable, via the std::bind command.
     *  \param currentTime Scalar representing the current time.
     *  \param currentStateVector Vector representing the current state.
     *  \return Vector representing the estimated measurement.
     */
    DependentVector createMeasurementFunction( const IndependentVariableType currentTime,
                                               const DependentVector& currentStateVector )
    {
        return inputMeasurementFunction_( currentTime, currentStateVector ) +
                mapOfSigmaPoints_[ currentSigmaPoint_ ].segment( 2 * stateDimension_, measurementDimension_ ); // add measurement noise
    }

    //! Function to clear the history of stored variables for derived class-specific variables.
    /*!
     *  Function to clear the history of stored variables for derived class-specific variables. This function adds to the list of
     *  variables to be cleared, the history of sigma point estimates over time.
     */
    void clearSpecificFilterHistory( ) { historyOfSigmaPoints_.clear( ); }

    //! Function to revert to the previous time step for derived class-specific variables.
    /*!
     *  Function to revert to the previous time step for derived class-specific variables. This function adds to the list of
     *  elements to be reverted, the latest sigma point estimates.
     *  \param timeToBeRemoved Double denoting the current time, i.e., the instant that has to be discarded.
     */
    virtual void specificRevertToPreviousTimeStep( const double timeToBeRemoved )
    {
        if ( historyOfSigmaPoints_.count( timeToBeRemoved ) != 0 )
        {
            historyOfSigmaPoints_.erase( timeToBeRemoved );
        }
    }

    //! Function to set the values of the constant parameters.
    /*!
     *  Function to set the values of the constant parameters, used by the unscented Kalman filter for various purposes.
     *  These parameters are saved in the constantParameters_ vector, and are retrieved by using the ConstantParameterIndices
     *  enumeration. The order and definition of the parameters is the following [Vittaldev, V.]:
     *      - alpha_index (\f$ \alpha \f$): used to distribute the sigma points around the a-priori estimate.
     *      - beta_index (\f$ \beta \f$): provides information about the probaility distribution function of the state.
     *      - gamma_index (\f$ \gamma \f$): abbreviation for \f$ \sqrt{ L + \lambda } \f$, where \f$ L \f$ is the length of the
     *          augmented state vector.
     *      - kappa_index (\f$ \kappa \f$): secondary scaling parameter, also used to distribute the sigma points around the
     *          a-priori estimate, but it has a smaller influence.
     *      - lambda_index (\f$ \lambda \f$): scaling parameter, also used to distribute the sigma points around the a-priori
     *          estimate.
     *  \param constantValueReference Reference to be used for the values of the alpha and kappa parameters. This variable has to
     *      be part of the ConstantParameterReferences enumeration (custom parameters are supported).
     *  \param customConstantParameters Values of the constant parameters alpha and kappa, in case the custom_parameters enumerate
     *      is used in the previous field.
     */
    void setConstantParameterValues( const ConstantParameterReferences constantValueReference,
                                     const std::pair< DependentVariableType, DependentVariableType >& customConstantParameters );

    //! Function to generate the weights for state and covariance estimation.
    /*!
     *  Function to generate the weights for state and covariance estimation, which will be used to determine the weighted average
     *  of the state and measurment vectors, and covariance matrix, based on the sigma points.
     */
    void generateEstimationWeights( );

    //! Function to compute the sigma points, based on the current state vector and covariance matrix.
    /*!
     *  Function to compute the sigma points, based on the current a-priori or previous step a-posteriori state vector and
     *  covariance matrix estimates. The sigma points are spread around the current a-priori state estimate, and their propagation
     *  is used to determine the sensitivity of the state model to changes in initial conditions. These offsets are then used to compute
     *  the new state and measurement vectors and covariance matrix estimates.
     */
    void computeSigmaPoints( const DependentVector& currentStateEstimate, const DependentMatrix& currentCovarianceEstimate )
    {
        // Update augmented state and covariance matrix to new values
        augmentedStateVector_.segment( 0, stateDimension_ ) = currentStateEstimate;
        augmentedCovarianceMatrix_.topLeftCorner( stateDimension_, stateDimension_ ) = currentCovarianceEstimate;

        // Pre-compute square root of augmented covariance matrix
        DependentMatrix augmentedCovarianceMatrixSquareRoot;
        if ( augmentedCovarianceMatrix_.determinant( ) != static_cast< DependentVariableType >( 0.0 ) )
        {
            // Matrix is invertible, so take matrix square-root directly
            augmentedCovarianceMatrixSquareRoot = augmentedCovarianceMatrix_.sqrt( );
        }
        else
        {
            // Check if matrix is diagonal
            DependentVector augmentedCovarianceMatrixDiagonal = augmentedCovarianceMatrix_.diagonal( );
            if ( augmentedCovarianceMatrix_.isApprox( DependentMatrix( augmentedCovarianceMatrixDiagonal.asDiagonal( ) ) ) )
            {
                augmentedCovarianceMatrixDiagonal = augmentedCovarianceMatrixDiagonal.array( ).sqrt( );
                augmentedCovarianceMatrixSquareRoot = augmentedCovarianceMatrixDiagonal.asDiagonal( );
            }
            else
            {
                // Otherwise, matrix is singular and the square-root cannot be computed
                throw std::runtime_error( "Error in unscented Kalman filter. Augemented covariance matrix is singular. "
                                          "Its square-root cannot be computed." );
            }
        }

        // Loop over sigma points and assign value
        for ( unsigned int i = 0; i < numberOfSigmaPoints_; i++ )
        {
            if ( i == 0 )
            {
                mapOfSigmaPoints_[ i ] = augmentedStateVector_;
            }
            else if ( i < ( augmentedCovarianceMatrixSquareRoot.cols( ) + 1 ) )
            {
                mapOfSigmaPoints_[ i ] = augmentedStateVector_ + constantParameters_.at( gamma_index ) *
                        augmentedCovarianceMatrixSquareRoot.col( i - 1 );
            }
            else
            {
                mapOfSigmaPoints_[ i ] = augmentedStateVector_ - constantParameters_.at( gamma_index ) *
                        augmentedCovarianceMatrixSquareRoot.col( ( i - 1 ) - augmentedCovarianceMatrixSquareRoot.cols( ) );
            }
        }
    }

    //! Function to compute the weighted average of the state and measurement vectors.
    /*!
     *  Function to compute the weighted average of the state and measurement vectors.
     *  \param weightedAverageVector Vector to which the weighted average is added (initially set to zero).
     *  \param sigmaPointEstimates Map of the sigma points generated by the computeSigmaPoints function.
     *  \return Weighted average of the state or measurement vector, i.e., the new a-priori state and the
     *      measurement estimates (returned by reference).
     */
    void computeWeightedAverageFromSigmaPointEstimates( DependentVector& weightedAverageVector,
                                                        const std::map< unsigned int, DependentVector >& sigmaPointEstimates )
    {
        // Loop over each sigma point
        for ( sigmaPointConstantIterator_ = sigmaPointEstimates.begin( );
              sigmaPointConstantIterator_ != sigmaPointEstimates.end( ); sigmaPointConstantIterator_++ )
        {
            weightedAverageVector += stateEstimationWeights_.at( sigmaPointConstantIterator_->first ) *
                    sigmaPointConstantIterator_->second;
        }
    }

    //! Function to compute the weighted average of the covariance and innovation matrices.
    /*!
     *  Function to compute the weighted average of the covariance and innovation matrices.
     *  \param weightedAverageMatrix Matrix to which the weighted average is added (initially set to zero).
     *  \param referenceVector Vector representing the a-priori state or measurement estimates.
     *  \param sigmaPointEstimates Map of the sigma points generated by the computeSigmaPoints function.
     *  \return Weighted average of the covariance and innovation matrices, i.e., the new a-priori covariance and the
     *      innovation estimates (returned by reference).
     */
    void computeWeightedAverageFromSigmaPointEstimates( DependentMatrix& weightedAverageMatrix,
                                                        const DependentVector& referenceVector,
                                                        const std::map< unsigned int, DependentVector >& sigmaPointEstimates )
    {
        // Loop over each sigma point
        for ( sigmaPointConstantIterator_ = sigmaPointEstimates.begin( );
              sigmaPointConstantIterator_ != sigmaPointEstimates.end( ); sigmaPointConstantIterator_++ )
        {
            weightedAverageMatrix += covarianceEstimationWeights_.at( sigmaPointConstantIterator_->first ) *
                    ( sigmaPointConstantIterator_->second - referenceVector ) *
                    ( sigmaPointConstantIterator_->second - referenceVector ).transpose( );
        }
    }

    //! Function to correct the covariance for the next time step.
    /*!
     *  Function to correct the covariance for the next time step.
     *  \param aPrioriCovarianceEstimate Matrix denoting the a-priori covariance estimate.
     *  \param innovationMatrix Matrix denoting the innovation, i.e., the correlation between system and measurement.
     *  \param kalmanGain Matrix denoting the Kalman gain, to be used to correct the state estimate with the external measurement data.
     */
    void correctCovariance( const DependentMatrix& aPrioriCovarianceEstimate,
                            const DependentMatrix& innovationMatrix, const DependentMatrix& kalmanGain )
    {
        this->aPosterioriCovarianceEstimate_ = aPrioriCovarianceEstimate - kalmanGain * innovationMatrix * kalmanGain.transpose( );
        this->historyOfCovarianceEstimates_[ this->currentTime_ ] = this->aPosterioriCovarianceEstimate_;
    }

    //! System function input by user.
    Function inputSystemFunction_;

    //! Measurement function input by user.
    Function inputMeasurementFunction_;

    //! Integer specifying length of state vector.
    unsigned int stateDimension_;

    //! Integer specifying length of measurement vector.
    unsigned int measurementDimension_;

    //! Integer specifying length of augmented state vector.
    unsigned int augmentedStateDimension_;

    //! Integer specifying number of sigma points.
    unsigned int numberOfSigmaPoints_;

    //! Vector of constant parameters.
    /*!
     *  Vector of constant parameters. See description of setConstantParameterValues function for information on the order and
     *  meaning of the constant parameters.
     */
    std::vector< DependentVariableType > constantParameters_;

    //! Vector of weights used for the computation of the weighted average of the state and measurement vectors.
    std::vector< DependentVariableType > stateEstimationWeights_;

    //! Vector of weights used for the computation of the weighted average of the covariance and innovation matrices.
    std::vector< DependentVariableType > covarianceEstimationWeights_;

    //! Augmented state vector.
    /*!
     *  Augmented state vector, which is defined by vertically concatenating the state vector and the expectations of the
     *  state and measurement noises (which for are defined to be zero in this application).
     */
    DependentVector augmentedStateVector_;

    //! Augmented covariance matrix.
    /*!
     *  Augmented covariance matrix, which is defined by diagonally concatenating the state covariance matrix and the
     *  state and measurement uncertainties (which are assumed to be constant and are provided by the user).
     */
    DependentMatrix augmentedCovarianceMatrix_;

    //! Map of sigma points.
    /*!
     *  Map of sigma points, as output by the computeSigmaPoints function. See the description of this function for more
     *  details of the sigma points and their use.
     */
    std::map< unsigned int, DependentVector > mapOfSigmaPoints_;

    //! Map of map of sigma points, used to store the history of sigma points.
    std::map< IndependentVariableType, std::map< unsigned int, DependentVector > > historyOfSigmaPoints_;

    //! Constant iterator to loop over sigma points (introduced for convenience).
    typename std::map< unsigned int, DependentVector >::const_iterator sigmaPointConstantIterator_;

    //! Integer specifying current sigma point.
    /*!
     *  Integer specifying current sigma point, while iterating over the mapOfSigmaPoints_. This parameter is specifically used
     *  when evaluating the systemFunction_ and measurementFunction_, such that the correct value of system and measurement
     *  noise can be added.
     */
    unsigned int currentSigmaPoint_;

};

//! Typedef for a filter with double data type.
typedef UnscentedKalmanFilter< > UnscentedKalmanFilterDouble;

//! Typedef for a shared-pointer to a filter with double data type.
typedef std::shared_ptr< UnscentedKalmanFilterDouble > UnscentedKalmanFilterDoublePointer;

//! Function to set the values of the constant parameters.
template< typename IndependentVariableType, typename DependentVariableType >
void UnscentedKalmanFilter< IndependentVariableType, DependentVariableType >::setConstantParameterValues(
        const ConstantParameterReferences constantValueReference,
        const std::pair< DependentVariableType, DependentVariableType >& customConstantParameters )
{
    // Assign size to vector
    constantParameters_.resize( 5 );

    // Set parameters based on input
    switch ( constantValueReference )
    {
    case reference_Wan_and_Van_der_Merwe:
    {
        constantParameters_.at( alpha_index ) = static_cast< DependentVariableType >( 0.003 );
        constantParameters_.at( kappa_index ) = static_cast< DependentVariableType >( 0.0 );
        break;
    }
    case reference_Lisano_and_Born_and_Axelrad:
    {
        constantParameters_.at( alpha_index ) = static_cast< DependentVariableType >( 1.0 );
        constantParameters_.at( kappa_index ) = static_cast< DependentVariableType >( 3.0 - stateDimension_ );
        break;
    }
    case reference_Challa_and_Moore_and_Rogers:
    {
        constantParameters_.at( alpha_index ) = static_cast< DependentVariableType >( 0.001 );
        constantParameters_.at( kappa_index ) = static_cast< DependentVariableType >( 1.0 );
        break;
    }
    case custom_parameters:
    {
        // Check that the values have been set
        if ( customConstantParameters.first == static_cast< DependentVariableType >( TUDAT_NAN ) ||
             customConstantParameters.second == static_cast< DependentVariableType >( TUDAT_NAN ) )
        {
            throw std::runtime_error( "Error in unscented Kalman filter. The value of the alpha and kappa parameters "
                                      "have not been specified, but the selected method is custom_parameters." );
        }

        // Assign values to parameters
        constantParameters_.at( alpha_index ) = customConstantParameters.first;
        constantParameters_.at( kappa_index ) = customConstantParameters.second;
        break;
    }
    default:
        throw std::runtime_error( "Error in unscented Kalman filter. The name of the reference for the alpha and kappa "
                                  "parameters is not recognized. To enter a custom pair of coefficients, "
                                  "use the value custom_parameters." );
    }

    // Set augmented state and sigma parameters
    augmentedStateDimension_ = 2 * stateDimension_ + measurementDimension_;
    numberOfSigmaPoints_ = 2 * augmentedStateDimension_ + 1;

    // Set remaining parameters
    constantParameters_.at( beta_index ) = 2.0;
    constantParameters_.at( lambda_index ) = std::pow( constantParameters_.at( alpha_index ), 2 ) *
            ( augmentedStateDimension_ + constantParameters_.at( kappa_index ) ) - augmentedStateDimension_;
    constantParameters_.at( gamma_index ) = std::sqrt( augmentedStateDimension_ + constantParameters_.at( lambda_index ) );
}

//! Function to generate the weights for state and covariance estimation.
template< typename IndependentVariableType, typename DependentVariableType >
void UnscentedKalmanFilter< IndependentVariableType, DependentVariableType >::generateEstimationWeights( )
{
    // Generate state and covariance estimation weights
    stateEstimationWeights_.push_back( constantParameters_.at( lambda_index ) /
                                       ( augmentedStateDimension_ + constantParameters_.at( lambda_index ) ) );
    for ( unsigned int i = 1; i < numberOfSigmaPoints_; i++ )
    {
        stateEstimationWeights_.push_back( 1.0 / ( 2.0 * ( augmentedStateDimension_ + constantParameters_.at( lambda_index ) ) ) );
    }
    covarianceEstimationWeights_ = stateEstimationWeights_;
    covarianceEstimationWeights_.at( 0 ) += 1.0 - std::pow( constantParameters_.at( alpha_index ), 2 ) +
            constantParameters_.at( beta_index );
}

} // namespace filters

} // namespace tudat

#endif // TUDAT_UNSCENTED_KALMAN_FILTER_H
