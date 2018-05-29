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
 */

#ifndef TUDAT_UNSCENTED_KALMAN_FILTER_H
#define TUDAT_UNSCENTED_KALMAN_FILTER_H

#include <unsupported/Eigen/MatrixFunctions>

#include "Tudat/Mathematics/Filters/kalmanFilter.h"

namespace tudat
{

namespace filters
{

//! Enumeration for value of contant parameters
enum ConstantParameterIndices
{
    alpha_index = 0,
    beta_index = 1,
    gamma_index = 2,
    kappa_index = 3,
    lambda_index = 4
};

//! Enumeration for value of contant parameters
enum ConstantParameterReferences
{
    reference_Wan_and_Van_der_Merwe = 0,        // reference [1]
    reference_Lisano_and_Born_and_Axelrad = 1,  // reference [2]
    reference_Challa_and_Moore_and_Rogers = 2,  // reference [3]
    custom_parameters = 3
};

template< typename IndependentVariableType = double, typename DependentVariableType = double >
void printMapContents( const std::map< IndependentVariableType,
                       Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 > >& mapToPrint )
{
    for ( typename std::map< IndependentVariableType,
          Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 > >::const_iterator mapIterator = mapToPrint.begin( );
          mapIterator != mapToPrint.end( ); mapIterator++ )
    {
        std::cout << mapIterator->first << ", " << mapIterator->second.transpose( ) << std::endl;
    }
}

//! Unscented Kalman filter.
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
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::SystemFunction SystemFunction;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::MeasurementFunction MeasurementFunction;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::IntegratorSettings IntegratorSettings;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::Integrator Integrator;

    //! Constructor.
    /*!
     *  Constructor.
     */
    UnscentedKalmanFilter( const SystemFunction& systemFunction,
                           const MeasurementFunction& measurementFunction,
                           const DependentMatrix& systemUncertainty,
                           const DependentMatrix& measurementUncertainty,
                           const IndependentVariableType initialTime,
                           const DependentVector& initialStateVector,
                           const DependentMatrix& initialCovarianceMatrix,
                           const boost::shared_ptr< IntegratorSettings > integratorSettings = NULL,
                           const ConstantParameterReferences constantValueReference = reference_Wan_and_Van_der_Merwe,
                           const std::pair< DependentVariableType, DependentVariableType > customConstantParameters =
            std::make_pair( TUDAT_NAN, TUDAT_NAN ) ) :
        KalmanFilterBase< IndependentVariableType, DependentVariableType >( systemUncertainty, measurementUncertainty,
                                                                            initialTime, initialStateVector,
                                                                            initialCovarianceMatrix, integratorSettings ),
        inputSystemFunction_( systemFunction ), inputMeasurementFunction_( measurementFunction )
    {
        // Get and set dimensions
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

        std::vector< std::vector< DependentVariableType > > vectorOfVectors = { constantParameters_,
                                                                                stateEstimationWeights_, covarianceEstimationWeights_ };
        for ( std::vector< DependentVariableType > vector: vectorOfVectors )
        {
            for ( DependentVariableType element: vector )
            {
                std::cout << element << ", ";
            }
            std::cout << std::endl;
        }
    }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~UnscentedKalmanFilter( ){ }

    //! Function to update the filter with the new step data.
    /*!
     *  Function to update the filter with the new step data.
     *  \param currentTime Scalar representing current time.
     *  \param currentControlVector Vector representing the current control input.
     *  \param currentMeasurementVector Vector representing current measurement.
     */
    void updateFilter( const IndependentVariableType currentTime, const DependentVector& currentControlVector,
                       const DependentVector& currentMeasurementVector )
    {
        // Update augmented state vector and covariance matrix
        augmentedStateVector_.segment( 0, stateDimension_ ) = this->aPosterioriStateEstimate_;
        augmentedCovarianceMatrix_.topLeftCorner( stateDimension_, stateDimension_ ) = this->aPosterioriCovarianceEstimate_;
        std::cout << "Augmented state: " << augmentedStateVector_.transpose( ) << std::endl;
        std::cout << "Augmented covariance: " << std::endl << augmentedCovarianceMatrix_ << std::endl;

        // Create sigma points
        generateSigmaPoints( );
        std::cout << "Sigma points: " << std::endl;
        printMapContents( mapOfSigmaPoints_ );

        // Prediction step
        // Compute series of state estimates based on sigma points
        std::map< unsigned int, DependentVector > sigmaPointsStateEstimates;
        if ( this->isStateToBeIntegrated_ )
        {
            // Propagate each sigma point
            for ( sigmaPointConstantIterator_ = mapOfSigmaPoints_.begin( );
                  sigmaPointConstantIterator_ != mapOfSigmaPoints_.end( ); sigmaPointConstantIterator_++ )
            {
                currentSigmaPoint_ = sigmaPointConstantIterator_->first;
                sigmaPointsStateEstimates[ sigmaPointConstantIterator_->first ] = this->propagateState(
                            currentTime, sigmaPointConstantIterator_->second.segment( 0, stateDimension_ ),
                            currentControlVector );
            }
        }
        else
        {
            // Compute each sigma point
            for ( sigmaPointConstantIterator_ = mapOfSigmaPoints_.begin( );
                  sigmaPointConstantIterator_ != mapOfSigmaPoints_.end( ); sigmaPointConstantIterator_++ )
            {
                currentSigmaPoint_ = sigmaPointConstantIterator_->first;
                sigmaPointsStateEstimates[ sigmaPointConstantIterator_->first ] = this->systemFunction_(
                            currentTime, sigmaPointConstantIterator_->second.segment( 0, stateDimension_ ),
                            currentControlVector );
            }
        }
        std::cout << "Propgated sigma points: " << std::endl;
        printMapContents( sigmaPointsStateEstimates );

        // Compute the weighted average to find the a-priori state vector
        DependentVector aPrioriStateEstimate = DependentVector::Zero( stateDimension_ );
        computeWeightedAverageFromSigmaPointEstimates( aPrioriStateEstimate, sigmaPointsStateEstimates );
        std::cout << "x_k_k1: " << aPrioriStateEstimate.transpose( ) << std::endl;

        // Compute the a-priori covariance matrix
        DependentMatrix aPrioriCovarianceEstimate = DependentMatrix::Zero( stateDimension_, stateDimension_ );
        computeWeightedAverageFromSigmaPointEstimates( aPrioriCovarianceEstimate, aPrioriStateEstimate, sigmaPointsStateEstimates );
        std::cout << "P_k_k1: " << std::endl << aPrioriCovarianceEstimate << std::endl;

        // Re-generate sigma points
        augmentedStateVector_.segment( 0, stateDimension_ ) = this->aPosterioriStateEstimate_;
        augmentedCovarianceMatrix_.topLeftCorner( stateDimension_, stateDimension_ ) = this->aPosterioriCovarianceEstimate_;
        generateSigmaPoints( );

        // Compute series of measurement estimates based on sigma points
        std::map< unsigned int, DependentVector > sigmaPointsMeasurementEstimates;
        for ( sigmaPointConstantIterator_ = mapOfSigmaPoints_.begin( );
              sigmaPointConstantIterator_ != mapOfSigmaPoints_.end( ); sigmaPointConstantIterator_++ )
        {
            currentSigmaPoint_ = sigmaPointConstantIterator_->first;
            sigmaPointsMeasurementEstimates[ sigmaPointConstantIterator_->first ] = this->measurementFunction_(
                        currentTime, sigmaPointConstantIterator_->second.segment( 0, stateDimension_ ) );
        }
        std::cout << "Measurement points: " << std::endl;
        printMapContents( sigmaPointsMeasurementEstimates );

        // Compute the weighted average to find the expected measurement vector
        DependentVector measurmentEstimate = DependentVector::Zero( measurementDimension_ );
        computeWeightedAverageFromSigmaPointEstimates( measurmentEstimate, sigmaPointsMeasurementEstimates );
        std::cout << "z_k_k1: " << measurmentEstimate.transpose( ) << std::endl;

        // Compute innovation and cross-correlation matrices
        DependentMatrix innovationMatrix = DependentMatrix::Zero( measurementDimension_, measurementDimension_ );
        computeWeightedAverageFromSigmaPointEstimates( innovationMatrix, measurmentEstimate, sigmaPointsMeasurementEstimates );
        DependentMatrix crossCorrelationMatrix = DependentMatrix::Zero( stateDimension_, measurementDimension_ );
        for ( sigmaPointConstantIterator_ = sigmaPointsStateEstimates.begin( );
              sigmaPointConstantIterator_ != sigmaPointsStateEstimates.end( ); sigmaPointConstantIterator_++ )
        {
            crossCorrelationMatrix += covarianceEstimationWeights_.at( sigmaPointConstantIterator_->first ) *
                    ( sigmaPointConstantIterator_->second - aPrioriStateEstimate ) *
                    ( sigmaPointsMeasurementEstimates[ sigmaPointConstantIterator_->first ] - measurmentEstimate ).transpose( );
        }
        std::cout << "P_zz: " << std::endl << innovationMatrix - Eigen::Vector1d::Constant( 1.01 ) << std::endl;
        std::cout << "P_xz: " << std::endl << crossCorrelationMatrix << std::endl;

        // Compute Kalman gain
        DependentMatrix kalmanGain = crossCorrelationMatrix * innovationMatrix.inverse( );
        std::cout << "K: " << std::endl << kalmanGain << std::endl;

        // Correction step
        this->correctStateAndCovariance( currentTime, aPrioriStateEstimate, aPrioriCovarianceEstimate,
                                         DependentMatrix::Zero( measurementDimension_, stateDimension_ ),
                                         currentMeasurementVector, measurmentEstimate, kalmanGain, true, innovationMatrix );
        std::cout << "x_k_k: " << this->aPosterioriStateEstimate_.transpose( ) << std::endl;
        std::cout << "P_k_k: " << std::endl << this->aPosterioriCovarianceEstimate_ << std::endl;
        std::cout << std::endl;
    }

private:

    //! Function to create the function that defines the system model.
    /*!
     *  Function to create the function that defines the system model. The output of this function is then bound
     *  to the systemFunction_ variable, via the boost::bind command.
     *  \param currentTime Scalar representing the current time.
     *  \param currentStateVector Vector representing the current state.
     *  \param currentControlVector Vector representing the current control input.
     *  \return Vector representing the estimated state.
     */
    DependentVector createSystemFunction( const IndependentVariableType currentTime,
                                          const DependentVector& currentStateVector,
                                          const DependentVector& currentControlVector )
    {
        return inputSystemFunction_( currentTime, currentStateVector, currentControlVector ) +
                mapOfSigmaPoints_[ currentSigmaPoint_ ].segment( stateDimension_, stateDimension_ ); // add system noise
    }

    //! Function to create the function that defines the system model.
    /*!
     *  Function to create the function that defines the system model. The output of this function is then bound
     *  to the measurementFunction_ variable, via the boost::bind command.
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

    //!
    /*!
     *
     */
    void setConstantParameterValues( const ConstantParameterReferences constantValueReference,
                                     const std::pair< DependentVariableType, DependentVariableType >& customConstantParameters )
    {
        // Set parameters based on input
        switch ( constantValueReference )
        {
        case reference_Wan_and_Van_der_Merwe:
            constantParameters_.at( alpha_index ) = 0.003;
            constantParameters_.at( kappa_index ) = 0.0;
            break;
        case reference_Lisano_and_Born_and_Axelrad:
            constantParameters_.at( alpha_index ) = 1.0;
            constantParameters_.at( kappa_index ) = 3.0 - stateDimension_;
            break;
        case reference_Challa_and_Moore_and_Rogers:
            constantParameters_.at( alpha_index ) = 0.001;
            constantParameters_.at( kappa_index ) = 1.0;
            break;
        case custom_parameters:
        {
            // Check that the values have been set
            if ( customConstantParameters.first == TUDAT_NAN || customConstantParameters.second == TUDAT_NAN )
            {
                throw std::runtime_error( "Error in unscented Kalman filter. The value of the alpha and kappa parameters "
                                          "have not been specified, but the selected method is custom_parameters." );
            }

            // Assign values to parameters
            constantParameters_.at( alpha_index ) = customConstantParameters.first;
            constantParameters_.at( kappa_index ) = customConstantParameters.second;
            break;
        }
        }

        // Set augmented state and sigma parameters
        augmentedStateDimension_ = 2 * stateDimension_ + measurementDimension_;
        numberOfSigmaPoints_ = 2.0 * augmentedStateDimension_ + 1.0;

        // Set remaining parameters
        constantParameters_.at( beta_index ) = 2.0;
        constantParameters_.at( lambda_index ) = std::pow( constantParameters_.at( alpha_index ), 2 ) *
                ( augmentedStateDimension_ + constantParameters_.at( kappa_index ) ) - augmentedStateDimension_;
        constantParameters_.at( gamma_index ) = std::sqrt( augmentedStateDimension_ + constantParameters_.at( lambda_index ) );
    }

    //!
    /*!
     *
     */
    void generateEstimationWeights( )
    {
        // Generate state and covariance estimation weights
        stateEstimationWeights_.push_back( constantParameters_.at( lambda_index ) /
                                           ( augmentedStateDimension_ + constantParameters_.at( lambda_index ) ) );
        covarianceEstimationWeights_.push_back( stateEstimationWeights_.back( ) + 1.0 -
                                                std::pow( constantParameters_.at( alpha_index ), 2 ) +
                                                constantParameters_.at( beta_index ) );
        for ( unsigned int i = 1; i < numberOfSigmaPoints_; i++ )
        {
            stateEstimationWeights_.push_back( 1.0 / ( 2.0 * ( augmentedStateDimension_ + constantParameters_.at( lambda_index ) ) ) );
            covarianceEstimationWeights_.push_back( stateEstimationWeights_.back( ) );
        }
    }

    //!
    /*!
     *
     */
    void generateSigmaPoints( )
    {
        // Pre-compute square root of augmented covariance matrix
        DependentMatrix augmentedCovarianceMatrixSquareRoot = augmentedCovarianceMatrix_.sqrt( );

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

    //!
    /*!
     *
     */
    void computeWeightedAverageFromSigmaPointEstimates( DependentVector& weightedAverageVector,
                                                        const std::map< unsigned int, DependentVector >& sigmaPointEstimates )
    {
        // Loop over each sigma point
        for ( sigmaPointConstantIterator_ = sigmaPointEstimates.begin( );
              sigmaPointConstantIterator_ != sigmaPointEstimates.end( ); sigmaPointConstantIterator_++ )
        {
            weightedAverageVector += stateEstimationWeights_.at( sigmaPointConstantIterator_->first ) * sigmaPointConstantIterator_->second;
        }
    }

    //!
    /*!
     *
     */
    void computeWeightedAverageFromSigmaPointEstimates( DependentMatrix& weightedAverageMatrix,
                                                        const DependentVector& referenceVector,
                                                        const std::map< unsigned int, DependentVector >& sigmaPointEstimates )
    {
        // Loop over each sigma point
        for ( sigmaPointConstantIterator_ = sigmaPointEstimates.begin( );
              sigmaPointConstantIterator_ != sigmaPointEstimates.end( ); sigmaPointConstantIterator_++ )
        {
            std::cout << covarianceEstimationWeights_.at( sigmaPointConstantIterator_->first ) *
                         ( sigmaPointConstantIterator_->second - referenceVector ) *
                         ( sigmaPointConstantIterator_->second - referenceVector ).transpose( ) << std::endl << std::endl;
            weightedAverageMatrix += covarianceEstimationWeights_.at( sigmaPointConstantIterator_->first ) *
                    ( sigmaPointConstantIterator_->second - referenceVector ) *
                    ( sigmaPointConstantIterator_->second - referenceVector ).transpose( );
        }
        std::cout << weightedAverageMatrix << std::endl;
    }

    //! System function input by user.
    SystemFunction inputSystemFunction_;

    //! Measurement function input by user.
    MeasurementFunction inputMeasurementFunction_;

    //! Dimension of state vector.
    unsigned int stateDimension_;

    //! Dimension of measurement vector.
    unsigned int measurementDimension_;

    //! Dimension of augmented state vector.
    /*!
     *
     */
    unsigned int augmentedStateDimension_;

    //! Number of sigma points.
    /*!
     *
     */
    unsigned int numberOfSigmaPoints_;

    //! Value of alpha parameter.
    /*!
     *
     */
    std::vector< DependentVariableType > constantParameters_ = std::vector< DependentVariableType >( 5, 0.0 );

    //!
    /*!
     *
     */
    std::vector< DependentVariableType > stateEstimationWeights_;

    //!
    /*!
     *
     */
    std::vector< DependentVariableType > covarianceEstimationWeights_;

    //!
    /*!
     *
     */
    DependentVector augmentedStateVector_;

    //!
    /*!
     *
     */
    DependentMatrix augmentedCovarianceMatrix_;

    //!
    /*!
     *
     */
    std::map< unsigned int, DependentVector > mapOfSigmaPoints_;

    //!
    /*!
     *
     */
    typename std::map< unsigned int, DependentVector >::const_iterator sigmaPointConstantIterator_;

    //!
    /*!
     *
     */
    unsigned int currentSigmaPoint_;

};

//! Typedef for a filter with double data type.
typedef UnscentedKalmanFilter< > UnscentedKalmanFilterDouble;

//! Typedef for a shared-pointer to a filter with double data type.
typedef boost::shared_ptr< UnscentedKalmanFilterDouble > UnscentedKalmanFilterDoublePointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_UNSCENTED_KALMAN_FILTER_H
