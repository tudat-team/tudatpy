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
 *      Ogata, K., Discrete-Time Control Systems, 2nd ed. Pearson Education Asia, 2002.
 */

#ifndef TUDAT_EXTENDED_KALMAN_FILTER_H
#define TUDAT_EXTENDED_KALMAN_FILTER_H

#include "Tudat/Mathematics/Filters/kalmanFilter.h"

namespace tudat
{

namespace filters
{

//! Extended Kalman filter class.
/*!
 *  Class for the set up and use of the extended Kalman filter.
 *  \tparam IndependentVariableType Type of independent variable. Default is double.
 *  \tparam DependentVariableType Type of dependent variable. Default is double.
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class ExtendedKalmanFilter: public KalmanFilterBase< IndependentVariableType, DependentVariableType >
{
public:

    //! Inherit typedefs from base class.
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::DependentVector DependentVector;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::DependentMatrix DependentMatrix;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::Function Function;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::MatrixFunction MatrixFunction;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::IntegratorSettings IntegratorSettings;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::Integrator Integrator;

    //! Default constructor.
    /*!
     *  Default constructor. This constructor takes state and measurement functions and their respective
     *  Jacobian functions as inputs. These functions can be a function of time and state vector.
     *  \param systemFunction Function returning the state as a function of time and state vector. Can be a differential
     *      equation if the integratorSettings is set (i.e., if it is not a nullptr).
     *  \param measurementFunction Function returning the measurement as a function of time and state.
     *  \param stateJacobianFunction Function returning the Jacobian of the system w.r.t. the state. The input values can
     *      be time and state vector.
     *  \param stateNoiseJacobianFunction Function returning the Jacobian of the system function w.r.t. the system noise. The
     *      input values can be time and state vector.
     *  \param measurementJacobianFunction Function returning the Jacobian of the measurement function w.r.t. the state. The input
     *      values can be time and state vector.
     *  \param measurementNoiseJacobianFunction Function returning the Jacobian of the measurement function w.r.t. the measurement
     *      noise. The input values can be time and state vector.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param filteringStepSize Scalar representing the value of the constant filtering time step.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param integratorSettings Pointer to integration settings defining the integrator to be used to propagate the state.
     */
    ExtendedKalmanFilter( const Function& systemFunction,
                          const Function& measurementFunction,
                          const MatrixFunction& stateJacobianFunction,
                          const MatrixFunction& stateNoiseJacobianFunction,
                          const MatrixFunction& measurementJacobianFunction,
                          const MatrixFunction& measurementNoiseJacobianFunction,
                          const DependentMatrix& systemUncertainty,
                          const DependentMatrix& measurementUncertainty,
                          const IndependentVariableType filteringStepSize,
                          const IndependentVariableType initialTime,
                          const DependentVector& initialStateVector,
                          const DependentMatrix& initialCovarianceMatrix,
                          const std::shared_ptr< IntegratorSettings > integratorSettings = nullptr ) :
        KalmanFilterBase< IndependentVariableType, DependentVariableType >( systemUncertainty, measurementUncertainty,
                                                                            filteringStepSize, initialTime, initialStateVector,
                                                                            initialCovarianceMatrix, integratorSettings ),
        inputSystemFunction_( systemFunction ), inputMeasurementFunction_( measurementFunction ),
        stateJacobianFunction_( stateJacobianFunction ), stateNoiseJacobianFunction_( stateNoiseJacobianFunction ),
        measurementJacobianFunction_( measurementJacobianFunction ), measurementNoiseJacobianFunction_( measurementNoiseJacobianFunction )
    {
        // Compute the discrete-time version of the system Jacobians
        if ( this->isStateToBeIntegrated_ )
        {
            discreteTimeStateJacobians_ = std::bind( &ExtendedKalmanFilter< IndependentVariableType,
                                                     DependentVariableType >::generateDiscreteTimeSystemJacobians, this,
                                                     std::placeholders::_1 );
        }
    }

    //! Destructor.
    ~ExtendedKalmanFilter( ){ }

    //! Function to update the filter with the new step data.
    /*!
     *  Function to update the filter with the new step data.
     *  \param currentMeasurementVector Vector representing current measurement.
     */
	void updateFilter( const DependentVector& currentMeasurementVector )
    {
        // Prediction step
        DependentVector aPrioriStateEstimate = this->predictState( );
        DependentMatrix currentStateJacobianMatrix;
        DependentMatrix currentStateNoiseJacobianMatrix;
        if ( this->isStateToBeIntegrated_ )
        {
            std::pair< DependentMatrix, DependentMatrix > discreteTimeJacobians =
                    discreteTimeStateJacobians_( aPrioriStateEstimate );
            currentStateJacobianMatrix = discreteTimeJacobians.first;
            currentStateNoiseJacobianMatrix = discreteTimeJacobians.second;
        }
        else
        {
            currentStateJacobianMatrix = stateJacobianFunction_( this->currentTime_, aPrioriStateEstimate );
            currentStateNoiseJacobianMatrix = stateNoiseJacobianFunction_( this->currentTime_, aPrioriStateEstimate );
        }
        DependentVector measurementEstimate = this->measurementFunction_( this->currentTime_, aPrioriStateEstimate );

        // Compute remaining Jacobians
        DependentMatrix currentMeasurementJacobianMatrix = measurementJacobianFunction_( this->currentTime_, aPrioriStateEstimate );
        DependentMatrix currentMeasurementNoiseJacobianMatrix = measurementNoiseJacobianFunction_( this->currentTime_, aPrioriStateEstimate );

        // Prediction step (continued)
        DependentMatrix aPrioriCovarianceEstimate = currentStateJacobianMatrix * this->aPosterioriCovarianceEstimate_ *
                currentStateJacobianMatrix.transpose( ) + currentStateNoiseJacobianMatrix * this->systemUncertainty_ *
                currentStateNoiseJacobianMatrix.transpose( );

        // Compute Kalman gain
        DependentMatrix kalmanGain = aPrioriCovarianceEstimate * currentMeasurementJacobianMatrix.transpose( ) * (
                    currentMeasurementJacobianMatrix * aPrioriCovarianceEstimate * currentMeasurementJacobianMatrix.transpose( ) +
                    currentMeasurementNoiseJacobianMatrix * this->measurementUncertainty_ *
                    currentMeasurementNoiseJacobianMatrix.transpose( ) ).inverse( );

        // Correction step
        this->currentTime_ += this->filteringStepSize_;
        this->correctState( aPrioriStateEstimate, currentMeasurementVector, measurementEstimate, kalmanGain );
        this->correctCovariance( aPrioriCovarianceEstimate, currentMeasurementJacobianMatrix, kalmanGain );
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
        return inputSystemFunction_( currentTime, currentStateVector );
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
        return inputMeasurementFunction_( currentTime, currentStateVector );
    }

    //! Function to generate the discrete-time version of the system Jacobians.
    /*!
     *  Function to generate the discrete-time version of the system Jacobians, from the continuous-time
     *  versions. The transformation is carried out by using the matrix exponential.
     *  \param currentStateVector Vector representing the current state.
     *  \return Pair of discrete-time state and noise Jacobians.
     */
    std::pair< DependentMatrix, DependentMatrix > generateDiscreteTimeSystemJacobians(
            const DependentVector& currentStateVector )
    {
        // Pre-compute Jacobians
        DependentMatrix stateJacobian = stateJacobianFunction_( this->currentTime_, currentStateVector );
        DependentMatrix noiseJacobian = stateNoiseJacobianFunction_( this->currentTime_, currentStateVector );

        // Get sizes
        unsigned int stateDimension = stateJacobian.rows( );
        unsigned int noiseCols = noiseJacobian.cols( );

        // Merge Jacobians in one matrix
        DependentMatrix continuousJacobians = DependentMatrix::Zero( stateDimension + noiseCols, stateDimension + noiseCols );
        continuousJacobians.block( 0, 0, stateDimension, stateDimension ) = stateJacobian;
        continuousJacobians.block( 0, stateDimension, stateDimension, noiseCols ) = noiseJacobian;

        // Generate discrete-time Jacobians
        DependentMatrix discreteJacobians = ( continuousJacobians * this->filteringStepSize_ ).exp( );

        // Extract state and noise Jacobians
        return std::make_pair( discreteJacobians.block( 0, 0, stateDimension, stateDimension ),
                               discreteJacobians.block( 0, stateDimension, stateDimension, noiseCols ) );
    }

    //! System function input by user.
    Function inputSystemFunction_;

    //! Measurement function input by user.
    Function inputMeasurementFunction_;

    //! State Jacobian matrix function.
    /*!
     *  Function returning the continuous-time Jacobian of the systemFunction_ w.r.t. the state vector. The input values
     *  can be time and state vector.
     */
    MatrixFunction stateJacobianFunction_;

    //! State noise Jacobian matrix function.
    /*!
     *  Function returning the continuous-time Jacobian of the systemFunction_ w.r.t. the state noise. The input values
     *  can be time and state vector.
     */
    MatrixFunction stateNoiseJacobianFunction_;

    //! Measurement Jacobian matrix function.
    /*!
     *  Function returning the continuous-time Jacobian of the measurementFunction_ w.r.t. the state vector.
     *  The input values can be time and state vector.
     */
    MatrixFunction measurementJacobianFunction_;

    //! Measurement noise Jacobian matrix function.
    /*!
     *  Function returning the continuous-time Jacobian of the measurementFunction_ w.r.t. the measurement noise.
     *  The input values can be time and state vector.
     */
    MatrixFunction measurementNoiseJacobianFunction_;

    //! Function to return the discrete-time state Jacobians.
    /*!
     *  Function to return the discrete-time state Jacobians, i.e., the Jacobians of the stateFunction_ w.r.t. the
     *  state itself and the state noise. Transformation to discrete-time is done with the function
     *  generateDiscreteTimeSystemJacobians.
     */
    std::function< std::pair< DependentMatrix, DependentMatrix >( const DependentVector& ) > discreteTimeStateJacobians_;

};

//! Typedef for a filter with double data type.
typedef ExtendedKalmanFilter< > ExtendedKalmanFilterDouble;

//! Typedef for a shared-pointer to a filter with double data type.
typedef std::shared_ptr< ExtendedKalmanFilterDouble > ExtendedKalmanFilterDoublePointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_EXTENDED_KALMAN_FILTER_H
