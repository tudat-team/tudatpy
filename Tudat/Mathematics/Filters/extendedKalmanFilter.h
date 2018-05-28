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

//! Extended Kalman filter.
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
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::SystemFunction SystemFunction;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::MeasurementFunction MeasurementFunction;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::SystemMatrixFunction SystemMatrixFunction;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::MeasurementMatrixFunction MeasurementMatrixFunction;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::IntegratorSettings IntegratorSettings;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::Integrator Integrator;

    //! Default constructor.
    /*!
     *  Default constructor. This constructor takes the state and measurement functions and their respective
     *  Jacobian functions as inputs. These functions can be a function of time, state and (for state) control vector.
     *  \param systemFunction Function returning the state as a function of time, state and control input. Can be a differential
     *      equation if the isStateToBeIntegrated boolean is set to true.
     *  \param measurementFunction Function returning the measurement as a function of time and state.
     *  \param stateJacobianFunction Function returning the Jacobian of the system w.r.t. the state. The input values can
     *      be time, state and control input.
     *  \param stateNoiseJacobianFunction Function returning the Jacobian of the system function w.r.t. the system noise. The input
     *      values can be time, state and control input.
     *  \param measurementJacobianFunction Function returning the Jacobian of the measurement function w.r.t. the state. The input
     *      values can be time and state.
     *  \param measurementNoiseJacobianFunction Function returning the Jacobian of the measurement function w.r.t. the measurement
     *      noise. The input values can be time and state.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param isStateToBeIntegrated Boolean defining whether the system function needs to be integrated.
     *  \param integrator Pointer to integrator to be used to propagate state.
     */
    ExtendedKalmanFilter( const SystemFunction& systemFunction,
                          const MeasurementFunction& measurementFunction,
                          const SystemMatrixFunction& stateJacobianFunction,
                          const SystemMatrixFunction& stateNoiseJacobianFunction,
                          const MeasurementMatrixFunction& measurementJacobianFunction,
                          const MeasurementMatrixFunction& measurementNoiseJacobianFunction,
                          const DependentMatrix& systemUncertainty,
                          const DependentMatrix& measurementUncertainty,
                          const IndependentVariableType initialTime,
                          const DependentVector& initialStateVector,
                          const DependentMatrix& initialCovarianceMatrix,
                          const bool isStateToBeIntegrated = false,
                          const boost::shared_ptr< IntegratorSettings > integratorSettings = NULL ) :
        KalmanFilterBase< IndependentVariableType, DependentVariableType >( systemUncertainty, measurementUncertainty, initialTime,
                                                                            initialStateVector, initialCovarianceMatrix,
                                                                            isStateToBeIntegrated, integratorSettings ),
        inputSystemFunction_( systemFunction ), inputMeasurementFunction_( measurementFunction ),
        stateJacobianFunction_( stateJacobianFunction ), stateNoiseJacobianFunction_( stateNoiseJacobianFunction ),
        measurementJacobianFunction_( measurementJacobianFunction ),
        measurementNoiseJacobianFunction_( measurementNoiseJacobianFunction )
    {
        // Compute the discrete-time version of the system Jacobians
        if ( isStateToBeIntegrated )
        {
            discreteTimeStateJacobians_ = boost::bind( &ExtendedKalmanFilter< IndependentVariableType,
                                                       DependentVariableType >::generateDiscreteTimeSystemJacobians,
                                                       this, _1, _2, _3 );
        }
    }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~ExtendedKalmanFilter( ){ }

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
        // Prediction step
        DependentVector aPrioriStateEstimate;
        DependentMatrix currentStateJacobianMatrix;
        DependentMatrix currentStateNoiseJacobianMatrix;
        if ( this->isStateToBeIntegrated_ )
        {
            aPrioriStateEstimate = this->integrateState( currentTime, currentControlVector );
            std::pair< DependentMatrix, DependentMatrix > discreteTimeJacobians =
                    discreteTimeStateJacobians_( currentTime, aPrioriStateEstimate, currentControlVector );
            currentStateJacobianMatrix = discreteTimeJacobians.first;
            currentStateNoiseJacobianMatrix = discreteTimeJacobians.second;
        }
        else
        {
            aPrioriStateEstimate = this->systemFunction_( currentTime, this->aPosterioriStateEstimate_,
                                                          currentControlVector );
            currentStateJacobianMatrix = stateJacobianFunction_( currentTime, this->aPosterioriStateEstimate_,
                                                                 currentControlVector );
            currentStateNoiseJacobianMatrix = stateNoiseJacobianFunction_( currentTime, this->aPosterioriStateEstimate_,
                                                                           currentControlVector );
        }
        DependentVector measurmentEstimate = this->measurementFunction_( currentTime, aPrioriStateEstimate );

        // Compute remaining Jacobians
        DependentMatrix currentMeasurementJacobianMatrix = measurementJacobianFunction_( currentTime, aPrioriStateEstimate );
        DependentMatrix currentMeasurementNoiseJacobianMatrix = measurementNoiseJacobianFunction_( currentTime,
                                                                                                   aPrioriStateEstimate );

        // Prediction step (continued)
        DependentMatrix aPrioriCovarianceEstimate = currentStateJacobianMatrix * this->aPosterioriCovarianceEstimate_ *
                currentStateJacobianMatrix.transpose( ) + currentStateNoiseJacobianMatrix * this->systemUncertainty_ *
                currentStateNoiseJacobianMatrix.transpose( );

        // Compute Kalman gain
        DependentMatrix kalmanGain = aPrioriCovarianceEstimate * currentMeasurementJacobianMatrix.transpose( ) * (
                    currentMeasurementJacobianMatrix * aPrioriCovarianceEstimate * currentMeasurementJacobianMatrix.transpose( ) +
                    currentMeasurementNoiseJacobianMatrix * this->measurementUncertainty_ *
                    currentMeasurementNoiseJacobianMatrix.transpose( ) ).inverse( );

        // Update step
        this->updateStateAndCovariance( currentTime, aPrioriStateEstimate, aPrioriCovarianceEstimate,
                                        currentMeasurementJacobianMatrix, currentMeasurementVector, measurmentEstimate,
                                        kalmanGain );
    }
//    ( const IndependentVariableType currentTime,
//    const DependentVector& aPrioriStateEstimate,
//    const DependentMatrix& aPrioriCovarianceEstimate,
//    const DependentMatrix& currentMeasurementMatrix,
//    const DependentVector& currentMeasurementVector,
//    const DependentVector& measurmentEstimate,
//    const DependentMatrix& kalmanGain )

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
        return inputSystemFunction_( currentTime, currentStateVector, currentControlVector );
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
        return inputMeasurementFunction_( currentTime, currentStateVector );
    }

    //! Function to generate the discrete-time version of the system Jacobians.
    /*!
     *  Function to generate the discrete-time version of the system Jacobians, from the continuous-time
     *  versions. The transformation is carried out by using the exponential of a matrix, which is expanded
     *  into [1]:
     * \f[
     *  \exp{ A * t } = I + \sum_{ k = 1 }^{ \infty } \frac{ A^k * t^k }{ k! }
     * \f]
     *  where only the first three terms of the expansion are used.
     *  \param currentTime Scalar representing the current time.
     *  \param currentStateVector Vector representing the current state.
     *  \param currentControlVector Vector representing the current control input.
     */
    std::pair< DependentMatrix, DependentMatrix > generateDiscreteTimeSystemJacobians(
            const IndependentVariableType currentTime,
            const DependentVector& currentStateVector,
            const DependentVector& currentControlVector )
    {
        // Pre-compute Jacobians
        DependentMatrix stateJacobian = stateJacobianFunction_( currentTime, currentStateVector, currentControlVector );
        DependentMatrix noiseJacobian = stateNoiseJacobianFunction_( currentTime, currentStateVector, currentControlVector );

        // Get sizes
        unsigned int stateRows = stateJacobian.rows( );
        unsigned int stateCols = stateJacobian.cols( );
        unsigned int noiseRows = noiseJacobian.rows( );
        unsigned int noiseCols = noiseJacobian.cols( );

        // Merge Jacobians in one matrix
        DependentMatrix continuousJacobians = DependentMatrix::Zero( stateRows + noiseCols, stateCols + noiseCols );
        continuousJacobians.block( 0, 0, stateJacobian.rows( ), stateJacobian.cols( ) ) = stateJacobian;
        continuousJacobians.block( 0, stateCols + noiseCols - 1, noiseRows, noiseCols ) = noiseJacobian;

        // Generate discrete-time Jacobians
        DependentMatrix identityMatrix = DependentMatrix::Identity( stateRows + noiseCols, stateCols + noiseCols );
        DependentMatrix discreteJacobians = identityMatrix + continuousJacobians * this->integrationStepSize_ +
                continuousJacobians * continuousJacobians * std::pow( this->integrationStepSize_, 2 ) / 2.0 +
                continuousJacobians * continuousJacobians * continuousJacobians *
                std::pow( this->integrationStepSize_, 3 ) / 6.0;

        // Extract state and noise Jacobians
        return std::make_pair( discreteJacobians.block( 0, 0, stateRows, stateCols ),
                               discreteJacobians.block( 0, stateCols + noiseCols - 1, noiseRows, noiseCols ) );
    }

    //! System function input by user.
    SystemFunction inputSystemFunction_;

    //! Measurement function input by user.
    MeasurementFunction inputMeasurementFunction_;

    //! State Jacobian matrix function.
    /*!
     *  Function returning the continuous-time Jacobian of the systemFunction_ w.r.t. the state vector. The input values
     *  can be time, state and control input.
     */
    SystemMatrixFunction stateJacobianFunction_;

    //! State noise Jacobian matrix function.
    /*!
     *  Function returning the continuous-time Jacobian of the systemFunction_ w.r.t. the state noise. The input values
     *  can be time, state and control input.
     */
    SystemMatrixFunction stateNoiseJacobianFunction_;

    //! Measurement Jacobian matrix function.
    /*!
     *  Function returning the continuous-time Jacobian of the measurementFunction_ w.r.t. the state vector.
     *  The input values can be time, state and control input.
     */
    MeasurementMatrixFunction measurementJacobianFunction_;

    //! Measurement noise Jacobian matrix function.
    /*!
     *  Function returning the continuous-time Jacobian of the measurementFunction_ w.r.t. the measurement noise.
     *  The input values can be time, state and control input.
     */
    MeasurementMatrixFunction measurementNoiseJacobianFunction_;

    //! Function to return the discrete-time state Jacobians.
    /*!
     *  Function to return the discrete-time state Jacobians, i.e., the Jacobians of the stateFunction_ w.r.t. the
     *  state itself and the state noise. Transformation to discrete-time is done with the function
     *  generateDiscreteTimeSystemJacobians.
     */
    boost::function< std::pair< DependentMatrix, DependentMatrix >(
            const IndependentVariableType, const DependentVector&, const DependentVector& ) > discreteTimeStateJacobians_;

};

//! Typedef for a filter with double data type.
typedef ExtendedKalmanFilter< > ExtendedKalmanFilterDouble;

//! Typedef for a shared-pointer to a filter with double data type.
typedef boost::shared_ptr< ExtendedKalmanFilterDouble > ExtendedKalmanFilterDoublePointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_EXTENDED_KALMAN_FILTER_H
