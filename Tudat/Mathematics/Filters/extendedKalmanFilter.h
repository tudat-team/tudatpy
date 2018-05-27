/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
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
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::IntegratorSettings IntegratorSettings;
    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::Integrator Integrator;

    //! Typedefs for matrix functions.
    typedef boost::function< DependentMatrix( const IndependentVariableType, const DependentVector&, const DependentVector& ) > SystemMatrixFunction;
    typedef boost::function< DependentMatrix( const IndependentVariableType, const DependentVector& ) > MeasurementMatrixFunction;

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
        stateJacobianFunction_( stateJacobianFunction ), stateNoiseJacobianFunction_( stateNoiseJacobianFunction ),
        measurementJacobianFunction_( measurementJacobianFunction ), measurementNoiseJacobianFunction_( measurementNoiseJacobianFunction )
    {
        // Create system and measurement functions based on input values
        this->systemFunction_ = systemFunction;
        this->measurementFunction_ = measurementFunction;

        // Create numerical integrator
        if ( this->isStateToBeIntegrated_ )
        {
            this->generateNumericalIntegrator( integratorSettings, initialStateVector );
        }
    }

    //! Constructor.
    /*!
     *  Constructor taking the constant state and measurement jacobian matrices as inputs.
     *  \param systemFunction
     *  \param measurementFunction
     *  \param stateJacobianMatrix
     *  \param stateNoiseJacobianMatrix
     *  \param measurementJacobianMatrix
     *  \param measurementNoiseJacobianMatrix
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param isStateToBeIntegrated Boolean defining whether the system function needs to be integrated.
     *  \param integrator Pointer to integrator to be used to propagate state.
     *  \paragraph integrationStepSize Step size for integration.
     */
    ExtendedKalmanFilter( const SystemFunction& systemFunction,
                          const MeasurementFunction& measurementFunction,
                          const DependentMatrix& stateJacobianMatrix,
                          const DependentMatrix& stateNoiseJacobianMatrix,
                          const DependentMatrix& measurementJacobianMatrix,
                          const DependentMatrix& measurementNoiseJacobianMatrix,
                          const DependentMatrix& systemUncertainty,
                          const DependentMatrix& measurementUncertainty,
                          const IndependentVariableType initialTime,
                          const DependentVector& initialStateVector,
                          const DependentMatrix& initialCovarianceMatrix,
                          const bool isStateToBeIntegrated = false,
                          const boost::shared_ptr< IntegratorSettings > integratorSettings = NULL ) :
        ExtendedKalmanFilter( systemFunction, measurementFunction, systemUncertainty, measurementUncertainty,
                              boost::lambda::constant( stateJacobianMatrix ), boost::lambda::constant( stateNoiseJacobianMatrix ),
                              boost::lambda::constant( measurementJacobianMatrix ), boost::lambda::constant( measurementNoiseJacobianMatrix ),
                              initialTime, initialStateVector, initialCovarianceMatrix, isStateToBeIntegrated, integratorSettings )
    { }

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
        // Compute variables for current step
        DependentMatrix currentStateJacobianMatrix = stateJacobianFunction_( currentTime, this->aPosterioriStateEstimate_, currentControlVector );
        DependentMatrix currentStateNoiseJacobianMatrix =
                stateNoiseJacobianFunction_( currentTime, this->aPosterioriStateEstimate_, currentControlVector );
        DependentMatrix currentMeasurementJacobianMatrix = measurementJacobianFunction_( currentTime, this->aPosterioriStateEstimate_ );
        DependentMatrix currentMeasurementNoiseJacobianMatrix =
                measurementNoiseJacobianFunction_( currentTime, this->aPosterioriStateEstimate_ );

        // Prediction step
        DependentVector aPrioriStateEstimate;
        if ( this->isStateToBeIntegrated_ )
        {
            aPrioriStateEstimate = this->integrateState( currentTime );
        }
        else
        {
            aPrioriStateEstimate = this->systemFunction_( currentTime, this->aPosterioriStateEstimate_, currentControlVector );
        }
        DependentVector measurmentEstimate = this->measurementFunction_( currentTime, aPrioriStateEstimate );
        DependentMatrix aPrioriCovarianceEstimate = currentStateJacobianMatrix * this->aPosterioriCovarianceEstimate_ *
                currentStateJacobianMatrix.transpose( ) + currentStateNoiseJacobianMatrix * this->systemUncertainty_ *
                currentStateNoiseJacobianMatrix.transpose( );

        // Compute Kalman gain
        DependentMatrix kalmanGain = aPrioriCovarianceEstimate * currentMeasurementJacobianMatrix.transpose( ) * (
                    currentMeasurementJacobianMatrix * aPrioriCovarianceEstimate * currentMeasurementJacobianMatrix.transpose( ) +
                    currentMeasurementNoiseJacobianMatrix * this->measurementUncertainty_ *
                    currentMeasurementNoiseJacobianMatrix.transpose( ) ).inverse( );

        // Update step
        this->updateStateAndCovariance( aPrioriStateEstimate, aPrioriCovarianceEstimate, currentMeasurementJacobianMatrix,
                                        currentMeasurementVector, measurmentEstimate, kalmanGain );
    }

private:

    //! State jacobian matrix function.
    /*!
     *
     */
    SystemMatrixFunction stateJacobianFunction_;

    //! State noise jacobian matrix function.
    /*!
     *
     */
    SystemMatrixFunction stateNoiseJacobianFunction_;

    //! Measurement jacobian matrix function.
    /*!
     *
     */
    MeasurementMatrixFunction measurementJacobianFunction_;

    //! Measurement noise jacobian matrix function.
    /*!
     *
     */
    MeasurementMatrixFunction measurementNoiseJacobianFunction_;

};

//! Typedef for a filter with double data type.
typedef ExtendedKalmanFilter< > ExtendedKalmanFilterDouble;

//! Typedef for a shared-pointer to a filter with double data type.
typedef boost::shared_ptr< ExtendedKalmanFilterDouble > ExtendedKalmanFilterDoublePointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_EXTENDED_KALMAN_FILTER_H
