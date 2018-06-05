/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      E. Mooij, AE4870B - Re-entry Systems, Lecture Notes, Delft University of Technology, 2016
 */

#ifndef TUDAT_LINEAR_KALMAN_FILTER_H
#define TUDAT_LINEAR_KALMAN_FILTER_H

#include "Tudat/Mathematics/Filters/kalmanFilter.h"

namespace tudat
{

namespace filters
{

//! Linear Kalman filter class.
/*!
 *  Class for the set up and use of the linear Kalman filter.
 *  \tparam IndependentVariableType Type of independent variable. Default is double.
 *  \tparam DependentVariableType Type of dependent variable. Default is double.
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class LinearKalmanFilter: public KalmanFilterBase< IndependentVariableType, DependentVariableType >
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
     *  Default constructor. This constructor takes state, control and measurement matrix functions as inputs.
     *  These functions can be a function of time, state and (for state and control matrices) control vector.
     *  \param stateTransitionMatrixFunction Function returning the state transition matrix, as a function
     *      of time, state and control input.
     *  \param controlMatrixFunction Function returning the control matrix as a function of time, state and control input.
     *  \param measurementMatrixFunction Function returning the measurement matrix as a function of time, state and
     *      control input.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param integrator Pointer to integrator to be used to propagate state.
     */
    LinearKalmanFilter( const MatrixFunction& stateTransitionMatrixFunction,
                        const MatrixFunction& controlMatrixFunction,
                        const MatrixFunction& measurementMatrixFunction,
                        const DependentMatrix& systemUncertainty,
                        const DependentMatrix& measurementUncertainty,
                        const IndependentVariableType initialTime,
                        const DependentVector& initialStateVector,
                        const DependentMatrix& initialCovarianceMatrix,
                        const boost::shared_ptr< IntegratorSettings > integratorSettings = NULL ) :
        KalmanFilterBase< IndependentVariableType, DependentVariableType >( systemUncertainty, measurementUncertainty,
                                                                            initialTime, initialStateVector,
                                                                            initialCovarianceMatrix, integratorSettings ),
        stateTransitionMatrixFunction_( stateTransitionMatrixFunction ), controlMatrixFunction_( controlMatrixFunction ),
        measurementMatrixFunction_( measurementMatrixFunction )
    {
        // Temporary block of integration
        if ( integratorSettings != NULL )
        {
            throw std::runtime_error( "Error in linear Kalman filter. Propagation of the state is "
                                      "not currently supported." );
        }
    }

    //! Constructor.
    /*!
     *  Constructor taking constant state, control and measurement matrices as inputs.
     *  \param stateTransitionMatrix Constant matrix representing the state transition matrix.
     *  \param controlMatrix Constant matrix representing the control matrix.
     *  \param measurementMatrix Constant matrix representing the measurement matrix.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param integratorSettings Pointer to integration settings defining the integrator to be used to propagate the state.
     */
    LinearKalmanFilter( const DependentMatrix& stateTransitionMatrix,
                        const DependentMatrix& controlMatrix,
                        const DependentMatrix& measurementMatrix,
                        const DependentMatrix& systemUncertainty,
                        const DependentMatrix& measurementUncertainty,
                        const IndependentVariableType initialTime,
                        const DependentVector& initialStateVector,
                        const DependentMatrix& initialCovarianceMatrix,
                        const boost::shared_ptr< IntegratorSettings > integratorSettings = NULL ) :
        LinearKalmanFilter( boost::lambda::constant( stateTransitionMatrix ),
                            boost::lambda::constant( controlMatrix ),
                            boost::lambda::constant( measurementMatrix ),
                            systemUncertainty, measurementUncertainty, initialTime, initialStateVector,
                            initialCovarianceMatrix, integratorSettings )
    { }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~LinearKalmanFilter( ){ }

    //! Function to update the filter with the new step data.
    /*!
     *  Function to update the filter with the new step data.
     *  \param currentTime Scalar representing current time.
     *  \param currentMeasurementVector Vector representing current measurement.
     */
    void updateFilter( const IndependentVariableType currentTime, const DependentVector& currentMeasurementVector )
    {
        // Compute variables for current step
        DependentMatrix currentSystemMatrix = stateTransitionMatrixFunction_( currentTime, this->aPosterioriStateEstimate_ );
        DependentMatrix currentMeasurementMatrix = measurementMatrixFunction_( currentTime, this->aPosterioriStateEstimate_ );

        // Prediction step
        DependentVector aPrioriStateEstimate = this->predictState( currentTime );
        DependentVector measurmentEstimate = this->measurementFunction_( currentTime, aPrioriStateEstimate );
        DependentMatrix aPrioriCovarianceEstimate = currentSystemMatrix * this->aPosterioriCovarianceEstimate_ *
                currentSystemMatrix.transpose( ) + this->systemUncertainty_;

        // Compute Kalman gain
        DependentMatrix kalmanGain = aPrioriCovarianceEstimate * currentMeasurementMatrix.transpose( ) * (
                    currentMeasurementMatrix * aPrioriCovarianceEstimate * currentMeasurementMatrix.transpose( ) +
                    this->measurementUncertainty_ ).inverse( );

        // Correction step
        this->correctState( currentTime, aPrioriStateEstimate, currentMeasurementVector, measurmentEstimate, kalmanGain );
        this->correctCovariance( currentTime, aPrioriCovarianceEstimate, currentMeasurementMatrix, kalmanGain );
    }

private:

    //! Function to create the function that defines the system model.
    /*!
     *  Function to create the function that defines the system model. The output of this function is then bound
     *  to the systemFunction_ variable, via the boost::bind command.
     *  \param currentTime Scalar representing the current time.
     *  \param currentStateVector Vector representing the current state.
     *  \return Vector representing the estimated state.
     */
    DependentVector createSystemFunction( const IndependentVariableType currentTime,
                                          const DependentVector& currentStateVector )
    {
        return stateTransitionMatrixFunction_( currentTime, currentStateVector ) * currentStateVector;
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
        return measurementMatrixFunction_( currentTime, currentStateVector ) * currentStateVector;
    }

    //! State transition matrix function.
    /*!
     *  State transition matrix function, which will be used to generate the system function, together with the controlMatrixFunction_,
     *  by multiplying the matrices with the state and control vectors, respectively.
     */
    MatrixFunction stateTransitionMatrixFunction_;

    //! Control matrix function.
    /*!
     *  Control matrix function, which will be used to generate the system function, together with the stateTransitionMatrixFunction_,
     *  by multiplying the matrices with the control and state vectors, respectively.
     */
    MatrixFunction controlMatrixFunction_;

    //! Measurement matrix function.
    /*!
     *  Measurement matrix function, which will be used to generate the measurement function, by multiplying the matrix with
     *  the state vector.
     */
    MatrixFunction measurementMatrixFunction_;

};

//! Typedef for a filter with double data type.
typedef LinearKalmanFilter< > LinearKalmanFilterDouble;

//! Typedef for a shared-pointer to a filter with double data type.
typedef boost::shared_ptr< LinearKalmanFilterDouble > LinearKalmanFilterDoublePointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_LINEAR_KALMAN_FILTER_H
