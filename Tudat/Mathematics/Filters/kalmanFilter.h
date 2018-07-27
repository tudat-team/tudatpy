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

#ifndef TUDAT_KALMAN_FILTER_H
#define TUDAT_KALMAN_FILTER_H

#include "Tudat/Mathematics/Filters/filter.h"

namespace tudat
{

namespace filters
{

//! Kalman filter class.
/*!
 *  Base class for the set up and use of Kalman filters.
 *  \tparam IndependentVariableType Type of independent variable. Default is double.
 *  \tparam DependentVariableType Type of dependent variable. Default is double.
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class KalmanFilterBase: public FilterBase< IndependentVariableType, DependentVariableType >
{
public:

    //! Inherit typedefs from base class.
    typedef typename FilterBase< IndependentVariableType, DependentVariableType >::DependentVector DependentVector;
    typedef typename FilterBase< IndependentVariableType, DependentVariableType >::DependentMatrix DependentMatrix;
    typedef typename FilterBase< IndependentVariableType, DependentVariableType >::Function Function;
    typedef typename FilterBase< IndependentVariableType, DependentVariableType >::MatrixFunction MatrixFunction;
    typedef typename FilterBase< IndependentVariableType, DependentVariableType >::IntegratorSettings IntegratorSettings;
    typedef typename FilterBase< IndependentVariableType, DependentVariableType >::Integrator Integrator;

    //! Constructor.
    /*!
     *  Constructor.
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
    KalmanFilterBase( const DependentMatrix& systemUncertainty,
                      const DependentMatrix& measurementUncertainty,
                      const IndependentVariableType initialTime,
                      const DependentVector& initialStateVector,
                      const DependentMatrix& initialCovarianceMatrix,
                      const boost::shared_ptr< IntegratorSettings > integratorSettings ) :
        FilterBase< IndependentVariableType, DependentVariableType >( systemUncertainty, measurementUncertainty,
                                                                      initialTime, initialStateVector,
                                                                      initialCovarianceMatrix, integratorSettings )
    { }

    //! Destructor.
    virtual ~KalmanFilterBase( ){ }

protected:

    //! Function to predict the state for the next time step.
    /*!
     *  Function to predict the state for the next time step, with the either the use of the integrator provided in
     *  the integratorSettings, or the systemFunction_ input by the user.
     *  \param currentTime Scalar representing the current time.
     *  \param currentControlVector Vector representing the current control input.
     *  \return Propagated state at the requested time.
     */
    DependentVector predictState( const IndependentVariableType currentTime )
    {
        return predictState( currentTime, this->aPosterioriStateEstimate_ );
    }

    //! Function to predict the state for the next time step, by overwriting previous state.
    /*!
     *  Function to predict the state for the next time step, by overwriting previous state, with the either the use of
     *  the integrator provided in the integratorSettings, or the systemFunction_ input by the user.
     *  \param currentTime Scalar representing the current time.
     *  \param currentStateVector Vector representing the current state (which overwrites the previous state).
     *  \return Propagated state at the requested time.
     */
    DependentVector predictState( const IndependentVariableType currentTime,
                                  const DependentVector& currentStateVector )
    {
        return this->isStateToBeIntegrated_ ? propagateState( currentTime, currentStateVector ) :
                                              this->systemFunction_( currentTime, currentStateVector );
    }

    //! Function to correct the covariance for the next time step.
    /*!
     *  Function to predict the state for the next time step, by overwriting previous state, with the either the use of
     *  the integrator provided in the integratorSettings, or the systemFunction_ input by the user.
     *  \param currentTime Scalar representing the current time.
     */
    virtual void correctCovariance( const IndependentVariableType currentTime, const DependentMatrix& aPrioriCovarianceEstimate,
                                    const DependentMatrix& currentMeasurementMatrix, const DependentMatrix& kalmanGain )
    {
        this->aPosterioriCovarianceEstimate_ = ( this->identityMatrix_ - kalmanGain * currentMeasurementMatrix ) *
                aPrioriCovarianceEstimate * ( this->identityMatrix_ - kalmanGain * currentMeasurementMatrix ).transpose( ) +
                kalmanGain * this->measurementUncertainty_ * kalmanGain.transpose( );
        this->historyOfCovarianceEstimates_[ currentTime ] = this->aPosterioriCovarianceEstimate_;
    }

private:

    //! Function to propagate state to the next time step, by overwriting previous state.
    /*!
     *  Function to propagate state to the next time step, by overwriting previous state, with the use of the integrator
     *  provided in the integratorSettings.
     *  \param currentTime Scalar representing the current time.
     *  \param currentStateVector Vector representing the current state (which overwrites the previous state).
     *  \return Propagated state at the requested time.
     */
    DependentVector propagateState( const IndependentVariableType currentTime,
                                    const DependentVector& currentStateVector )
    {
        // Reset time and state
        this->integrator_->modifyCurrentState( currentStateVector, currentTime );

        // Integrate equations
        return this->integrator_->performIntegrationStep( this->integrationStepSize_ );
    }

};

//! Typedef for a filter with double data type.
typedef KalmanFilterBase< > KalmanFilterDouble;

//! Typedef for a shared-pointer to a filter with double data type.
typedef boost::shared_ptr< KalmanFilterDouble > KalmanFilterDoublePointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_KALMAN_FILTER_H
