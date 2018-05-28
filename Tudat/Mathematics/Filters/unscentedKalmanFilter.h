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

#ifndef TUDAT_UNSCENTED_KALMAN_FILTER_H
#define TUDAT_UNSCENTED_KALMAN_FILTER_H

#include "Tudat/Mathematics/Filters/kalmanFilter.h"

namespace tudat
{

namespace filters
{

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
//    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::SystemMatrixFunction SystemMatrixFunction;
//    typedef typename KalmanFilterBase< IndependentVariableType, DependentVariableType >::MeasurementMatrixFunction MeasurementMatrixFunction;
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
                           const boost::shared_ptr< IntegratorSettings > integratorSettings = NULL ) :
        KalmanFilterBase< IndependentVariableType, DependentVariableType >( systemUncertainty, measurementUncertainty,
                                                                            initialTime, initialStateVector,
                                                                            initialCovarianceMatrix, integratorSettings ),
        inputSystemFunction_( systemFunction ), inputMeasurementFunction_( measurementFunction )
    { }

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

    //! System function input by user.
    SystemFunction inputSystemFunction_;

    //! Measurement function input by user.
    MeasurementFunction inputMeasurementFunction_;

};

} // namespace filters

} // namespace tudat

#endif // TUDAT_UNSCENTED_KALMAN_FILTER_H
