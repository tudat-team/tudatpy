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

#ifndef TUDAT_LINEAR_KALMAN_FILTER_H
#define TUDAT_LINEAR_KALMAN_FILTER_H

#include "Tudat/Mathematics/Filters/kalmanFilter.h"

namespace tudat
{

namespace filters
{

//! Linear Kalman filter.
/*!
 *
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class LinearKalmanFilter: public KalmanFilterCore< IndependentVariableType, DependentVariableType >
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     */
    LinearKalmanFilter( const DependentMatrix& systemMatrix,
                        const DependentMatrix& inputMatrix,
                        const DependentMatrix& measurementMatrix,
                        const DependentMatrix& systemUncertainty,
                        const DependentMatrix& measurementUncertainty,
                        const bool isStateToBeIntegrated = false,
                        const IntegratorPointer integrator = NULL ) :
        KalmanFilterCore( isStateToBeIntegrated, integrator ),
        systemMatrix_( systemMatrix ), inputMatrix_( inputMatrix ), measurementMatrix_( measurementMatrix ),
        systemUncertainty_( systemUncertainty ), measurementUncertainty_( measurementUncertainty )
    {
        // Create system and measurement functions based on input parameters
        systemFunction_ = boost::bind( &createSystemFunction, _1, _2, _3 );
        measurementFunction_ = boost::bind( &createMeasurementFunction, _1, _2 );
    }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~LinearKalmanFilter( ){ }

    //!
    DependentVector updateFilter( const DependentVector& currentStateVector, const DependentVector& currentControlVector )
    {

    }

private:

    //! System function.
    DependentVector createSystemFunction( const IndependentVariableType independentVariable,
                                          const DependentVector& stateVector,
                                          const DependentVector& controlVector )
    {
        return systemMatrix_ * stateVector + inputMatrix_ * controlVector + this->produceSystemNoise( );
    }

    //! Measurement function.
    DependentVector createMeasurementFunction( const IndependentVariableType independentVariable,
                                               const DependentVector& stateVector )
    {
        return measurementMatrix_ * stateVector + this->produceMeasurementNoise( );
    }

    //!
    DependentMatrix systemMatrix_;

    //!
    DependentMatrix inputMatrix_;

    //!
    DependentMatrix measurementMatrix_;

};

} // namespace filters

} // namespace tudat

#endif // TUDAT_LINEAR_KALMAN_FILTER_H
