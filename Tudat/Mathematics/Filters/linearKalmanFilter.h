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

//! Linear Kalman filter.
/*!
 *
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class LinearKalmanFilterCore: public KalmanFilterCore< IndependentVariableType, DependentVariableType >
{
public:

    //! Inherit typedefs from base class.
    typedef typename KalmanFilterCore< IndependentVariableType, DependentVariableType >::DependentVector DependentVector;
    typedef typename KalmanFilterCore< IndependentVariableType, DependentVariableType >::DependentMatrix DependentMatrix;
    typedef typename KalmanFilterCore< IndependentVariableType, DependentVariableType >::SystemFunction SystemFunction;
    typedef typename KalmanFilterCore< IndependentVariableType, DependentVariableType >::MeasurementFunction MeasurementFunction;
    typedef typename KalmanFilterCore< IndependentVariableType, DependentVariableType >::Integrator Integrator;

    //! Typedefs for matrix functions.
    typedef boost::function< DependentMatrix( const IndependentVariableType, const DependentVector&, const DependentVector& ) > SystemMatrixFunction;
    typedef boost::function< DependentMatrix( const IndependentVariableType, const DependentVector& ) > MeasurementMatrixFunction;

    //! Default constructor.
    /*!
     *  Default constructor.
     */
    LinearKalmanFilterCore( const SystemMatrixFunction& systemMatrixFunction,
                            const SystemMatrixFunction& inputMatrixFunction,
                            const MeasurementMatrixFunction& measurementMatrixFunction,
                            const DependentMatrix& systemUncertainty,
                            const DependentMatrix& measurementUncertainty,
                            const IndependentVariableType initialTime,
                            const DependentVector& initialStateVector,
                            const bool isStateToBeIntegrated = false,
                            const boost::shared_ptr< Integrator > integrator = NULL ) :
        KalmanFilterCore< IndependentVariableType, DependentVariableType >( systemUncertainty, measurementUncertainty,
                                                                            isStateToBeIntegrated, integrator ),
        systemMatrixFunction_( systemMatrixFunction ), inputMatrixFunction_( inputMatrixFunction ),
        measurementMatrixFunction_( measurementMatrixFunction ), initialTime_( initialTime ),
        aPosterioriStateEstimate_( initialStateVector ), aPosterioriCovarianceEstimate_( systemUncertainty )
    {
        // Create system and measurement functions based on input parameters
        this->systemFunction_ = boost::bind( &LinearKalmanFilterCore< IndependentVariableType, DependentVariableType >::createSystemFunction,
                                             this, _1, _2, _3 );
        this->measurementFunction_ = boost::bind( &LinearKalmanFilterCore< IndependentVariableType, DependentVariableType >::createMeasurementFunction,
                                                  this, _1, _2 );

        // Generate identity matrix
        identityMatrix_ = DependentMatrix::Identity( systemUncertainty.rows( ), systemUncertainty.cols( ) );
    }

    //! Constructor.
    /*!
     *  Constructor.
     */
    LinearKalmanFilterCore( const DependentMatrix& systemMatrix,
                            const DependentMatrix& inputMatrix,
                            const DependentMatrix& measurementMatrix,
                            const DependentMatrix& systemUncertainty,
                            const DependentMatrix& measurementUncertainty,
                            const IndependentVariableType initialTime,
                            const DependentVector& initialStateVector,
                            const bool isStateToBeIntegrated = false,
                            const boost::shared_ptr< Integrator > integrator = NULL ) :
        LinearKalmanFilterCore( boost::lambda::constant( systemMatrix ),
                                boost::lambda::constant( inputMatrix ),
                                boost::lambda::constant( measurementMatrix ),
                                systemUncertainty, measurementUncertainty, initialTime, initialStateVector, isStateToBeIntegrated, integrator )
    { }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~LinearKalmanFilterCore( ){ }

    //!
    /*!
     *
     */
    Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 > updateFilter(
            const IndependentVariableType currentTime, const DependentVector& currentControlVector,
            const DependentVector& currentMeasurementVector )
    {
        // Compute variables for current step
        DependentMatrix currentSystemMatrix = systemMatrixFunction_( currentTime, aPosterioriStateEstimate_, currentControlVector );
        DependentMatrix currentMeasurementMatrix = measurementMatrixFunction_( currentTime, aPosterioriStateEstimate_ );
//        std::cout << std::endl << "State Matrix: " << std::endl << currentSystemMatrix << std::endl
//                  << "Measurement Matrix: " << std::endl << currentMeasurementMatrix << std::endl << std::endl;

        // Prediction step
        DependentVector aPrioriStateEstimate = this->systemFunction_( currentTime, aPosterioriStateEstimate_, currentControlVector );
        DependentVector measurmentEstimate = this->measurementFunction_( currentTime, aPosterioriStateEstimate_ );
        DependentMatrix aPrioriCovarianceEstimate = currentSystemMatrix * aPosterioriCovarianceEstimate_ * currentSystemMatrix.transpose( ) +
                this->systemUncertainty_;
//        std::cout << "State Estimate: " << aPrioriStateEstimate.transpose( ) << std::endl
//                  << "Measurement Estimate: " << measurmentEstimate.transpose( ) << std::endl
//                  << "Covariance Estimate: " << std::endl << aPrioriCovarianceEstimate << std::endl << std::endl;

        // Compute Kalman gain
        DependentMatrix kalmanGain = aPrioriCovarianceEstimate * currentMeasurementMatrix.transpose( ) * (
                    currentMeasurementMatrix * aPrioriCovarianceEstimate * currentMeasurementMatrix.transpose( ) +
                    this->measurementUncertainty_ ).inverse( );
//        std::cout << "Kalman Gain: " << std::endl << kalmanGain << std::endl << std::endl;

        // Update step
        aPosterioriStateEstimate_ = aPrioriStateEstimate + kalmanGain * ( currentMeasurementVector - measurmentEstimate );
        aPosterioriCovarianceEstimate_ = ( identityMatrix_ + kalmanGain * currentMeasurementMatrix ) * aPrioriCovarianceEstimate;
//        std::cout << "New State Estimate: " << aPosterioriStateEstimate_.transpose( ) << std::endl
//                  << "New Covariance Estimate: " << std::endl << aPosterioriCovarianceEstimate_ << std::endl << std::endl;

        // Give output
        return aPosterioriStateEstimate_;
    }

    //!
    /*!
     *
     */
    virtual Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 > getCurrentStateEstimate( )
    {
        return aPosterioriStateEstimate_;
    }

private:

    //! System function.
    DependentVector createSystemFunction( const IndependentVariableType independentVariable,
                                          const DependentVector& stateVector,
                                          const DependentVector& controlVector )
    {
        return systemMatrixFunction_( independentVariable, stateVector, controlVector ) * stateVector +
                inputMatrixFunction_( independentVariable, stateVector, controlVector ) * controlVector + this->produceSystemNoise( );
    }

    //! Measurement function.
    DependentVector createMeasurementFunction( const IndependentVariableType independentVariable,
                                               const DependentVector& stateVector )
    {
        return measurementMatrixFunction_( independentVariable, stateVector ) * stateVector + this->produceMeasurementNoise( );
    }

    //!
    SystemMatrixFunction systemMatrixFunction_;

    //!
    SystemMatrixFunction inputMatrixFunction_;

    //!
    MeasurementMatrixFunction measurementMatrixFunction_;

    //!
    IndependentVariableType initialTime_;

    //!
    DependentVector aPosterioriStateEstimate_;

    //!
    DependentMatrix aPosterioriCovarianceEstimate_;

    //!
    DependentMatrix identityMatrix_;

};

//! Typedef for a filter with double data type.
typedef LinearKalmanFilterCore< > LinearKalmanFilter;

//! Typedef for a shared-pointer to a filter with double data type.
typedef boost::shared_ptr< LinearKalmanFilter > LinearKalmanFilterPointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_LINEAR_KALMAN_FILTER_H
