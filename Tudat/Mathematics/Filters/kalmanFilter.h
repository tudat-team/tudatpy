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

#include <limits>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"

namespace tudat
{

namespace filters
{

//! Kalman Filter class.
/*!
 *  Base class for the set up and use of Kalman filters.
 *  \tparam IndependentVariableType Type of independent variable. Default is double.
 *  \tparam DependentVariableType Type of dependent variable. Default is double.
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class KalmanFilterBase
{
public:

    //! Typedef of the state and measurement vectors.
    typedef Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 > DependentVector;

    //! Typedef of the state and measurement matrices.
    typedef Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic > DependentMatrix;

    //! Typedef of the function describing the system.
    typedef boost::function< DependentVector( const IndependentVariableType, const DependentVector&, const DependentVector& ) > SystemFunction;

    //! Typedef of the function describing the measurements.
    typedef boost::function< DependentVector( const IndependentVariableType, const DependentVector& ) > MeasurementFunction;

    //! Typedef of the integrator settings.
    typedef numerical_integrators::IntegratorSettings< IndependentVariableType > IntegratorSettings;

    //! Typedef of the integrator.
    typedef numerical_integrators::NumericalIntegrator< IndependentVariableType, DependentVector,
    DependentVector, IndependentVariableType > Integrator;

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
                      const bool isStateToBeIntegrated,
                      const boost::shared_ptr< IntegratorSettings > integratorSettings ) :
        systemUncertainty_( systemUncertainty ), measurementUncertainty_( measurementUncertainty ), initialTime_( initialTime ),
        aPosterioriStateEstimate_( initialStateVector ), aPosterioriCovarianceEstimate_( initialCovarianceMatrix ),
        isStateToBeIntegrated_( isStateToBeIntegrated )
    {
        // Check that uncertainty matrices are square
        if ( systemUncertainty_.rows( ) != systemUncertainty_.cols( ) )
        {
            throw std::runtime_error( "Error in setting up filter. The system uncertainty matrix has to be square." );
        }
        if ( measurementUncertainty_.rows( ) != measurementUncertainty_.cols( ) )
        {
            throw std::runtime_error( "Error in setting up filter. The measurement uncertainty matrix has to be square." );
        }

        // Create noise distributions
        generateNoiseDistributions( );

        // Generate identity matrix
        identityMatrix_ = DependentMatrix::Identity( systemUncertainty_.rows( ), systemUncertainty_.cols( ) );
    }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    virtual ~KalmanFilterBase( ){ }

    //! Function to update the filter with the data from the new time step.
    /*!
     *  Function to update the filter with the new step data.
     *  \param currentTime Scalar representing current time.
     *  \param currentControlVector Vector representing the current control input.
     *  \param currentMeasurementVector Vector representing current measurement.
     */
    virtual void updateFilter( const IndependentVariableType currentTime, const DependentVector& currentControlVector,
                               const DependentVector& currentMeasurementVector ) = 0;

    //! Function to update the state vector and covariance matrix with the data from the new time step.
    /*!
     *  Function to update the state vector and covariance matrix with the data from the new time step.
     *  \param currentTime Scalar representing current time.
     *  \param currentControlVector Vector representing the current control input.
     *  \param currentMeasurementVector Vector representing current measurement.
     */
    void updateStateAndCovariance( const DependentVector& aPrioriStateEstimate,
                                   const DependentMatrix& aPrioriCovarianceEstimate,
                                   const DependentMatrix& currentMeasurementMatrix,
                                   const DependentVector& currentMeasurementVector,
                                   const DependentVector& measurmentEstimate,
                                   const DependentMatrix& kalmanGain )
    {
        // Update step
        aPosterioriStateEstimate_ = aPrioriStateEstimate + kalmanGain * ( currentMeasurementVector - measurmentEstimate );
        aPosterioriCovarianceEstimate_ = ( identityMatrix_ - kalmanGain * currentMeasurementMatrix ) * aPrioriCovarianceEstimate;

        // Store values in history
        estimatedStateHistory_.push_back( aPosterioriStateEstimate_ );
        estimatedCovarianceHistory_.push_back( aPosterioriCovarianceEstimate_ );
    }

    //! Function to retrieve current state estimate.
    /*!
     *  Function to retrieve current state estimate. The state estimate needs to first be computed by updating the filter with the
     *  updateFilter function.
     *  \return Current state estimate.
     */
    Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 > getCurrentStateEstimate( )
    {
        return aPosterioriStateEstimate_;
    }

    //! Function to retrieve the history of estimated states.
    /*!
     *  Function to retrieve the history of estimated states.
     *  \return History of estimated states.
     */
    std::vector< DependentVector > getEstimatedStateHistory( )
    {
        return estimatedStateHistory_;
    }

    //! Function to retrieve the history of estimated covariance matrices.
    /*!
     *  Function to retrieve the history of estimated covariance matrices.
     *  \return History of estimated covariance matrices.
     */
    std::vector< DependentMatrix > getEstimatedCovarianceHistory( )
    {
        return estimatedCovarianceHistory_;
    }

    //! Function to retrieve the history of system and measurement noise used by the updateFilter function.
    /*!
     *  Function to retrieve the history of system and measurement noise.
     *  \return History of system and measurement noise, output as a std::pair.
     */
    std::pair< std::vector< DependentVector >, std::vector< DependentVector > > getNoiseHistory( )
    {
        return std::make_pair( systemNoiseHistory_, measurementNoiseHistory_ );
    }

protected:

    //! Function to generate the numerical integrator to be used for propagation of the state.
    /*!
     *
     */
    void generateNumericalIntegrator(
            const boost::shared_ptr< IntegratorSettings > integratorSettings,
            const DependentVector& initialStateVector )
    {
        // Get time step information
        integrationStepSize_ = integratorSettings->initialTimeStep_;

        // Generate integrator
        switch ( integratorSettings->integratorType_ )
        {
        case numerical_integrators::euler:
        {
            integrator_ = boost::make_shared< numerical_integrators::EulerIntegrator< IndependentVariableType, DependentVector,
                    DependentVector, IndependentVariableType > >( systemFunction_, initialTime_, initialStateVector );
            break;
        }
        default:
            throw std::runtime_error( "Error in Kalman filtering. Only Euler integration is currently supported." );
        }
    }

    //! Function to propagate state to the next time step.
    /*!
     *  Function to propagate state to the next time step, by using the integrator settings provided
     *  in the input.
     *  \param intervalEnd Time to propagate to.
     *  \param initialStepSize Initial (or constant) step size for propagation.
     *  \param finalTimeTolerance Tolerance used to define whether final time has been reached.
     *  \return Propagated state at the requested time.
     */
    DependentVector integrateState( const IndependentVariableType intervalEnd,
                                    const IndependentVariableType finalTimeTolerance = std::numeric_limits< IndependentVariableType >::epsilon( ) )
    {
        return integrator_->integrateTo( intervalEnd, integrationStepSize_, finalTimeTolerance );
    }

    //! Function to produce system noise.
    /*!
     *  Function to produce system noise, based on a Gaussian distribution, with zero mean and standard
     *  deviation given by the diagonal elements of the input system uncertainty matrix.
     *  \return Vector representing system noise.
     */
    DependentVector produceSystemNoise( )
    {
        // Declare system noise vector
        DependentVector systemNoise;
        systemNoise.resize( systemUncertainty_.rows( ), 1 );

        // Loop over dimensions and add noise
        for ( int i = 0; i < systemUncertainty_.rows( ); i++ )
        {
            if ( systemNoiseDistribution_.at( i ) != NULL )
            {
                systemNoise[ i ] = static_cast< DependentVariableType >( systemNoiseDistribution_.at( i )->getRandomVariableValue( ) );
            }
            else
            {
                systemNoise[ i ] = static_cast< DependentVariableType >( 0.0 );
            }
        }

        // Give back noise
        systemNoiseHistory_.push_back( systemNoise );
        return systemNoise;
    }

    //! Function to produce measurement noise.
    /*!
     *  Function to produce measurement noise, based on a Gaussian distribution, with zero mean and standard
     *  deviation given by the diagonal elements of the input measurement uncertainty matrix.
     *  \return Vector representing measurement noise.
     */
    DependentVector produceMeasurementNoise( )
    {
        // Declare system noise vector
        DependentVector measurementNoise;
        measurementNoise.resize( measurementUncertainty_.rows( ), 1 );

        // Loop over dimensions and add noise
        for ( int i = 0; i < measurementUncertainty_.rows( ); i++ )
        {
            if ( measurementNoiseDistribution_.at( i ) != NULL )
            {
                measurementNoise[ i ] = static_cast< DependentVariableType >( measurementNoiseDistribution_.at( i )->getRandomVariableValue( ) );
            }
            else
            {
                measurementNoise[ i ] = static_cast< DependentVariableType >( 0.0 );
            }
        }

        // Give back noise
        measurementNoiseHistory_.push_back( measurementNoise );
        return measurementNoise;
    }

    //! System function.
    /*!
     *  System function that will be used to retrieve the a-priori estimated state for the next step.
     */
    SystemFunction systemFunction_;

    //! Measurement function.
    /*!
     *  Measurement function that will be used to retrieve the estimated measurement for the next step,
     *  based on the current state.
     */
    MeasurementFunction measurementFunction_;

    //! Matrix representing the uncertainty in system modeling.
    DependentMatrix systemUncertainty_;

    //! Matrix representing the uncertainty in measurement modeling.
    DependentMatrix measurementUncertainty_;

    //! Scalar representing the initial time.
    IndependentVariableType initialTime_;

    //! Vector representing the a-posteriori estimated state.
    /*!
     *  Vector representing the a-posteriori estimated state, i.e., the state after the prediction and
     *  update steps of the Kalman filter.
     */
    DependentVector aPosterioriStateEstimate_;

    //! Matrix representing the a-posteriori estimated covariance.
    /*!
     *  Matrix representing the a-posteriori estimated covariance, i.e., the covariance after the prediction and
     *  update steps of the Kalman filter.
     */
    DependentMatrix aPosterioriCovarianceEstimate_;

    //! Boolean specifying whether the state needs to be integrated.
    bool isStateToBeIntegrated_;

    //! Vector of estimated states vectors.
    std::vector< DependentVector > estimatedStateHistory_;

    //! Vector of estimated covariance matrices.
    std::vector< DependentMatrix > estimatedCovarianceHistory_;

private:

    //! Function to generate the noise distributions for both system and measurement modeling.
    /*!
     *  Function to generate the noise distributions for both system and measurement modeling, which uses
     *  a Gaussian distribution, with zero mean and standard deviation given by the diagonal elements of the
     *  input system and measurement uncertainty matrices.
     */
    void generateNoiseDistributions( )
    {
        using namespace tudat::statistics;

        // Create system noise
        for ( int i = 0; i < systemUncertainty_.rows( ); i++ )
        {
            if ( static_cast< double >( systemUncertainty_( i, i ) ) != 0.0 )
            {
                systemNoiseDistribution_.push_back(
                            createBoostContinuousRandomVariableGenerator(
                                normal_boost_distribution, { 0.0, static_cast< double >(
                                                             std::sqrt( systemUncertainty_( i, i ) ) ) }, 12345.0 ) );
            }
            else
            {
                systemNoiseDistribution_.push_back( NULL );
            }
        }

        // Create measurement noise
        for ( int i = 0; i < measurementUncertainty_.rows( ); i++ )
        {
            if ( static_cast< double >( measurementUncertainty_( i, i ) ) != 0.0 )
            {
                measurementNoiseDistribution_.push_back(
                            createBoostContinuousRandomVariableGenerator(
                                normal_boost_distribution, { 0.0, static_cast< double >(
                                                             std::sqrt( measurementUncertainty_( i, i ) ) ) }, 54321.0 ) );
            }
            else
            {
                measurementNoiseDistribution_.push_back( NULL );
            }
        }
    }

    //! Indentity matrix.
    /*!
     *  Indentity matrix with the correct dimensions for the specific application.
     */
    DependentMatrix identityMatrix_;

    //! Vector where the system noise generators are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > systemNoiseDistribution_;

    //! Vector where the measurement noise generators are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > measurementNoiseDistribution_;

    //! Pointer to the integrator settings.
    /*!
     *  Pointer to the integrator settings, which is used to propagate the state to the new time step.
     */
    boost::shared_ptr< Integrator > integrator_;

    //! Scalar representing step size for integration.
    /*!
     *  Scalar representing step size for integration. If integrator_ points to a constant step size integrator, then
     *  this will be the constant step size, otherwise it will be the initial step size.
     */
    IndependentVariableType integrationStepSize_;

    //! Vector of system noise.
    std::vector< DependentVector > systemNoiseHistory_;

    //! Vector of measurement noise.
    std::vector< DependentVector > measurementNoiseHistory_;

};

//! Typedef for a filter with double data type.
typedef KalmanFilterBase< > KalmanFilterDouble;

//! Typedef for a shared-pointer to a filter with double data type.
typedef boost::shared_ptr< KalmanFilterDouble > KalmanFilterDoublePointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_KALMAN_FILTER_H
