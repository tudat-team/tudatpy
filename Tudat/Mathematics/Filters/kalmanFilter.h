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
#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"
#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"

namespace tudat
{

namespace filters
{

//! Kalman Filter class.
/*!
 *
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class KalmanFilterCore
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

    //! Typedef of the integrator.
    typedef numerical_integrators::NumericalIntegrator< IndependentVariableType, DependentVector > Integrator;

    //! Constructor.
    /*!
     *  Constructor.
     */
    KalmanFilterCore( const DependentMatrix& systemUncertainty,
                      const DependentMatrix& measurementUncertainty,
                      const bool isStateToBeIntegrated,
                      const boost::shared_ptr< Integrator > integrator ) :
        systemUncertainty_( systemUncertainty ), measurementUncertainty_( measurementUncertainty ),
      isStateToBeIntegrated_( isStateToBeIntegrated ), integrator_( integrator )
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
    }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    virtual ~KalmanFilterCore( ){ }

    //!
    /*!
     *
     */
    virtual Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 > updateFilter(
            const IndependentVariableType currentTime, const DependentVector& currentControlVector,
            const DependentVector& currentMeasurementVector ) = 0;

    //!
    /*!
     *
     */
    virtual Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 > getCurrentStateEstimate( ) = 0;

    //!
    /*!
     *
     */
    std::pair< std::vector< DependentVector >, std::vector< DependentVector > > getNoiseHistory( )
    {
        return std::make_pair( systemNoiseHistory_, measurementNoiseHistory_ );
    }

    //!
    /*!
     *
     */
    DependentVector integrateState( const IndependentVariableType intervalEnd,
                                    const IndependentVariableType initialStepSize,
                                    const IndependentVariableType finalTimeTolerance = std::numeric_limits< IndependentVariableType >::epsilon( ) )
    {
        return integrator_->integrateTo( intervalEnd, initialStepSize, finalTimeTolerance );
    }

protected:

    //!
    /*!
     *
     */
    virtual DependentVector createSystemFunction( const IndependentVariableType independentVariable,
                                                  const DependentVector& stateVector,
                                                  const DependentVector& controlVector ) = 0;

    //!
    /*!
     *
     */
    virtual DependentVector createMeasurementFunction( const IndependentVariableType independentVariable,
                                                       const DependentVector& stateVector ) = 0;

    //!
    /*!
     *
     */
    DependentVector produceSystemNoise( )
    {
        // Declare system noise vector
        DependentVector systemNoise;
        systemNoise.resize( systemUncertainty_.rows( ), 1 );

        // Loop over dimensions and add noise
        for ( int i = 0; i < systemUncertainty_.rows( ); i++ )
        {
            systemNoise[ i ] = static_cast< DependentVariableType >( systemNoiseDistribution_.at( i )->getRandomVariableValue( ) );
        }

        // Give back noise
        systemNoiseHistory_.push_back( systemNoise );
        return systemNoise;
    }

    //!
    /*!
     *
     */
    DependentVector produceMeasurementNoise( )
    {
        // Declare system noise vector
        DependentVector measurementNoise;
        measurementNoise.resize( measurementUncertainty_.rows( ), 1 );

        // Loop over dimensions and add noise
        for ( int i = 0; i < measurementUncertainty_.rows( ); i++ )
        {
            measurementNoise[ i ] = static_cast< DependentVariableType >( measurementNoiseDistribution_.at( i )->getRandomVariableValue( ) );
        }

        // Give back noise
        measurementNoiseHistory_.push_back( measurementNoise );
        return measurementNoise;
    }

    //!
    SystemFunction systemFunction_;

    //!
    MeasurementFunction measurementFunction_;

    //!
    DependentMatrix systemUncertainty_;

    //!
    DependentMatrix measurementUncertainty_;

    //!
    bool isStateToBeIntegrated_;

    //!
    boost::shared_ptr< Integrator > integrator_;

private:

    //!
    /*!
     *
     */
    void generateNoiseDistributions( )
    {
        using namespace tudat::statistics;

        // Create system noise
        for ( int i = 0; i < systemUncertainty_.rows( ); i++ )
        {
            systemNoiseDistribution_.push_back( createBoostContinuousRandomVariableGenerator(
                        normal_boost_distribution, { 0.0, static_cast< double >( systemUncertainty_( i, i ) ) }, 12345 ) );
        }

        // Create measurement noise
        for ( int i = 0; i < measurementUncertainty_.rows( ); i++ )
        {
            measurementNoiseDistribution_.push_back( createBoostContinuousRandomVariableGenerator(
                        normal_boost_distribution, { 0.0, static_cast< double >( measurementUncertainty_( i, i ) ) }, 12345 ) );
        }
    }

    //!
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > systemNoiseDistribution_;

    //!
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > measurementNoiseDistribution_;

    //!
    std::vector< DependentVector > systemNoiseHistory_;

    //!
    std::vector< DependentVector > measurementNoiseHistory_;

};

//! Typedef for a filter with double data type.
typedef KalmanFilterCore< > KalmanFilter;

//! Typedef for a shared-pointer to a filter with double data type.
typedef boost::shared_ptr< KalmanFilter > KalmanFilterPointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_KALMAN_FILTER_H
