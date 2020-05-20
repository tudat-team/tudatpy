/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/simulation/simulation.h>

#include "tudat/io/applicationOutput.h"

// Constant parameters for example
const double gravitationalParameter = 32.2;

//! Function describing the system model.
Eigen::Vector3d stateFunction( const double time, const Eigen::Vector3d& state, const Eigen::Vector3d& control )
{
    TUDAT_UNUSED_PARAMETER( time );
    TUDAT_UNUSED_PARAMETER( control );
    Eigen::Vector3d stateDerivative = Eigen::Vector3d::Zero( );
    stateDerivative[ 0 ] = state[ 1 ];
    stateDerivative[ 1 ] = 0.0034 * gravitationalParameter * std::exp( - state[ 0 ] / 22000.0 ) *
            std::pow( state[ 1 ], 2 ) / ( 2.0 * state[ 2 ] ) - gravitationalParameter;
    return stateDerivative;
}

//! Function describing the measurement model.
Eigen::Vector1d measurementFunction( const double time, const Eigen::Vector3d& state )
{
    TUDAT_UNUSED_PARAMETER( time );
    Eigen::Vector1d measurement;
    measurement[ 0 ] = state[ 0 ];
    return measurement;
}

//! Function producing the system Jacobian.
Eigen::Matrix3d stateJacobianFunction( const double time, const Eigen::Vector3d& state, const Eigen::Vector3d& control )
{
    TUDAT_UNUSED_PARAMETER( time );
    TUDAT_UNUSED_PARAMETER( control );
    Eigen::Matrix3d stateJacobian = Eigen::Matrix3d::Zero( );
    stateJacobian( 0, 1 ) = 1.0;
    stateJacobian( 1, 0 ) = - 0.0034 * gravitationalParameter * std::exp( - state[ 0 ] / 22000.0 ) *
            std::pow( state[ 1 ], 2 ) / ( 44000.0 * state[ 2 ] );
    stateJacobian( 1, 1 ) = 0.0034 * gravitationalParameter * std::exp( - state[ 0 ] / 22000.0 ) *
            state[ 1 ] / state[ 2 ];
    stateJacobian( 1, 2 ) = - 0.0034 * gravitationalParameter * std::exp( - state[ 0 ] / 22000.0 ) *
            std::pow( state[ 1 ], 2 ) / ( 2.0 * std::pow( state[ 2 ], 2 ) );
    return stateJacobian;
}

//! Function producing the measurement Jacobian.
Eigen::RowVector3d measurementJacobianFunction( const double time, const Eigen::Vector3d& state )
{
    TUDAT_UNUSED_PARAMETER( time );
    TUDAT_UNUSED_PARAMETER( state );
    Eigen::RowVector3d measurementJacobian = Eigen::RowVector3d::Zero( );
    measurementJacobian[ 0 ] = 1.0;
    return measurementJacobian;
}

//! Class for control vector.
template< typename IndependentVariableType, typename DependentVariableType, int NumberOfElements >
class ControlSystem
{
public:

    //! Typedef of the control vector.
    typedef Eigen::Matrix< DependentVariableType, NumberOfElements, 1 > DependentVector;

    //! Typedef of the function describing the system.
    typedef std::function< DependentVector( const IndependentVariableType,
                                            const DependentVector& ) > ControlFunction;

    //! Default constructor.
    /*!
     *  Default constructor.
     *  \param controlFunction Function to compute the control vector.
     */
    ControlSystem( const ControlFunction& controlFunction ) :
        controlFunction_( controlFunction )
    {
        // Set control vector to zero
        currentControlVector_.setZero( );
    }

    //! Default destructor.
    ~ControlSystem( ) { }

    //! Function to retireve the current control vector.
    /*!
     *  Function to retireve the current control vector. The function setCurrentControlVector needs to be called before the control
     *  vector is retrieved.
     *  \return Vector denoting the current control.
     */
    DependentVector getCurrentControlVector( )
    {
        return currentControlVector_;
    }

    //! Function to set the control vector.
    /*!
     *  Function to set the control vector, based on the current time and state. The control vector is then
     *  computed based on the input control function.
     *  \param currentTime Double denoting the current time.
     *  \param currentStateVector Vector denoting the current state.
     */
    void setCurrentControlVector( const IndependentVariableType currentTime,
                                  const DependentVector& currentStateVector )
    {
        currentControlVector_ = controlFunction_( currentTime, currentStateVector );
    }

private:

    //! Function to compute the control vector.
    ControlFunction controlFunction_;

    //! Vector denoting the current control.
    DependentVector currentControlVector_;

};

//! Execute examples on tabulated atmosphere.
int main( )
{
    using namespace tudat;
    using namespace tudat::filters;
    using namespace tudat_applications;

    // Set initial conditions
    const double initialTime = 0;
    const double timeStepSize = 0.1;
    const unsigned int numberOfTimeSteps = 300;

    Eigen::Vector3d initialStateVector;
    initialStateVector[ 0 ] = 200000.0;
    initialStateVector[ 1 ] = -6000.0;
    initialStateVector[ 2 ] = 500.0;

    Eigen::Vector3d initialEstimatedStateVector;
    initialEstimatedStateVector[ 0 ] = 200025.0;
    initialEstimatedStateVector[ 1 ] = -6150.0;
    initialEstimatedStateVector[ 2 ] = 800.0;

    Eigen::Matrix3d initialEstimatedStateCovarianceMatrix = Eigen::Matrix3d::Zero( );
    initialEstimatedStateCovarianceMatrix( 0, 0 ) = std::pow( 1000.0, 2 );
    initialEstimatedStateCovarianceMatrix( 1, 1 ) = 20000.0;
    initialEstimatedStateCovarianceMatrix( 2, 2 ) = std::pow( 300.0, 2 );

    // Set system and measurement uncertainty
    Eigen::Matrix3d systemUncertainty = Eigen::Matrix3d::Zero( );
    Eigen::Vector1d measurementUncertainty = Eigen::Vector1d::Zero( );
    systemUncertainty( 0, 0 ) = std::pow( 100.0, 2 );
    systemUncertainty( 1, 1 ) = std::pow( 10.0, 2 );
    systemUncertainty( 2, 2 ) = std::pow( 1.0, 2 );
    measurementUncertainty[ 0 ] = std::pow( 25.0, 2 );

    // Set integrator settings
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< > > (
                numerical_integrators::euler, initialTime, timeStepSize );

    // Create control classes
    // These are only included to show how a control system would need to be implemented, but they have no effect at all on the
    // results of this example.
    std::shared_ptr< ControlSystem< double, double, 3 > > extendedControl =
            std::make_shared< ControlSystem< double, double, 3 > >(
                [ ]( const double, const Eigen::Vector3d& ){ return Eigen::Vector3d::Zero( ); } );
    std::shared_ptr< ControlSystem< double, double, 3 > > unscentedControl =
            std::make_shared< ControlSystem< double, double, 3 > >(
                [ ]( const double, const Eigen::Vector3d& ){ return Eigen::Vector3d::Zero( ); } );

    // Create filter settings
    std::shared_ptr< FilterSettings< double, double > > extendedFilterSettings =
            std::make_shared< ExtendedKalmanFilterSettings< double, double > >(
                systemUncertainty,
                measurementUncertainty,
                timeStepSize,
                initialTime,
                initialEstimatedStateVector,
                initialEstimatedStateCovarianceMatrix,
                integratorSettings );
    std::shared_ptr< FilterSettings< double, double > > unscentedFilterSettings =
            std::make_shared< UnscentedKalmanFilterSettings< double, double > >(
                systemUncertainty,
                measurementUncertainty,
                timeStepSize,
                initialTime,
                initialEstimatedStateVector,
                initialEstimatedStateCovarianceMatrix,
                integratorSettings );

    // Create filter objects
    std::shared_ptr< FilterBase< double, double > > extendedFilter = createFilter< double, double >(
                extendedFilterSettings,
                std::bind( &stateFunction, std::placeholders::_1, std::placeholders::_2,
                           std::bind( &ControlSystem< double, double, 3 >::getCurrentControlVector, extendedControl ) ),
                std::bind( &measurementFunction, std::placeholders::_1, std::placeholders::_2 ),
                std::bind( &stateJacobianFunction, std::placeholders::_1, std::placeholders::_2,
                           std::bind( &ControlSystem< double, double, 3 >::getCurrentControlVector, extendedControl ) ),
                ( [ ]( double, const Eigen::VectorXd& ){ return  Eigen::Matrix3d::Identity( ); } ),
                std::bind( &measurementJacobianFunction, std::placeholders::_1, std::placeholders::_2 ),
                ( [ ]( double, const Eigen::VectorXd& ){ return Eigen::Vector1d::Identity( ); } ) );

    std::shared_ptr< FilterBase< double, double > > unscentedFilter = createFilter< double, double >(
                unscentedFilterSettings,
                std::bind( &stateFunction, std::placeholders::_1, std::placeholders::_2,
                           std::bind( &ControlSystem< double, double, 3 >::getCurrentControlVector, unscentedControl ) ),
                std::bind( &measurementFunction, std::placeholders::_1, std::placeholders::_2 ) );

    // Loop over each time step
    const bool showProgress = false;
    double currentTime = extendedFilter->getCurrentTime( );;
    Eigen::Vector3d currentActualStateVector = initialStateVector;
    Eigen::Vector3d currentNoisyStateVector;
    Eigen::Vector3d currentControlVector = Eigen::Vector3d::Zero( );
    Eigen::Vector1d currentMeasurementVector;
    std::map< double, Eigen::Vector3d > actualStateVectorHistory;
    std::map< double, Eigen::Vector1d > measurementVectorHistory;
    actualStateVectorHistory[ initialTime ] = initialStateVector;
    for ( unsigned int i = 0; i < numberOfTimeSteps; i++ )
    {
        // Compute actual values and perturb them
        currentActualStateVector += ( stateFunction( currentTime, currentActualStateVector, currentControlVector ) +
                                      unscentedFilter->produceSystemNoise( ) ) * timeStepSize;
        currentMeasurementVector = measurementFunction( currentTime, currentActualStateVector ) +
                unscentedFilter->produceMeasurementNoise( );

        // Update control classes
        extendedControl->setCurrentControlVector( currentTime, extendedFilter->getCurrentStateEstimate( ) );
        unscentedControl->setCurrentControlVector( currentTime, unscentedFilter->getCurrentStateEstimate( ) );

        // Update filters
        extendedFilter->updateFilter( currentMeasurementVector );
        unscentedFilter->updateFilter( currentMeasurementVector );

        // Update time
        currentTime = extendedFilter->getCurrentTime( );

        // Store values
        actualStateVectorHistory[ currentTime ] = currentActualStateVector;
        measurementVectorHistory[ currentTime ] = currentMeasurementVector;

        // Print progress
        if ( showProgress )
        {
            std::cout << "Time: " << currentTime << std::endl
                      << "Measurement: " << currentMeasurementVector.transpose( ) << std::endl
                      << "EKF Estimated State: " << extendedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl
                      << "UKF Estimated State: " << unscentedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl << std::endl;
        }
    }

    // Save actual state history
    input_output::writeDataMapToTextFile( actualStateVectorHistory, "actualStateHistory.dat", getOutputPath( "FilterEstimation" ) );

    // Extract and save Kalman filter state histories
    input_output::writeDataMapToTextFile( extendedFilter->getEstimatedStateHistory( ),
                                          "EKFEstimatedStateHistory.dat", getOutputPath( "FilterEstimation" ) );
    input_output::writeDataMapToTextFile( unscentedFilter->getEstimatedStateHistory( ),
                                          "UKFEstimatedStateHistory.dat", getOutputPath( "FilterEstimation" ) );

    // Extract and save Kalman filter covariance histories
    input_output::writeDataMapToTextFile( extendedFilter->getEstimatedCovarianceHistory( ),
                                          "EKFEstimatedCovarianceHistory.dat", getOutputPath( "FilterEstimation" ) );
    input_output::writeDataMapToTextFile( unscentedFilter->getEstimatedCovarianceHistory( ),
                                          "UKFEstimatedCovarianceHistory.dat", getOutputPath( "FilterEstimation" ) );

    // Extract and save noise history
    std::pair< std::vector< Eigen::VectorXd >, std::vector< Eigen::VectorXd > > noiseHistory = unscentedFilter->getNoiseHistory( );
    Eigen::MatrixXd systemNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.first );
    Eigen::MatrixXd measurementNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.second );
    input_output::writeMatrixToFile( systemNoise, "systemNoise.dat", 16, getOutputPath( "FilterEstimation" ) );
    input_output::writeMatrixToFile( measurementNoise, "measurementNoise.dat", 16, getOutputPath( "FilterEstimation" ) );

    return EXIT_SUCCESS;
}
