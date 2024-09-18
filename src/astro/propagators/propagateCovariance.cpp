
#include<tudat/astro/propagators/propagateCovariance.h>

namespace tudat
{

namespace propagators
{

//! Function to retrueve full state transition and sensitivity matrices at epochs
void getFullVariationalEquationsSolutionHistory(
        std::map< double, Eigen::MatrixXd >& fullVariationalEquationsSolutionHistory,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    for( unsigned int i = 0; i < evaluationTimes.size( ); i++ )
    {
        fullVariationalEquationsSolutionHistory[ evaluationTimes.at( i ) ] =
                stateTransitionInterface->getFullCombinedStateTransitionAndSensitivityMatrix( evaluationTimes.at( i ), false  );
    }
}

//! Function to retrueve full state transition and sensitivity matrices at epochs
void getFullVariationalEquationsSolutionHistory(
        std::map< double, Eigen::MatrixXd >& fullVariationalEquationsSolutionHistory,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime )
{
    std::vector< double > evaluationTimes;
    double currentTime = initialTime;
    while( currentTime < finalTime )
    {
        evaluationTimes.push_back( currentTime );
        currentTime += timeStep;
    }

    getFullVariationalEquationsSolutionHistory(
                fullVariationalEquationsSolutionHistory, stateTransitionInterface, evaluationTimes );
}

//! Function to propagate full covariance at the initial time to state covariance at later times
void propagateCovariance(
        std::map< double, Eigen::MatrixXd >& propagatedCovariance,
        const Eigen::MatrixXd& initialCovariance,
        const std::map< double, Eigen::MatrixXd >& fullVariationalEquationsSolutionHistory )
{
    for( auto resultIterator : fullVariationalEquationsSolutionHistory )
    {
        propagatedCovariance[ resultIterator.first ] =
                resultIterator.second * initialCovariance *
                resultIterator.second.transpose( );
    }
}

//! Function to propagate full covariance at the initial time to state covariance at later times
void propagateCovariance(
        std::map< double, Eigen::MatrixXd >& propagatedCovariance,
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    if( initialCovariance.rows( ) != stateTransitionInterface->getFullParameterVectorSize( ) )
    {
        throw std::runtime_error( "Error when propagating single-arc covariance, sizes are incompatible" );
    }

    std::map< double, Eigen::MatrixXd > fullVariationalEquationsSolutionHistory;
    getFullVariationalEquationsSolutionHistory(
                fullVariationalEquationsSolutionHistory, stateTransitionInterface, evaluationTimes );
    propagateCovariance( propagatedCovariance, initialCovariance, fullVariationalEquationsSolutionHistory );
}

std::map< double, Eigen::MatrixXd > propagateCovariance(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    propagateCovariance(
                propagatedCovariance, initialCovariance, stateTransitionInterface, evaluationTimes );
    return propagatedCovariance;
}

//! Function to propagate full covariance at the initial time to state covariance at later times
void propagateCovariance(
        std::map< double, Eigen::MatrixXd >& propagatedCovariance,
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime )
{
    std::vector< double > evaluationTimes;
    double currentTime = initialTime;
    while( currentTime < finalTime )
    {
        evaluationTimes.push_back( currentTime );
        currentTime += timeStep;
    }

    propagateCovariance( propagatedCovariance, initialCovariance, stateTransitionInterface, evaluationTimes );

}

Eigen::MatrixXd convertCovarianceToFrame(
        const Eigen::MatrixXd inputCovariance,
        const Eigen::VectorXd inertialCartesianRelativeState,
        const reference_frames::SatelliteReferenceFrames inputFrame,
        const reference_frames::SatelliteReferenceFrames outputFrame )
{
    Eigen::MatrixXd outputCovariance;
    if( inertialCartesianRelativeState.rows( ) % 6 != 0 )
    {
        throw std::runtime_error( "Error when converting Cartesian state covariance to alternative frame; input size is not 6N" );
    }
    else if( ( inputCovariance.rows( ) != inertialCartesianRelativeState.rows( ) ) ||
             ( inputCovariance.cols( ) != inertialCartesianRelativeState.rows( ) ) )
    {
        throw std::runtime_error( "Error when converting Cartesian state covariance to alternative frame, state and covariance size don't match" );
    }
    else
    {
        int numberOfBodes = inertialCartesianRelativeState.rows( ) / 6;
        Eigen::MatrixXd covarianceTransformation =
                Eigen::MatrixXd::Zero( 6 * numberOfBodes, 6 * numberOfBodes );
        for( int i = 0; i < numberOfBodes; i++ )
        {
            Eigen::Matrix3d rotationMatrix = reference_frames::getRotationBetweenSatelliteFrames(
                        inertialCartesianRelativeState.segment( i * 6, 6 ), inputFrame, outputFrame );
            covarianceTransformation.block( 6 * i, 6 * i, 3, 3 ) = rotationMatrix;
            covarianceTransformation.block( 6 * i + 3, 6 * i + 3, 3, 3 ) = rotationMatrix;
        }
        outputCovariance = covarianceTransformation * inputCovariance * covarianceTransformation.transpose( );

    }
    return outputCovariance;
}

std::map< double, Eigen::MatrixXd > convertCovarianceHistoryToFrame(
        const std::map< double, Eigen::MatrixXd > inputCovariances,
        const std::map< double, Eigen::VectorXd > inertialCartesianRelativeStates,
        const reference_frames::SatelliteReferenceFrames inputFrame,
        const reference_frames::SatelliteReferenceFrames outputFrame )
{
    std::map< double, Eigen::MatrixXd > outputCovariances;

    auto covarianceIterator = inputCovariances.begin( );
    auto stateIterator = inertialCartesianRelativeStates.begin( );

    while( covarianceIterator != inputCovariances.end( ) )
    {
        if( stateIterator == inertialCartesianRelativeStates.end( ) )
        {
            throw std::runtime_error( "Error when converting covariance history, state history has ended before covariance history" );
        }
        if( covarianceIterator->first != stateIterator->first )
        {
            throw std::runtime_error( "Error when converting covariance history, input times are not consistent for state and covariance" );
        }
        outputCovariances[ covarianceIterator->first ] = convertCovarianceToFrame(
                    covarianceIterator->second, stateIterator->second,
                    inputFrame, outputFrame );
        covarianceIterator++;
        stateIterator++;
    }

    return outputCovariances;
}


void convertCovarianceHistoryToFormalErrorHistory(
        std::map< double, Eigen::VectorXd >& propagatedFormalErrors,
        std::map< double, Eigen::MatrixXd >& propagatedCovariance )
{
    propagatedFormalErrors.clear( );
    for( auto covarianceIterator : propagatedCovariance )
    {
        propagatedFormalErrors[ covarianceIterator.first ] =
                Eigen::VectorXd( covarianceIterator.second.diagonal( ).array( ).sqrt( ) );
    }
}

//! Function to propagate full covariance at the initial time to state formal errors at later times
void propagateFormalErrors(
        std::map< double, Eigen::VectorXd >& propagatedFormalErrors,
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    propagateCovariance(
                propagatedCovariance, initialCovariance, stateTransitionInterface, evaluationTimes );
    convertCovarianceHistoryToFormalErrorHistory( propagatedFormalErrors, propagatedCovariance );
}

//! Function to propagate full covariance at the initial time to state formal errors at later times
std::map< double, Eigen::VectorXd > propagateFormalErrors(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::VectorXd > propagatedFormalErrors;
    propagateFormalErrors( propagatedFormalErrors, initialCovariance, stateTransitionInterface, evaluationTimes );
    return propagatedFormalErrors;
}

//! Function to propagate full covariance at the initial time to state formal errors at later times
void propagateFormalErrors(
        std::map< double, Eigen::VectorXd >& propagatedFormalErrors,
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    propagateCovariance(
                propagatedCovariance, initialCovariance, stateTransitionInterface, timeStep, initialTime, finalTime );

    for( auto covarianceIterator : propagatedCovariance )
    {
        propagatedFormalErrors[ covarianceIterator.first ] =
                Eigen::VectorXd( covarianceIterator.second.diagonal( ).array( ).sqrt( ) );
    }
}

} // namespace propagators

} // namespace tudat

