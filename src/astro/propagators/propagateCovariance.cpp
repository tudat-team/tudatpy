
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
                stateTransitionInterface->getFullCombinedStateTransitionAndSensitivityMatrix( evaluationTimes.at( i )  );
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

    for( auto covarianceIterator : propagatedCovariance )
    {
        propagatedFormalErrors[ covarianceIterator.first ] =
                Eigen::VectorXd( covarianceIterator.second.diagonal( ).array( ).sqrt( ) );
    }
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

