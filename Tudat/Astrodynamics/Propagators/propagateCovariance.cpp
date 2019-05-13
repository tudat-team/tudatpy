
#include<Tudat/Astrodynamics/Propagators/propagateCovariance.h>

namespace tudat
{

namespace propagators
{


std::map< double, Eigen::MatrixXd > propagateCovariance(
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime )
{
    if( initialCovariance.rows( ) != stateTransitionInterface->getFullParameterVectorSize( ) )
    {
        throw std::runtime_error( "Error when propagating single-arc covariance, sizes are incompatible" );
    }

    std::map< double, Eigen::MatrixXd > parameterCovariance;

    Eigen::MatrixXd currentCombinedFullStateTransitionSensitivityMatrix;
    double currentTime = initialTime;
    while( currentTime < finalTime )
    {
        currentCombinedFullStateTransitionSensitivityMatrix =
                stateTransitionInterface->getFullCombinedStateTransitionAndSensitivityMatrix( currentTime );

        parameterCovariance[ currentTime ] =
                currentCombinedFullStateTransitionSensitivityMatrix * initialCovariance *
                currentCombinedFullStateTransitionSensitivityMatrix.transpose( );
        currentTime += timeStep;
    }

    return parameterCovariance;
}

std::map< double, Eigen::VectorXd > propagateFormalErrors(
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance = propagateCovariance(
                initialCovariance, stateTransitionInterface, timeStep, initialTime, finalTime );

    std::map< double, Eigen::VectorXd > propagatedFormalErrors;
    for( auto covarianceIterator : propagatedCovariance )
    {
        propagatedFormalErrors[ covarianceIterator.first ] =
            Eigen::VectorXd( covarianceIterator.second.diagonal( ).array( ).sqrt( ) );
    }
    return propagatedFormalErrors;
}

} // namespace propagators

} // namespace tudat

