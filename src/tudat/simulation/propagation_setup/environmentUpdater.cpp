#include "tudat/simulation/propagation_setup/environmentUpdater.h"

namespace tudat
{

namespace propagators
{

// template class EnvironmentUpdater< double, double >;

std::vector< StateDerivativeDependency > getTorqueStateDerivativeDependencies( const basic_astrodynamics::AvailableTorque torqueType )
{
    std::vector< StateDerivativeDependency > dependencyList;
    switch( torqueType )
    {
        case basic_astrodynamics::inertial_torque: {
            dependencyList.push_back( inertia_tensor_derivative_dependency );
            break;
        }
        default: {
            break;
        }
    }
    return dependencyList;
}

std::vector< StateDerivativeDependency > getGravityStateDerivativeDependencies(
        const basic_astrodynamics::GravityDeformationType gravityDeformationType )
{
    std::vector< StateDerivativeDependency > dependencyList;
    switch( gravityDeformationType )
    {
        case basic_astrodynamics::maxwell_deformation: {
            dependencyList.push_back( rotation_rate_derivative_dependency );
            break;
        }
        default: {
            break;
        }
    }
    return dependencyList;
}

IntegratedStateType getStateTypeForDependency( StateDerivativeDependency dependency )
{
    switch( dependency )
    {
        case inertia_tensor_derivative_dependency:
            return gravity_deformation_state;
        case rotation_rate_derivative_dependency:
            return rotational_state;
        default:
            throw std::runtime_error( "state derivative dependency not recognised." );
    }
}

}  // namespace propagators

}  // namespace tudat
