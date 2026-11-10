#include "tudat/astro/orbit_determination/rotational_dynamics_partials/constantTorquePartial.h"
#include "tudat/simulation/estimation_setup/createTorquePartials.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to retrieve the gravitational parameter function of a body with rigid body properties derived from its gravity field
std::function< double( ) > getGravitationalParameterFunctionOfGravityDerivedBody( const std::shared_ptr< Body > body )
{
    std::function< double( ) > gravitationalParameterFunction;
    if( std::dynamic_pointer_cast< FromGravityFieldRigidBodyProperties >( body->getMassProperties( ) ) != nullptr )
    {
        gravitationalParameterFunction =
                std::bind( &gravitation::GravityFieldModel::getGravitationalParameter, body->getGravityFieldModel( ) );
    }
    return gravitationalParameterFunction;
}

//! Function to create torque partial to be used for constant torques in angular acceleration
std::shared_ptr< acceleration_partials::TorquePartial > createConstantTorqueRotationalDynamicsPartial(
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const basic_astrodynamics::SingleBodyTorqueModelMap& torqueVector )
{
    std::function< Eigen::Vector3d( ) > angularVelocityFunction =
            std::bind( &Body::getCurrentAngularVelocityVectorInLocalFrame, acceleratedBody.second );
    std::function< Eigen::Matrix3d( ) > inertiaTensorFunction = std::bind( &Body::getBodyInertiaTensor, acceleratedBody.second );

    // The inertia tensor only depends on the gravity field if it is derived from it (gravitational parameter function is empty otherwise)
    std::function< double( ) > gravitationalParameterFunction =
            getGravitationalParameterFunctionOfGravityDerivedBody( acceleratedBody.second );

    std::function< double( ) > inertiaTensorNormalizationFunction;
    if( gravitationalParameterFunction != nullptr &&
        std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( acceleratedBody.second->getGravityFieldModel( ) ) !=
                nullptr )
    {
        inertiaTensorNormalizationFunction = std::bind( &gravitation::SphericalHarmonicsGravityField::getInertiaTensorNormalizationFactor,
                                                        std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                                                                acceleratedBody.second->getGravityFieldModel( ) ) );
    }

    return std::make_shared< acceleration_partials::ConstantTorquePartial >( angularVelocityFunction,
                                                                             inertiaTensorFunction,
                                                                             inertiaTensorNormalizationFunction,
                                                                             gravitationalParameterFunction,
                                                                             torqueVector,
                                                                             acceleratedBody.first );
}

}  // namespace simulation_setup

}  // namespace tudat
