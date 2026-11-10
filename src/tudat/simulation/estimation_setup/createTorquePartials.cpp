#include "tudat/astro/orbit_determination/rotational_dynamics_partials/constantTorquePartial.h"
#include "tudat/simulation/estimation_setup/createTorquePartials.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create torque partial to be used for constant torques in angular acceleration
std::shared_ptr< acceleration_partials::TorquePartial > createConstantTorqueRotationalDynamicsPartial(
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const basic_astrodynamics::SingleBodyTorqueModelMap& torqueVector )
{
    std::function< Eigen::Vector3d( ) > angularVelocityFunction =
            std::bind( &Body::getCurrentAngularVelocityVectorInLocalFrame, acceleratedBody.second );
    std::function< Eigen::Matrix3d( ) > inertiaTensorFunction = std::bind( &Body::getBodyInertiaTensor, acceleratedBody.second );

    return std::make_shared< acceleration_partials::ConstantTorquePartial >( angularVelocityFunction,
                                                                             inertiaTensorFunction,
                                                                             acceleratedBody.second->getMassProperties( ),
                                                                             torqueVector,
                                                                             acceleratedBody.first );
}

}  // namespace simulation_setup

}  // namespace tudat
