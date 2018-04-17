
#include "Tudat/SimulationSetup/EstimationSetup/createTorquePartials.h"

namespace tudat
{

namespace simulation_setup
{

boost::shared_ptr< acceleration_partials::TorquePartial > createInertialTorquePartial(
        const std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > acceleratedBody )
{
    boost::function< Eigen::Vector3d( ) > angularVelocityFunction =
            boost::bind( &Body::getCurrentAngularVelocityVectorInLocalFrame, acceleratedBody.second );
    boost::function< Eigen::Matrix3d( ) > inertiaTensorFunction =
            boost::bind( &Body::getBodyInertiaTensor, acceleratedBody.second );

    return boost::make_shared< acceleration_partials::InertialTorquePartial >(
                angularVelocityFunction, inertiaTensorFunction, acceleratedBody.first );
}

} // namespace simulation_setup

} // namespace tudat

