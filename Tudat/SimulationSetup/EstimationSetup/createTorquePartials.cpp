#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/constantTorquePartial.h"
#include "Tudat/SimulationSetup/EstimationSetup/createTorquePartials.h"

namespace tudat
{

namespace simulation_setup
{

boost::shared_ptr< acceleration_partials::TorquePartial > createConstantTorqueRotationalDynamicsPartial(
        const std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const basic_astrodynamics::SingleBodyTorqueModelMap& torqueVector )
{
    boost::function< Eigen::Vector3d( ) > angularVelocityFunction =
            boost::bind( &Body::getCurrentAngularVelocityVectorInLocalFrame, acceleratedBody.second );
    boost::function< Eigen::Matrix3d( ) > inertiaTensorFunction =
            boost::bind( &Body::getBodyInertiaTensor, acceleratedBody.second );

    boost::function< double( ) > inertiaTensorNormalizationFunction;
    if( boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                acceleratedBody.second->getGravityFieldModel( ) ) != NULL )
    {
        inertiaTensorNormalizationFunction =
                boost::bind( &gravitation::SphericalHarmonicsGravityField::getInertiaTensorNormalizationFactor,
                             boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                                 acceleratedBody.second->getGravityFieldModel( ) ) );
    }

    boost::function< double( ) > gravitationalParameterFunction;
    if( acceleratedBody.second->getGravityFieldModel( ) != NULL )
    {
        gravitationalParameterFunction =
                boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                             acceleratedBody.second->getGravityFieldModel( ) );
    }

    return boost::make_shared< acceleration_partials::ConstantTorquePartial >(
                angularVelocityFunction, inertiaTensorFunction, inertiaTensorNormalizationFunction, gravitationalParameterFunction,
                torqueVector, acceleratedBody.first );
}

} // namespace simulation_setup

} // namespace tudat

