#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"

namespace tudat
{

namespace basic_astrodynamics
{

Eigen::Vector3d updateAndGetTorque(
        const std::shared_ptr< TorqueModel > torqueModel,
        const double currentTime )
{
    // Update members.
    torqueModel->updateMembers( currentTime );

    // Evaluate and return acceleration.
    return torqueModel->getTorque( );
}

}

}
