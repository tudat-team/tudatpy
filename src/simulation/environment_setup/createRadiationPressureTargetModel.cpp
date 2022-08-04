#include "tudat/simulation/environment_setup/createRadiationPressureTargetModel.h"


namespace tudat
{
namespace simulation_setup
{

std::shared_ptr<electromagnetism::RadiationPressureTargetModel> createRadiationPressureTargetModel(
        std::shared_ptr<RadiationPressureTargetModelSettings> modelSettings, const std::string &body)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<electromagnetism::RadiationPressureTargetModel> radiationPressureTargetModel;

    switch(modelSettings->getRadiationPressureTargetModelType())
    {
        case cannonball_target:
        {
            auto cannonballTargetModelSettings =
                    std::dynamic_pointer_cast< CannonballRadiationPressureTargetModelSettings >(modelSettings);

            if(cannonballTargetModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected cannonball radiation pressure target for body " + body );
            }

            radiationPressureTargetModel = std::make_shared<CannonballRadiationPressureTargetModel>(
                    cannonballTargetModelSettings->getArea(),
                    cannonballTargetModelSettings->getCoefficient()
            );
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize radiation pressure target model settings for " + body );
    }

    return radiationPressureTargetModel;
}
} // tudat
} // electromagnetism
