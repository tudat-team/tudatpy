#include "tudat/simulation/environment_setup/createRadiationSourceModel.h"

namespace tudat
{
namespace simulation_setup
{

std::shared_ptr<electromagnetism::LuminosityModel>
createLuminosityModel(const std::shared_ptr<LuminosityModelSettings> modelSettings,
                        const std::string &body)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<electromagnetism::LuminosityModel> luminosityModel;

    switch(modelSettings->getLuminosityModelType())
    {
        case constant_radiant_power:
        {
            auto constantLuminosityModelSettings =
                    std::dynamic_pointer_cast< ConstantLuminosityModelSettings >(modelSettings);
            if(constantLuminosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected constant luminosity model for body " + body );
            }
            luminosityModel = std::make_shared<ConstantLuminosityModel>(
                constantLuminosityModelSettings->getLuminosity()
            );
            break;
        }
        case irradiance_based_radiant_power:
        {
            auto irradianceBasedLuminosityModelSettings =
                    std::dynamic_pointer_cast< IrradianceBasedLuminosityModelSettings >(modelSettings);
            if(irradianceBasedLuminosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected irradiance-based luminosity model for body " + body );
            }
            luminosityModel = std::make_shared<IrradianceBasedLuminosityModel>(
                irradianceBasedLuminosityModelSettings->getIrradianceAtDistanceFunction(),
                irradianceBasedLuminosityModelSettings->getDistance()
            );
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize radiation source model settings for " + body );
    }

    return luminosityModel;
}

std::shared_ptr<electromagnetism::RadiationSourceModel>
createRadiationSourceModel(
        const std::shared_ptr<RadiationSourceModelSettings> modelSettings,
        const std::string &body,
        const SystemOfBodies& bodies)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<electromagnetism::RadiationSourceModel> radiationSourceModel;

    switch(modelSettings->getRadiationSourceModelType())
    {
    case isotropic_point_source:
    {
        auto isotropicPointModelSettings =
                std::dynamic_pointer_cast< IsotropicPointRadiationSourceModelSettings >(modelSettings);

        if(isotropicPointModelSettings == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected isotropic point radiation source for body " + body );
        }
        if(isotropicPointModelSettings->getLuminosityModelSettings() == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected isotropic point radiation source to have a luminosity model for body " + body);
        }

        auto sourceBody = bodies.getBody(body);
        auto luminosityModel = createLuminosityModel(
                isotropicPointModelSettings->getLuminosityModelSettings(), body);

        radiationSourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);
        break;
    }
    default:
        throw std::runtime_error( "Error, do not recognize radiation source model settings for " + body );
    }

    return radiationSourceModel;
}

} // tudat
} // electromagnetism
