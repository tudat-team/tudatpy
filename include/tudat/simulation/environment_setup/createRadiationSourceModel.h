#ifndef TUDAT_CREATERADIATIONSOURCEMODEL_H
#define TUDAT_CREATERADIATIONSOURCEMODEL_H

#include <memory>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"

namespace tudat
{
namespace simulation_setup
{

enum LuminosityModelType
{
    constant_radiant_power,
    irradiance_based_radiant_power
};

class LuminosityModelSettings
{
public:
    explicit LuminosityModelSettings(
            const LuminosityModelType luminosityModelType) :
            luminosityModelType_(luminosityModelType) {}

    virtual ~LuminosityModelSettings() = default;

    LuminosityModelType getLuminosityModelType() const
    {
        return luminosityModelType_;
    }

private:
    LuminosityModelType luminosityModelType_;
};

class ConstantLuminosityModelSettings : public LuminosityModelSettings
{
public:
    explicit ConstantLuminosityModelSettings(
            const double luminosity) :
            LuminosityModelSettings(constant_radiant_power),
            luminosity_(luminosity) {}

    double getLuminosity() const
    {
        return luminosity_;
    }

private:
    double luminosity_;
};

class IrradianceBasedLuminosityModelSettings : public LuminosityModelSettings
{
public:
    explicit IrradianceBasedLuminosityModelSettings(
            const std::function<double()> irradianceAtDistanceFunction,
            const double distance) :
            LuminosityModelSettings(irradiance_based_radiant_power),
            irradianceAtDistanceFunction_(irradianceAtDistanceFunction),
            distance_(distance) {}

    explicit IrradianceBasedLuminosityModelSettings(
            const double irradianceAtDistance,
            const double distance) :
            LuminosityModelSettings(irradiance_based_radiant_power),
            irradianceAtDistanceFunction_([=] () { return distance; }),
            distance_(distance) {}

    const std::function<double()> getIrradianceAtDistanceFunction() const
    {
        return irradianceAtDistanceFunction_;
    }

    double getDistance() const
    {
        return distance_;
    }

private:
    std::function<double()> irradianceAtDistanceFunction_;
    double distance_;
};


//! List of radiation source models available in simulations
/*!
 *  List of radiation source models available in simulations. Radiation source models
 *  not defined by this given enum cannot be used for automatic model setup.
 */
enum RadiationSourceModelType
{
    isotropic_point_source,
    paneled_source
};

class RadiationSourceModelSettings
{
public:
    explicit RadiationSourceModelSettings(
            const RadiationSourceModelType radiationSourceModelType) :
            radiationSourceModelType_(radiationSourceModelType) {}

    virtual ~RadiationSourceModelSettings() = default;

    RadiationSourceModelType getRadiationSourceModelType() const
    {
        return radiationSourceModelType_;
    }

private:
    RadiationSourceModelType radiationSourceModelType_;
};

class IsotropicPointRadiationSourceModelSettings : public RadiationSourceModelSettings
{
public:
    explicit IsotropicPointRadiationSourceModelSettings(
            const std::shared_ptr<LuminosityModelSettings> luminosityModelSettings) :
            RadiationSourceModelSettings(isotropic_point_source),
            luminosityModelSettings_(luminosityModelSettings) {}

    const std::shared_ptr<LuminosityModelSettings> getLuminosityModelSettings() const
    {
        return luminosityModelSettings_;
    }

private:
    std::shared_ptr<LuminosityModelSettings> luminosityModelSettings_;
};

std::shared_ptr<electromagnetism::LuminosityModel> createLuminosityModel(
        std::shared_ptr< LuminosityModelSettings > modelSettings,
        const std::string& body);

std::shared_ptr<electromagnetism::RadiationSourceModel> createRadiationSourceModel(
        std::shared_ptr< RadiationSourceModelSettings > modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies);


} // tudat
} // simulation_setup

#endif //TUDAT_CREATERADIATIONSOURCEMODEL_H
