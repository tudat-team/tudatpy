#ifndef TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H
#define TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H

#include <memory>

#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"

namespace tudat
{
namespace simulation_setup
{

enum RadiationPressureTargetModelType
{
    cannonball_target,
    paneled_target
};

class RadiationPressureTargetModelSettings
{
public:
    explicit RadiationPressureTargetModelSettings(
            RadiationPressureTargetModelType radiationPressureTargetModelType) :
            radiationPressureTargetModelType_(radiationPressureTargetModelType) {}

    virtual ~RadiationPressureTargetModelSettings() = default;

    RadiationPressureTargetModelType getRadiationPressureTargetModelType() const
    {
        return radiationPressureTargetModelType_;
    }

private:
    RadiationPressureTargetModelType radiationPressureTargetModelType_;
};

class CannonballRadiationPressureTargetModelSettings : public RadiationPressureTargetModelSettings
{
public:
    explicit CannonballRadiationPressureTargetModelSettings(
            double area, double coefficient) :
            RadiationPressureTargetModelSettings(cannonball_target),
            area_(area),
            coefficient_(coefficient) {}

    double getArea() const
    {
        return area_;
    }

    double getCoefficient() const
    {
        return coefficient_;
    }

private:
    double area_;

    double coefficient_;
};

std::shared_ptr<electromagnetism::RadiationPressureTargetModel> createRadiationPressureTargetModel(
        std::shared_ptr< RadiationPressureTargetModelSettings > modelSettings,
        const std::string& body);

} // tudat
} // simulation_setup

#endif //TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H
