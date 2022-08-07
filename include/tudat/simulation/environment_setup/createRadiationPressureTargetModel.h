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

class PaneledRadiationPressureTargetModelSettings : public RadiationPressureTargetModelSettings
{
public:
    class Panel;

    explicit PaneledRadiationPressureTargetModelSettings(
            std::vector<Panel> panels) :
            RadiationPressureTargetModelSettings(paneled_target),
            panels_(panels) {}

    const std::vector<Panel> getPanels() const
    {
        return panels_;
    }

private:
    std::vector<Panel> panels_;
};

class PaneledRadiationPressureTargetModelSettings::Panel
{
public:
    explicit Panel(double area,
                   double specularReflectivity,
                   double diffuseReflectivity,
                   const std::function<Eigen::Vector3d()> surfaceNormalFunction) :
            area_(area),
            specularReflectivity_(specularReflectivity),
            diffuseReflectivity_(diffuseReflectivity),
            surfaceNormalFunction_(surfaceNormalFunction)  {}

    explicit Panel(double area,
                   double specularReflectivity,
                   double diffuseReflectivity,
                   const Eigen::Vector3d surfaceNormal) :
            Panel(area,
                  specularReflectivity,
                  diffuseReflectivity,
                  [=] () { return surfaceNormal; }) {}

    double getArea() const
    {
        return area_;
    }

    std::function<Eigen::Vector3d()> getSurfaceNormalFunction() const
    {
        return surfaceNormalFunction_;
    }

    double getSpecularReflectivity() const
    {
        return specularReflectivity_;
    }

    double getDiffuseReflectivity() const
    {
        return diffuseReflectivity_;
    }

private:
    double area_;
    double specularReflectivity_;
    double diffuseReflectivity_;
    std::function<Eigen::Vector3d()> surfaceNormalFunction_;
};

inline std::shared_ptr<PaneledRadiationPressureTargetModelSettings>
        paneledRadiationPressureTargetModelSettings(std::initializer_list<PaneledRadiationPressureTargetModelSettings::Panel> panels)
{
    return std::make_shared< PaneledRadiationPressureTargetModelSettings >(
            std::vector<PaneledRadiationPressureTargetModelSettings::Panel>(panels));
}

std::shared_ptr<electromagnetism::RadiationPressureTargetModel> createRadiationPressureTargetModel(
        std::shared_ptr< RadiationPressureTargetModelSettings > modelSettings,
        const std::string& body);

typedef PaneledRadiationPressureTargetModelSettings::Panel TargetPanelSettings;

} // tudat
} // simulation_setup

#endif //TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H
