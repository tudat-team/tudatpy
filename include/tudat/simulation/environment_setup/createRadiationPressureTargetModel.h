#ifndef TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H
#define TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H

#include <memory>

#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/simulation/environment_setup/body.h"

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
            double area,
            double coefficient) :
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
            const std::vector<Panel>& panels) :
            RadiationPressureTargetModelSettings(paneled_target),
            panels_(panels) {}

    const std::vector<Panel>& getPanels() const
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
                   const std::function<Eigen::Vector3d()>& surfaceNormalFunction) :
            area_(area),
            specularReflectivity_(specularReflectivity),
            diffuseReflectivity_(diffuseReflectivity),
            surfaceNormalFunction_(surfaceNormalFunction)  {}

    explicit Panel(double area,
                   double specularReflectivity,
                   double diffuseReflectivity,
                   const Eigen::Vector3d& surfaceNormal) :
            Panel(area,
                  specularReflectivity,
                  diffuseReflectivity,
                  [=] () { return surfaceNormal; }) {}

    explicit Panel(double area,
                   double specularReflectivity,
                   double diffuseReflectivity,
                   const std::string& bodyToTrack,
                   const bool towardsTrackedBody = true) :
            area_(area),
            specularReflectivity_(specularReflectivity),
            diffuseReflectivity_(diffuseReflectivity),
            bodyToTrack_(bodyToTrack),
            towardsTrackedBody_(towardsTrackedBody) {}

    double getArea() const
    {
        return area_;
    }

    const std::function<Eigen::Vector3d()>& getSurfaceNormalFunction() const
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

    const std::string &getBodyToTrack() const
    {
        return bodyToTrack_;
    }

    bool isTowardsTrackedBody() const
    {
        return towardsTrackedBody_;
    }

private:
    double area_;
    double specularReflectivity_;
    double diffuseReflectivity_;
    std::function<Eigen::Vector3d()> surfaceNormalFunction_;
    std::string bodyToTrack_;
    bool towardsTrackedBody_{true};
};

inline std::shared_ptr<CannonballRadiationPressureTargetModelSettings>
        cannonballRadiationPressureTargetModelSettings(double area, double coefficient)
{
    return std::make_shared< CannonballRadiationPressureTargetModelSettings >(area, coefficient);
}

inline std::shared_ptr<PaneledRadiationPressureTargetModelSettings>
        paneledRadiationPressureTargetModelSettings(std::initializer_list<PaneledRadiationPressureTargetModelSettings::Panel> panels)
{
    return std::make_shared< PaneledRadiationPressureTargetModelSettings >(
            std::vector<PaneledRadiationPressureTargetModelSettings::Panel>(panels));
}

std::shared_ptr<electromagnetism::RadiationPressureTargetModel> createRadiationPressureTargetModel(
        const std::shared_ptr< RadiationPressureTargetModelSettings >& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies);

typedef PaneledRadiationPressureTargetModelSettings::Panel TargetPanelSettings;

} // tudat
} // simulation_setup

#endif //TUDAT_CREATERADIATIONPRESSURETARGETMODEL_H
