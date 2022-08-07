#ifndef TUDAT_RADIATIONPRESSURETARGETMODEL_H
#define TUDAT_RADIATIONPRESSURETARGETMODEL_H

#include <vector>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/electromagnetism/reflectionLaw.h"

namespace tudat
{
namespace electromagnetism
{

// All calculations in body-fixed frame
class RadiationPressureTargetModel
{
public:
    virtual ~RadiationPressureTargetModel() = default;

    virtual Eigen::Vector3d evaluateRadiationPressureForce(
            double sourceIrradiance, Eigen::Vector3d sourceToTargetDirection) const = 0;
};


class CannonballRadiationPressureTargetModel : public RadiationPressureTargetModel
{
public:
    CannonballRadiationPressureTargetModel(double area, double coefficient)
            : area_(area), coefficient_(coefficient) {}

    Eigen::Vector3d evaluateRadiationPressureForce(
            double sourceIrradiance,
            Eigen::Vector3d sourceToTargetDirection) const override;

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


class PaneledRadiationPressureTargetModel : public RadiationPressureTargetModel
{
public:
    class Panel;

    explicit PaneledRadiationPressureTargetModel(
            std::vector<Panel> panels)
            : panels_(panels) {}

    PaneledRadiationPressureTargetModel(
            std::initializer_list<Panel> panels)
            : panels_(panels) {}

    Eigen::Vector3d evaluateRadiationPressureForce(
            double sourceIrradiance,
            Eigen::Vector3d sourceToTargetDirection) const override;

    const std::vector<Panel> getPanels() const
    {
        return panels_;
    }

private:
    std::vector<Panel> panels_;
};

class PaneledRadiationPressureTargetModel::Panel
{
public:
    explicit Panel(double area,
                   const std::function<Eigen::Vector3d()> surfaceNormalFunction,
                   const std::shared_ptr<ReflectionLaw> reflectionLaw) :
            area_(area),
            surfaceNormalFunction_(surfaceNormalFunction),
            reflectionLaw_(reflectionLaw) {}

    explicit Panel(double area,
                   const Eigen::Vector3d surfaceNormal,
                   const std::shared_ptr<ReflectionLaw> reflectionLaw) :
            Panel(area, [=] () { return surfaceNormal; }, reflectionLaw) {}

    double getArea() const
    {
        return area_;
    }

    Eigen::Vector3d getSurfaceNormal() const
    {
        return surfaceNormalFunction_();
    }

    std::shared_ptr<ReflectionLaw> getReflectionLaw() const
    {
        return reflectionLaw_;
    }

private:
    double area_;
    std::function<Eigen::Vector3d()> surfaceNormalFunction_;
    std::shared_ptr<ReflectionLaw> reflectionLaw_;
};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONPRESSURETARGETMODEL_H
