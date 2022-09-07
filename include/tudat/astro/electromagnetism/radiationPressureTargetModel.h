/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

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

    void updateMembers(double currentTime);

    virtual Eigen::Vector3d evaluateRadiationPressureForce(
            double sourceIrradiance, Eigen::Vector3d sourceToTargetDirection) const = 0;

protected:
    virtual void updateMembers_(const double currentTime) {};

    double currentTime_{TUDAT_NAN};
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
            const std::vector<Panel>& panels)
            : panels_(panels) {}

    PaneledRadiationPressureTargetModel(
            std::initializer_list<Panel> panels)
            : panels_(panels) {}

    Eigen::Vector3d evaluateRadiationPressureForce(
            double sourceIrradiance,
            Eigen::Vector3d sourceToTargetDirection) const override;

    const std::vector<Panel>& getPanels() const
    {
        return panels_;
    }

private:
    void updateMembers_(double currentTime) override;

    std::vector<Panel> panels_;
};

class PaneledRadiationPressureTargetModel::Panel
{
    friend class PaneledRadiationPressureTargetModel;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    explicit Panel(double area,
                   const std::function<Eigen::Vector3d()>& surfaceNormalFunction,
                   const std::shared_ptr<ReflectionLaw>& reflectionLaw) :
            surfaceNormalFunction_(surfaceNormalFunction),
            reflectionLaw_(reflectionLaw),
            area_(area) {}

    explicit Panel(double area,
                   const Eigen::Vector3d& surfaceNormal,
                   const std::shared_ptr<ReflectionLaw>& reflectionLaw) :
            Panel(area, [=] () { return surfaceNormal; }, reflectionLaw) {}

    double getArea() const
    {
        return area_;
    }

    Eigen::Vector3d getSurfaceNormal() const
    {
        return surfaceNormal_;
    }

    std::shared_ptr<ReflectionLaw> getReflectionLaw() const
    {
        return reflectionLaw_;
    }

private:
    void updateMembers();

    Eigen::Vector3d surfaceNormal_;
    std::function<Eigen::Vector3d()> surfaceNormalFunction_;
    std::shared_ptr<ReflectionLaw> reflectionLaw_;
    double area_;
};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONPRESSURETARGETMODEL_H
