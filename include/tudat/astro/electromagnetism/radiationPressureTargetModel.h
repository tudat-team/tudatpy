/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Oliver Montenbruck, et al. Satellite Orbits. Springer Berlin Heidelberg, 2000.
 */

#ifndef TUDAT_RADIATIONPRESSURETARGETMODEL_H
#define TUDAT_RADIATIONPRESSURETARGETMODEL_H

#include <map>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/electromagnetism/reflectionLaw.h"


namespace tudat
{
namespace electromagnetism
{

/*!
 * Class modeling a target that is accelerated by radiation pressure.
 *
 * Parameters to all functions are in local target-fixed coordinates. Proper translation and rotation is ensured by the
 * RadiationPressureAcceleration class.
 */
class RadiationPressureTargetModel
{
public:
    explicit RadiationPressureTargetModel(
        const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies = {}) :
        sourceToTargetOccultingBodies_(sourceToTargetOccultingBodies) {}

    virtual ~RadiationPressureTargetModel() = default;

    void updateMembers(double currentTime);

    /*!
     * Calculate radiation pressure force from incident radiation using geometrical/optical target properties.
     *
     * @param sourceIrradiance Incident irradiance magnitude [W/m²]
     * @param sourceToTargetDirection Direction of incoming radiation in local (i.e. target-fixed) coordinates
     * @return Radiation pressure force vector in local (i.e. target-fixed) coordinates [N]
     */
    virtual Eigen::Vector3d evaluateRadiationPressureForce(
            double sourceIrradiance, Eigen::Vector3d sourceToTargetDirection) const = 0;
    
    std::map<std::string, std::vector<std::string>> getSourceToTargetOccultingBodies() const
    {
        return sourceToTargetOccultingBodies_;
    }

protected:
    virtual void updateMembers_(const double currentTime) {};

    double currentTime_{TUDAT_NAN};
    // Only needed to transfer occultation settings from body setup to acceleration setup
    std::map<std::string, std::vector<std::string>> sourceToTargetOccultingBodies_;
};

/*!
 * Class modeling a target as sphere ("cannonball"). The sphere has an isotropic radiation pressure coefficient and
 * the same cross-sectional area from any direction.
 */
class CannonballRadiationPressureTargetModel : public RadiationPressureTargetModel
{
public:
    /*!
     * Constructor.
     *
     * @param area Cross-sectional area (i.e. projected area of the sphere)
     * @param coefficient Radiation pressure coefficient (between 1 [pure absorption] and 2 [pure specular reflection])
     * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
     *      to occult sources as seen from this target
     */
    CannonballRadiationPressureTargetModel(
        double area, double coefficient,
        const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies = {}) :
        RadiationPressureTargetModel(sourceToTargetOccultingBodies),
        area_(area), coefficient_(coefficient) {}

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

/*!
 * Class modeling a target as collection of panels, e.g., representing the box body and solar panels.
 */
class PaneledRadiationPressureTargetModel : public RadiationPressureTargetModel
{
public:
    class Panel;

    /*!
     * Constructor.
     *
     * @param panels Panels comprising this paneled target
     * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
     *      to occult sources as seen from this target
     */
    explicit PaneledRadiationPressureTargetModel(
            const std::vector<Panel>& panels,
            const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies = {}) :
            RadiationPressureTargetModel(sourceToTargetOccultingBodies), panels_(panels) {}

    /*!
     * Constructor.
     *
     * @param panels Panels comprising this paneled target
     * @param sourceToTargetOccultingBodies Map (source name -> list of occulting body names) of bodies
     *      to occult sources as seen from this target
     */
    PaneledRadiationPressureTargetModel(
            std::initializer_list<Panel> panels,
            const std::map<std::string, std::vector<std::string>>& sourceToTargetOccultingBodies = {}) :
            RadiationPressureTargetModel(sourceToTargetOccultingBodies), panels_(panels) {}

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

/*!
 * Class modeling a panel in a paneled radiation pressure target. A panel is defined by its area, orientation and
 * reflection law.
 */
class PaneledRadiationPressureTargetModel::Panel
{
    friend class PaneledRadiationPressureTargetModel;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /*!
     * Constructor.
     *
     * @param area Area of the panel [m²]
     * @param surfaceNormalFunction Function returning panel's normal vector in local (i.e. target-fixed) coordinates
     * @param reflectionLaw Reflection law governing how force magnitude and direction due to radiation
     */
    explicit Panel(double area,
                   const std::function<Eigen::Vector3d()>& surfaceNormalFunction,
                   const std::shared_ptr<ReflectionLaw>& reflectionLaw) :
            surfaceNormalFunction_(surfaceNormalFunction),
            reflectionLaw_(reflectionLaw),
            area_(area) {}

    /*!
     * Constructor.
     *
     * @param area Area of the panel [m²]
     * @param surfaceNormal Constant normal vector of the panel in local (i.e. target-fixed) coordinates
     * @param reflectionLaw Reflection law governing how force magnitude and direction due to radiation
     */
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
