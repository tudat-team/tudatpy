/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"


namespace tudat
{
namespace electromagnetism
{

void RadiationPressureTargetModel::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        updateMembers_(currentTime);
    }
}

Eigen::Vector3d CannonballRadiationPressureTargetModel::evaluateRadiationPressureForce(double sourceIrradiance,
                                                                                       Eigen::Vector3d sourceToTargetDirection) const
{
    // From Montenbruck (2000), Sec. 3.4
    const auto radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
    const auto forceMagnitude = coefficient_ * area_ * radiationPressure;
    Eigen::Vector3d force = forceMagnitude * sourceToTargetDirection;
    return force;
}

Eigen::Vector3d PaneledRadiationPressureTargetModel::evaluateRadiationPressureForce(
        double sourceIrradiance,
        Eigen::Vector3d sourceToTargetDirection) const
{
    Eigen::Vector3d force = Eigen::Vector3d::Zero();
    for (auto& panel : panels_)
    {
        const auto surfaceNormal = panel.getSurfaceNormal();
        const double cosBetweenNormalAndIncoming = (-sourceToTargetDirection).dot(surfaceNormal);
        if (cosBetweenNormalAndIncoming >= 0)
        {
            const double effectiveArea = panel.getArea() * cosBetweenNormalAndIncoming;
            const auto radiationPressure = sourceIrradiance / physical_constants::SPEED_OF_LIGHT;
            const auto reactionVector =
                    panel.getReflectionLaw()->evaluateReactionVector(surfaceNormal, sourceToTargetDirection);
            force += radiationPressure * effectiveArea * reactionVector;
        }
    }
    return force;
}

void PaneledRadiationPressureTargetModel::updateMembers_(double currentTime)
{
    for (auto& panel : panels_)
    {
        panel.updateMembers();
    }
}

void PaneledRadiationPressureTargetModel::Panel::updateMembers()
{
    // Evaluate only once per timestep since surface normal function could be expensive to evaluate
    surfaceNormal_ = surfaceNormalFunction_();
}
} // tudat
} // electromagnetism
