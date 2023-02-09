/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/electromagnetism/reflectionLaw.h"

#include <Eigen/Core>

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/mathematicalConstants.h"


namespace tudat
{
namespace electromagnetism
{

double SpecularDiffuseMixReflectionLaw::evaluateReflectedFraction(const Eigen::Vector3d& surfaceNormal,
                                                                  const Eigen::Vector3d& incomingDirection,
                                                                  const Eigen::Vector3d& observerDirection) const
{
    // Check if any reflected radiation would reach observer
    const double cosBetweenNormalAndIncoming = surfaceNormal.dot(-incomingDirection);
    const double cosBetweenNormalAndObserver = surfaceNormal.dot(observerDirection);
    if (cosBetweenNormalAndIncoming <= 0 || cosBetweenNormalAndObserver <= 0)
    {
        // Radiation is incident on backside, or observer is on backside
        return 0;
    }

    // Wetterer (2014) Eq. 4
    const auto diffuseReflectance = diffuseReflectivity_ / mathematical_constants::PI;

    // TODO-DOMINIK is this correct? Should follow from Montenbruck (2014)
    double instantaneousReradiationReflectance = 0;
    if (withInstantaneousReradiation_)
    {
        instantaneousReradiationReflectance = absorptivity_ / mathematical_constants::PI;
    }

    double specularReflectance = 0;
    if (specularReflectivity_ > 0) {
        const auto mirrorOfIncomingDirection = computeMirrorlikeReflection(incomingDirection, surfaceNormal);
        if (observerDirection.isApprox(mirrorOfIncomingDirection))
        {
            // Observer only receives specular reflection if it is in mirrored path of incident radiation
            // Wetterer (2014) Eq. 4
            specularReflectance = specularReflectivity_ / cosBetweenNormalAndIncoming;
        }
    }
    
    return diffuseReflectance + specularReflectance + instantaneousReradiationReflectance;
}

Eigen::Vector3d SpecularDiffuseMixReflectionLaw::evaluateReactionVector(const Eigen::Vector3d& surfaceNormal,
                                                                        const Eigen::Vector3d& incomingDirection) const
{
    const double cosBetweenNormalAndIncoming = surfaceNormal.dot(-incomingDirection);
    if (cosBetweenNormalAndIncoming <= 0)
    {
        // Radiation is incident on backside of surface
        return Eigen::Vector3d::Zero();
    }

    // Montenbruck (2014) Eq. 5
    const Eigen::Vector3d reactionFromIncidence = (absorptivity_ + diffuseReflectivity_) * incomingDirection;
    const Eigen::Vector3d reactionFromReflection =
            -(2. / 3 * diffuseReflectivity_ + 2 * specularReflectivity_ * cosBetweenNormalAndIncoming) * surfaceNormal;

    Eigen::Vector3d reactionFromInstantaneousReradiation;
    if (withInstantaneousReradiation_)
    {
        // Montenbruck (2014) Eq. 6
        // Instantaneous Lambertian reradiation behaves like diffuse Lambertian reflection
        reactionFromInstantaneousReradiation = -(2. / 3 * absorptivity_) * surfaceNormal;
    }
    else
    {
        reactionFromInstantaneousReradiation = Eigen::Vector3d::Zero();
    }

    return reactionFromIncidence + reactionFromReflection + reactionFromInstantaneousReradiation;
}

void SpecularDiffuseMixReflectionLaw::validateCoefficients() const
{
    auto sumOfCoeffs = absorptivity_ + specularReflectivity_ + diffuseReflectivity_;
    if (std::fabs(1 - sumOfCoeffs) >= 20 * std::numeric_limits<double>::epsilon())
    {
        std::cerr << "Warning, coefficients of specular-diffuse-mix reflection law, " <<
                "should sum to 1" << std::endl;
    }
}

Eigen::Vector3d computeMirrorlikeReflection(
        const Eigen::Vector3d& vectorToMirror,
        const Eigen::Vector3d& surfaceNormal)
{
    const double vectorDotNormal = vectorToMirror.dot(surfaceNormal);
    if (vectorDotNormal >= 0)
    {
        // Vector is incident on backside of surface
        return Eigen::Vector3d::Zero();
    }
    else
    {
        return vectorToMirror - 2 * vectorDotNormal * surfaceNormal;
    }
}

} // tudat
} // electromagnetism
