#include "tudat/astro/electromagnetism/reflectionLaw.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace electromagnetism
{

double SpecularDiffuseMixReflectionLaw::evaluateReflectedFraction(Eigen::Vector3d surfaceNormal,
                                                                  Eigen::Vector3d incomingDirection,
                                                                  Eigen::Vector3d observerDirection) const
{
    // Check if any reflected radiation would reach observer
    const auto cosBetweenNormalAndIncoming =
            linear_algebra::computeCosineOfAngleBetweenVectors(surfaceNormal, -incomingDirection);
    const auto cosBetweenNormalAndObserver =
            linear_algebra::computeCosineOfAngleBetweenVectors(surfaceNormal, observerDirection);
    if (cosBetweenNormalAndIncoming <= 0 || cosBetweenNormalAndObserver <= 0)
    {
        // Radiation is incident on backside, or observer is on backside
        return 0;
    }

    // Wetterer (2014) Eq. 4
    const auto diffuseReflectance = diffuseReflectivity / mathematical_constants::PI;

    const auto mirrorOfIncomingDirection = computeMirrorlikeReflection(incomingDirection, surfaceNormal);
    if (observerDirection.isApprox(mirrorOfIncomingDirection))
    {
        // Observer only receives specular reflection if it is in mirrored path of incident radiation
        const auto specularReflectance = specularReflectivity_ / cosBetweenNormalAndIncoming;
        return diffuseReflectance + specularReflectance;
    }
    else
    {
        return diffuseReflectance;
    }
}

Eigen::Vector3d SpecularDiffuseMixReflectionLaw::evaluateReactionVector(Eigen::Vector3d surfaceNormal,
                                                                          Eigen::Vector3d incomingDirection) const
{
    const auto cosBetweenNormalAndIncoming =
            linear_algebra::computeCosineOfAngleBetweenVectors(surfaceNormal, -incomingDirection);
    if (cosBetweenNormalAndIncoming <= 0)
    {
        // Radiation is incident on backside of surface
        return Eigen::Vector3d::Zero();
    }

    // Montenbruck (2014) Eq. 5
    const auto reactionFromIncidence = (absorptivity_ + diffuseReflectivity) * incomingDirection;
    const auto reactionFromReflection =
            -(2. / 3 * diffuseReflectivity + 2 * specularReflectivity_ * cosBetweenNormalAndIncoming) * surfaceNormal;
    Eigen::Vector3d reactionFromInstantaneousReradiation;
    if (withInstantaneousLambertianReradiation_)
    {
        // Montenbruck (2014) Eq. 6
        // Instantaneous Lambertian reradiation behaves similarly to diffuse Lambertian reflection
        reactionFromInstantaneousReradiation = -(2. / 3 * absorptivity_) * surfaceNormal;
    }
    else
    {
        reactionFromInstantaneousReradiation = Eigen::Vector3d::Zero();
    }

    return reactionFromIncidence + reactionFromReflection + reactionFromInstantaneousReradiation;
}

Eigen::Vector3d computeMirrorlikeReflection(
        const Eigen::Vector3d& vectorToMirror,
        const Eigen::Vector3d& surfaceNormal)
{
    const auto vectorDotNormal = vectorToMirror.dot(surfaceNormal);
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
