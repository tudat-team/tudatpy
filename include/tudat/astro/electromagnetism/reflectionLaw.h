/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Charles J. Wetterer, et al. "Refining Space Object Radiation Pressure Modeling with Bidirectional Reflectance
 *          Distribution Functions". Journal of Guidance, Control, and Dynamics 37. 1(2014): 185–196.
 *      O. Montenbruck, et al. "Enhanced solar radiation pressure modeling for Galileo satellites".
 *          Journal of Geodesy 89. 3(2014): 283–297.
 */

#ifndef TUDAT_REFLECTIONLAW_H
#define TUDAT_REFLECTIONLAW_H

#include <iostream>
#include <memory>

#include <Eigen/Core>

namespace tudat
{
namespace electromagnetism
{

class ReflectionLaw
{
public:
    ReflectionLaw() = default;

    virtual ~ReflectionLaw() = default;

    /*!
     * Evaluate reflected fraction (per steradian) of incoming light in a certain observer direction. Such a function
     * is called Bidirectional Reflectance Distribution Function (BRDF). For example, this is important for albedo
     * calculations.
     *
     * All vectors must be normalized.
     *
     * @param surfaceNormal Surface normal unit vector
     * @param incomingDirection Incoming radiation unit vector from source to target
     * @param observerDirection Unit vector from surface in direction of observer
     * @return The reflected fraction of radiation per steradian [1/sr]
     */
    virtual double evaluateReflectedFraction(
            const Eigen::Vector3d& surfaceNormal,
            const Eigen::Vector3d& incomingDirection,
            const Eigen::Vector3d& observerDirection) const = 0;

    /*!
     * Evaluate direction of reaction force due to incident and reflected radiation. If the incoming vector is incident
     * on the backside of the surface, the zero vector is returned.
     *
     * @param surfaceNormal Surface normal unit vector
     * @param incomingDirection Incoming radiation unit vector from source to target
     * @return The vector containing direction and magnitude of reaction due to surface properties [-]
     */
    virtual Eigen::Vector3d evaluateReactionVector(
            const Eigen::Vector3d& surfaceNormal,
            const Eigen::Vector3d& incomingDirection) const = 0;
};

/*!
 * A reflection law for a mix of purely specular and purely diffuse Lambertian reflection.
 */
class SpecularDiffuseMixReflectionLaw : public ReflectionLaw
{
public:
    explicit SpecularDiffuseMixReflectionLaw(
            double absorptivity,
            double specularReflectivity,
            double diffuseReflectivity,
            bool withInstantaneousLambertianReradiation = false) :
        absorptivity_(absorptivity),
        specularReflectivity_(specularReflectivity),
        diffuseReflectivity(diffuseReflectivity),
        withInstantaneousLambertianReradiation_(withInstantaneousLambertianReradiation)
    {
        auto sumOfCoeffs = absorptivity + specularReflectivity + diffuseReflectivity;
        if (std::fabs(1 - sumOfCoeffs) >= 20 * std::numeric_limits<double>::epsilon())
        {
            std::cerr << "Warning when creating specular-diffuse-mix reflection law, " <<
                         "coefficients should sum to 1" << std::endl;
        }
    }

    static std::shared_ptr<SpecularDiffuseMixReflectionLaw> fromAbsorptivityAndSpecularReflectivity(
            double absorptivity,
            double specularReflectivity)
    {
        const auto diffuseReflectivity = 1 - absorptivity - specularReflectivity;
        return std::make_shared<SpecularDiffuseMixReflectionLaw>(
                absorptivity,
                specularReflectivity,
                diffuseReflectivity);
    }

    static std::shared_ptr<SpecularDiffuseMixReflectionLaw> fromAbsorptivityAndDiffuseReflectivity(
            double absorptivity,
            double diffuseReflectivity)
    {
        const auto specularReflectivity = 1 - absorptivity - diffuseReflectivity;
        return std::make_shared<SpecularDiffuseMixReflectionLaw>(
                absorptivity,
                specularReflectivity,
                diffuseReflectivity);
    }

    static std::shared_ptr<SpecularDiffuseMixReflectionLaw> fromSpecularAndDiffuseReflectivity(
            double specularReflectivity,
            double diffuseReflectivity)
    {
        const auto absorptivity = 1 - specularReflectivity - diffuseReflectivity;
        return std::make_shared<SpecularDiffuseMixReflectionLaw>(
                absorptivity,
                specularReflectivity,
                diffuseReflectivity);
    }

    double evaluateReflectedFraction(
            const Eigen::Vector3d& surfaceNormal,
            const Eigen::Vector3d& incomingDirection,
            const Eigen::Vector3d& observerDirection) const override;

    Eigen::Vector3d evaluateReactionVector(
            const Eigen::Vector3d& surfaceNormal,
            const Eigen::Vector3d& incomingDirection) const override;

    double getAbsorptivity() const
    {
        return absorptivity_;
    }

    double getSpecularReflectivity() const
    {
        return specularReflectivity_;
    }

    double getDiffuseReflectivity() const
    {
        return diffuseReflectivity;
    }

private:
    double absorptivity_;

    double specularReflectivity_;

    double diffuseReflectivity;

    bool withInstantaneousLambertianReradiation_;
};

/*!
 * A reflection law for purely diffuse Lambertian reflection.
 */
class LambertianReflectionLaw : public SpecularDiffuseMixReflectionLaw
{
public:
    explicit LambertianReflectionLaw(
            double diffuseReflectivity) :
        SpecularDiffuseMixReflectionLaw(1-diffuseReflectivity, 0, diffuseReflectivity) {}
};

/*!
 * Computed the mirrorlike reflection of an incoming vector on a surface. If the incoming vector is incident on the
 * backside of the surface, the zero vector is returned.
 *
 * @param vectorToMirror Incoming vector to mirror, pointing towards the surface
 * @param surfaceNormal Surface normal unit vector
 * @return Mirrorlike reflection of vectorToMirror, pointing away from surface
 */
Eigen::Vector3d computeMirrorlikeReflection(
        const Eigen::Vector3d& vectorToMirror,
        const Eigen::Vector3d& surfaceNormal);

} // tudat
} // electromagnetism

#endif //TUDAT_REFLECTIONLAW_H
