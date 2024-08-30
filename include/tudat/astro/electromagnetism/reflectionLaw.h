/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      O. Montenbruck, et al. "Enhanced solar radiation pressure modeling for Galileo satellites".
 *          Journal of Geodesy 89. 3(2014): 283–297.
 *      Kristin Vielberg, et al. "Extended forward and inverse modeling of radiation pressure accelerations for
 *          LEO satellites". Journal of Geodesy 94. 4(2020).
 *      Charles J. Wetterer, et al. "Refining Space Object Radiation Pressure Modeling with Bidirectional Reflectance
 *          Distribution Functions". Journal of Guidance, Control, and Dynamics 37. 1(2014): 185–196.
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

/*!
 * Class modeling reflective properties of a surface. The reflection law determines direction and magnitude of reflected
 * radiation (for sources) as well as the reaction force of the surface (for targets).
 */
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
     * @return Reflected fraction of radiation per steradian [1/sr]
     */
    virtual double evaluateReflectedFraction(
            const Eigen::Vector3d& surfaceNormal,
            const Eigen::Vector3d& incomingDirection,
            const Eigen::Vector3d& observerDirection) const = 0;

    /*!
     * Evaluate direction of reaction force due to incident and reflected radiation. For pure absorption, the force is
     * along the incoming radiation. For pure specular reflection, the force is normal to the surface. If the incoming
     * vector is incident on the backside of the surface, the zero vector is returned.
     *
     * @param surfaceNormal Surface normal unit vector
     * @param incomingDirection Incoming radiation unit vector from source to target
     * @return Reaction vector
     */
    virtual Eigen::Vector3d evaluateReactionVector(
            const Eigen::Vector3d& surfaceNormal,
            const Eigen::Vector3d& incomingDirection) const = 0;

    virtual Eigen::Vector3d evaluateReactionVectorPartialWrtSpecularReflectivity(
            const Eigen::Vector3d& surfaceNormal,
            const Eigen::Vector3d& incomingDirection) const = 0;

    virtual Eigen::Vector3d evaluateReactionVectorPartialWrtDiffuseReflectivity(
            const Eigen::Vector3d& surfaceNormal,
            const Eigen::Vector3d& incomingDirection) const = 0;

    virtual Eigen::Matrix3d evaluateReactionVectorDerivativeWrtTargetPosition(
        const Eigen::Vector3d& surfaceNormal,
        const Eigen::Vector3d& incomingDirection,
        const double cosineOfAngleBetweenVectors,
        const Eigen::Vector3d& currentReactionVector,
        const Eigen::Matrix3d& sourceUnitVectorPartial,
        const Eigen::Matrix< double, 1, 3 >& cosineAnglePartial ) = 0;
};

/*!
 * Reflection law for a mix of purely specular and purely diffuse Lambertian reflection. All radiation that is not
 * reflected is absorbed. Absorbed radiation can (optionally) be instantaneously diffusely reradiated, e.g., as
 * simplified model of a satellite's own thermal radiation (Vielberg, 2020). Instantaneous reradiation only makes
 * sense for targets, for sources use the dedicated thermal panel radiosity models.
 * 
 * Note that, to receive any specular reflection in evaluateReflectedFraction(), the observer has to be exactly
 * in the mirror-like reflection of the incoming ray, which is unlikely.
 * 
 * Model is described by Wetterer (2014) Eq. 3.
 */
class SpecularDiffuseMixReflectionLaw : public ReflectionLaw
{
public:
    /*!
     * Constructor.
     *
     * @param absorptivity Absorptivity (between 0 and 1)
     * @param specularReflectivity Specular reflectivity (between 0 and 1)
     * @param diffuseReflectivity Diffuse reflectivity (between 0 and 1)
     * @param withInstantaneousReradiation Whether to instantaneously reradiate absorbed radiation
     */
    explicit SpecularDiffuseMixReflectionLaw(
            double absorptivity,
            double specularReflectivity,
            double diffuseReflectivity,
            bool withInstantaneousReradiation = false) :
            absorptivity_(absorptivity),
            specularReflectivity_(specularReflectivity),
            diffuseReflectivity_(diffuseReflectivity),
            withInstantaneousReradiation_(withInstantaneousReradiation)
    {
        validateCoefficients();
    }

    double evaluateReflectedFraction(
            const Eigen::Vector3d& surfaceNormal,
            const Eigen::Vector3d& incomingDirection,
            const Eigen::Vector3d& observerDirection) const override;

    Eigen::Vector3d evaluateReactionVector(
        const Eigen::Vector3d& surfaceNormal,
        const Eigen::Vector3d& incomingDirection) const override;

    Eigen::Vector3d evaluateReactionVectorPartialWrtSpecularReflectivity(
        const Eigen::Vector3d& surfaceNormal,
        const Eigen::Vector3d& incomingDirection) const override;

    Eigen::Vector3d evaluateReactionVectorPartialWrtDiffuseReflectivity(
        const Eigen::Vector3d& surfaceNormal,
        const Eigen::Vector3d& incomingDirection) const override;

    Eigen::Matrix3d evaluateReactionVectorDerivativeWrtTargetPosition(
        const Eigen::Vector3d& surfaceNormal,
        const Eigen::Vector3d& incomingDirection,
        const double cosineOfAngleBetweenVectors,
        const Eigen::Vector3d& currentReactionVector,
        const Eigen::Matrix3d& sourceUnitVectorPartial,
        const Eigen::Matrix< double, 1, 3 >& cosineAnglePartial ) override;

    double getAbsorptivity() const
    {
        return absorptivity_;
    }

    void setAbsorptivity( const double absorptivity )
    {
        absorptivity_ = absorptivity;
    }

    double getSpecularReflectivity() const
    {
        return specularReflectivity_;
    }

    void setSpecularReflectivity( const double specularReflectivity)
    {
        specularReflectivity_ = specularReflectivity;
    }

    double getDiffuseReflectivity() const
    {
        return diffuseReflectivity_;
    }

    void setDiffuseReflectivity( const double diffuseReflectivity)
    {
        diffuseReflectivity_ = diffuseReflectivity;
    }

    bool isWithInstantaneousReradiation() const
    {
        return withInstantaneousReradiation_;
    }

protected:
    void validateCoefficients() const;

    double absorptivity_;
    double specularReflectivity_;
    double diffuseReflectivity_;
    bool withInstantaneousReradiation_;
};

/*!
 * A reflection law for purely diffuse Lambertian reflection.
 */
class LambertianReflectionLaw : public SpecularDiffuseMixReflectionLaw
{
public:
    /*!
     * Constructor.

     * @param diffuseReflectivity Diffuse reflectivity, e.g., albedo (between 0 and 1)
     */
    explicit LambertianReflectionLaw(
            double diffuseReflectivity) :
        SpecularDiffuseMixReflectionLaw(
            1-diffuseReflectivity,
            0,
            diffuseReflectivity,
            false) {}

    void setDiffuseReflectivity(double diffuseReflectivity)
    {
        absorptivity_ = 1 - diffuseReflectivity;
        diffuseReflectivity_ = diffuseReflectivity;
    }
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


inline std::shared_ptr<SpecularDiffuseMixReflectionLaw> reflectionLawFromAbsorptivityAndSpecularReflectivity(
        double absorptivity,
        double specularReflectivity,
        bool withInstantaneousReradiation = false)
{
    const auto diffuseReflectivity = 1.0 - absorptivity - specularReflectivity;
    return std::make_shared<SpecularDiffuseMixReflectionLaw>(
            absorptivity,
            specularReflectivity,
            diffuseReflectivity,
            withInstantaneousReradiation);
}

inline std::shared_ptr<SpecularDiffuseMixReflectionLaw> reflectionLawFromAbsorptivityAndDiffuseReflectivity(
        double absorptivity,
        double diffuseReflectivity,
        bool withInstantaneousReradiation = false)
{
    const auto specularReflectivity = 1.0 - absorptivity - diffuseReflectivity;
    return std::make_shared<SpecularDiffuseMixReflectionLaw>(
            absorptivity,
            specularReflectivity,
            diffuseReflectivity,
            withInstantaneousReradiation);
}

inline std::shared_ptr<SpecularDiffuseMixReflectionLaw> reflectionLawFromSpecularAndDiffuseReflectivity(
        double specularReflectivity,
        double diffuseReflectivity,
        bool withInstantaneousReradiation = false)
{
    const auto absorptivity = 1.0 - specularReflectivity - diffuseReflectivity;
    return std::make_shared<SpecularDiffuseMixReflectionLaw>(
            absorptivity,
            specularReflectivity,
            diffuseReflectivity,
            withInstantaneousReradiation);
}


} // tudat
} // electromagnetism

#endif //TUDAT_REFLECTIONLAW_H
