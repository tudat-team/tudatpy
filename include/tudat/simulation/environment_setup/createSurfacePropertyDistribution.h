/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATBUNDLE_CREATESURFACEPROPERTYDISTRIBUTION_H
#define TUDATBUNDLE_CREATESURFACEPROPERTYDISTRIBUTION_H

#include <memory>

#include "tudat/astro/electromagnetism/surfacePropertyDistribution.h"


namespace tudat
{
namespace simulation_setup
{

/*!
 * Types of surface property distributions.
 */
enum class SurfacePropertyDistributionType
{
    constant,
    spherical_harmonics,
    second_degree_zonal_periodic,
    custom_surface_distribution
};

/*!
 * Settings for distribution of a property on the surface of a sphere such as albedo or emissivity.
 *
 * @see SurfacePropertyDistribution
 */
class SurfacePropertyDistributionSettings
{
public:
    explicit SurfacePropertyDistributionSettings(
            const SurfacePropertyDistributionType surfacePropertyDistributionType) :
            surfacePropertyDistributionType_(surfacePropertyDistributionType) {}

    virtual ~SurfacePropertyDistributionSettings() = default;

    SurfacePropertyDistributionType getSurfacePropertyDistributionType() const
    {
        return surfacePropertyDistributionType_;
    }

private:
    SurfacePropertyDistributionType surfacePropertyDistributionType_;
};

/*!
 * Settings for a constant surface property distribution.
 *
 * @see ConstantSurfacePropertyDistribution
 */
class ConstantSurfacePropertyDistributionSettings : public SurfacePropertyDistributionSettings
{
public:
    /*!
     * Constructor.
     *
     * @param constantValue Constant value
     */
    explicit ConstantSurfacePropertyDistributionSettings(
            const double constantValue) :
            SurfacePropertyDistributionSettings(SurfacePropertyDistributionType::constant),
            constantValue_(constantValue) {}

    double getConstantValue() const
    {
        return constantValue_;
    }

private:
    double constantValue_;
};

enum class SphericalHarmonicsSurfacePropertyDistributionModel
{
    custom,
    albedo_dlam1 /**< DLAM-1 lunar albedo model: Floberghagen, R. et al. "Lunar Albedo Force Modeling and its Effect on Low Lunar Orbit and Gravity Field Determination". ASR 23. 4(1999): 733-738. */
};

/*!
 * Settings for a surface property distribution described by a spherical harmonics expansion. The reference frame of the
 * body and spherical harmonics must be identical.
 *
 * @see SphericalHarmonicsSurfacePropertyDistribution
 */
class SphericalHarmonicsSurfacePropertyDistributionSettings : public SurfacePropertyDistributionSettings
{
public:
    /*!
     * Constructor with custom model.
     *
     * @param cosineCoefficients Cosine spherical harmonic coefficients (not normalized)
     * @param sineCoefficients Sine spherical harmonic coefficients (not normalized)
     */
    explicit SphericalHarmonicsSurfacePropertyDistributionSettings(
            const Eigen::MatrixXd& cosineCoefficients,
            const Eigen::MatrixXd& sineCoefficients) :
            SurfacePropertyDistributionSettings(SurfacePropertyDistributionType::spherical_harmonics),
            model_(SphericalHarmonicsSurfacePropertyDistributionModel::custom),
            cosineCoefficients_(cosineCoefficients),
            sineCoefficients_(sineCoefficients) {}

    /*!
    * Constructor with model included in Tudat.
    *
    * @param model Spherical harmonics model to be used
    */
    explicit SphericalHarmonicsSurfacePropertyDistributionSettings(
            SphericalHarmonicsSurfacePropertyDistributionModel model);

    const Eigen::MatrixXd& getCosineCoefficients() const
    {
        return cosineCoefficients_;
    }

    const Eigen::MatrixXd& getSineCoefficients() const
    {
        return sineCoefficients_;
    }

private:
    SphericalHarmonicsSurfacePropertyDistributionModel model_;

    // Cosine spherical harmonic coefficients (not normalized)
    Eigen::MatrixXd cosineCoefficients_;

    // Sine spherical harmonic coefficients (not normalized)
    Eigen::MatrixXd sineCoefficients_;
};

enum class KnockeTypeSurfacePropertyDistributionModel
{
    custom,
    albedo_knocke, /**< Knocke Earth albedo model: Knocke, Philip et al. "Earth radiation pressure effects on satellites." Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988. */
    emissivity_knocke /**< Knocke Earth emissivity model: Knocke, Philip et al. "Earth radiation pressure effects on satellites." Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988. */
};

/*!
 * Settings for a surface property distribution described by a second-degree zonal periodic spherical harmonics expansion. The
 * reference frame of the body and spherical harmonics must be identical.
 *
 * @see SphericalHarmonicsSurfacePropertyDistribution
 */
class SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings : public SurfacePropertyDistributionSettings
{
public:
    /*!
     * Constructor with custom model.
     *
     * @param a0 Zeroeth-degree coefficient
     * @param c0 Constant term of first-degree zonal coefficient
     * @param c1 Cosine coefficient of first-degree zonal coefficient
     * @param c2 Sine coefficient of first-degree zonal coefficient
     * @param a2 Second-degree zonal coefficient
     * @param referenceEpoch Reference epoch for periodicity [seconds]
     * @param period Period of periodicity [days]
     */
    explicit SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings(
            const double a0,
            const double c0,
            const double c1,
            const double c2,
            const double a2,
            const double referenceEpoch,
            const double period) :
            SurfacePropertyDistributionSettings(SurfacePropertyDistributionType::second_degree_zonal_periodic),
            model_(KnockeTypeSurfacePropertyDistributionModel::custom),
            a0(a0),
            c0(c0),
            c1(c1),
            c2(c2),
            a2(a2),
            referenceEpoch(referenceEpoch),
            period(period) {}

    /*!
    * Constructor with model included in Tudat.
    *
    * @param model Model to be used
    */
    explicit SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings(
            KnockeTypeSurfacePropertyDistributionModel model);

    double getA0() const
    {
        return a0;
    }

    double getC0() const
    {
        return c0;
    }

    double getC1() const
    {
        return c1;
    }

    double getC2() const
    {
        return c2;
    }

    double getA2() const
    {
        return a2;
    }

    double getReferenceEpoch() const
    {
        return referenceEpoch;
    }

    double getPeriod() const
    {
        return period;
    }

private:
    KnockeTypeSurfacePropertyDistributionModel model_;

    double a0;
    double c0;
    double c1;
    double c2;
    double a2;
    double referenceEpoch;
    double period;
};


class CustomSurfacePropertyDistributionSettings : public SurfacePropertyDistributionSettings
{
public:
    explicit CustomSurfacePropertyDistributionSettings(
        const std::function< double( const double, const double, const double ) > customFunction ) :
        SurfacePropertyDistributionSettings(SurfacePropertyDistributionType::custom_surface_distribution),
        customFunction_( customFunction ){ }

    std::function< double( const double, const double, const double ) > getCustomFunction( )
    {
        return customFunction_;
    }

private:
    std::function< double( const double, const double, const double ) > customFunction_;

};
/*!
 * Create settings for constant surface property distribution.
 *
 * @param constantValue Constant value
 * @return Shared pointer to settings for a constant surface property distribution.
 */
inline std::shared_ptr<SurfacePropertyDistributionSettings>
        constantSurfacePropertyDistributionSettings(double constantValue)
{
    return std::make_shared< ConstantSurfacePropertyDistributionSettings >(constantValue);
}

/*!
 * Create settings for spherical harmonics surface property distribution from coefficients.
 *
 * @param cosineCoefficients Cosine spherical harmonic coefficients (not normalized)
 * @param sineCoefficients Sine spherical harmonic coefficients (not normalized)
 * @return Shared pointer to settings for a spherical harmonics surface property distribution.
 */
inline std::shared_ptr<SurfacePropertyDistributionSettings>
        sphericalHarmonicsSurfacePropertyDistributionSettings(
                const Eigen::MatrixXd& cosineCoefficients,
                const Eigen::MatrixXd& sineCoefficients)
{
    return std::make_shared< SphericalHarmonicsSurfacePropertyDistributionSettings >(
            cosineCoefficients, sineCoefficients);
}

/*!
 * Create settings for spherical harmonics surface property distribution from model included in Tudat.
 *
 * @param model Spherical harmonics model to be used
 * @return Shared pointer to settings for a spherical harmonics surface property distribution.
 */
inline std::shared_ptr<SurfacePropertyDistributionSettings>
        sphericalHarmonicsSurfacePropertyDistributionSettings(
                SphericalHarmonicsSurfacePropertyDistributionModel model)
{
    return std::make_shared< SphericalHarmonicsSurfacePropertyDistributionSettings >(model);
}

/*!
 * Create settings for second-degree zonal surface property distribution from model included in Tudat.
 *
 * @param model Model to be used
 * @return Shared pointer to settings for a second-degree zonal surface property distribution.
 */
inline std::shared_ptr<SurfacePropertyDistributionSettings>
        secondDegreeZonalPeriodicSurfacePropertyDistributionSettings(
                KnockeTypeSurfacePropertyDistributionModel model)
{
    return std::make_shared< SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings >(model);
}

inline std::shared_ptr<SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings>
manualSecondDegreeZonalPeriodicSurfacePropertyDistributionSettings(
    const double a0,
    const double c0,
    const double c1,
    const double c2,
    const double a2,
    const double referenceEpoch,
    const double period )
{
    return std::make_shared< SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings >(
        a0, c0, c1, c2, a2, referenceEpoch, period
        );
}

inline std::shared_ptr<SurfacePropertyDistributionSettings>
customSurfacePropertyDistributionSettings(
    const std::function< double( const double, const double, const double ) > customFunction )
{
    return std::make_shared< CustomSurfacePropertyDistributionSettings >(customFunction);
}

/*!
 * Create surface property distribution from its settings.
 *
 * @param distributionSettings Settings of the surface property distribution
 * @param body Body to which the surface property distribution belongs
 * @return Shared pointer to surface property distribution
 */
std::shared_ptr<electromagnetism::SurfacePropertyDistribution> createSurfacePropertyDistribution(
        const std::shared_ptr<SurfacePropertyDistributionSettings>& distributionSettings,
        const std::string& body);

} // tudat
} // electromagnetism

#endif //TUDATBUNDLE_CREATESURFACEPROPERTYDISTRIBUTION_H
