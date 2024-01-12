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
 *      Knocke, Philip et al. "Earth radiation pressure effects on satellites."
 *          Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988.
 */

#ifndef TUDAT_SURFACEPROPERTYDISTRIBUTION_H
#define TUDAT_SURFACEPROPERTYDISTRIBUTION_H

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/sphericalHarmonics.h"


namespace tudat
{
namespace electromagnetism
{

using mathematical_constants::PI;

/*!
 * Class modeling the distribution of a property on the surface of a sphere such as albedo or emissivity. The
 * distribution can vary with latitude, longitude and time.
 */
class SurfacePropertyDistribution
{
public:
    explicit SurfacePropertyDistribution() = default;

    virtual ~SurfacePropertyDistribution() = default;

    /*!
     * Update class members.
     *
     * @param currentTime Current simulation time
     */
    void updateMembers(double currentTime);

    virtual double getValue(double latitude, double longitude) = 0;

    virtual bool isTimeInvariant() = 0;

protected:
    virtual void updateMembers_(const double currentTime) {};

    double currentTime_{TUDAT_NAN};
};

/*!
 * Class modeling the distribution of a property on the surface of a sphere, such as albedo or emissivity. The
 * distribution is constant with respect to latitude, longitude and time.
 */
class ConstantSurfacePropertyDistribution : public SurfacePropertyDistribution
{
public:
    explicit ConstantSurfacePropertyDistribution(const double constantValue) :
            constantValue_(constantValue) {}

    double getValue(double latitude, double longitude) override
    {
        return constantValue_;
    }

    bool isTimeInvariant() override
    {
        return true;
    }

private:
    double constantValue_;
};

/*!
 * Class modeling the distribution of a property on the surface of a sphere, such as albedo or emissivity. The
 * distribution is constant with respect to time. Spatial variations are given by a spherical harmonics expansion.
 */
class SphericalHarmonicsSurfacePropertyDistribution : public SurfacePropertyDistribution
{
public:
    /*!
     * Constructor.
     *
     * @param cosineCoefficients Cosine spherical harmonic coefficients (not normalized)
     * @param sineCoefficients Sine spherical harmonic coefficients (not normalized)
     */
    explicit SphericalHarmonicsSurfacePropertyDistribution(
            const Eigen::MatrixXd& cosineCoefficients,
            const Eigen::MatrixXd& sineCoefficients) :
            cosineCoefficients_(cosineCoefficients),
            sineCoefficients_(sineCoefficients),
            maximumDegree_(cosineCoefficients.rows() - 1),
            maximumOrder_( cosineCoefficients.cols() - 1),
            sphericalHarmonicsCache_(maximumDegree_, maximumOrder_, false),
            legendreCache_(sphericalHarmonicsCache_.getLegendreCache())
            {
                if((cosineCoefficients.rows() != sineCoefficients.rows()) ||
                    (cosineCoefficients.cols() != sineCoefficients.cols()))
                {
                    throw std::runtime_error(
                            "Error when creating spherical harmonics surface property distribution; sine and cosine sizes are incompatible" );
                }

                legendreCache_->setComputeFirstDerivatives(false);
            }

    double getValue(double latitude, double longitude) override;

    bool isTimeInvariant() override
    {
        return true;
    }

    const Eigen::MatrixXd& getCosineCoefficients() const
    {
        return cosineCoefficients_;
    }

    const Eigen::MatrixXd& getSineCoefficients() const
    {
        return sineCoefficients_;
    }

    int getMaximumDegree() const
    {
        return maximumDegree_;
    }

    int getMaximumOrder() const
    {
        return maximumOrder_;
    }

private:
    // Cosine spherical harmonic coefficients (not normalized)
    Eigen::MatrixXd cosineCoefficients_;

    // Sine spherical harmonic coefficients (not normalized)
    Eigen::MatrixXd sineCoefficients_;

    const int maximumDegree_;

    const int maximumOrder_;

    basic_mathematics::SphericalHarmonicsCache sphericalHarmonicsCache_;

    std::shared_ptr<basic_mathematics::LegendreCache> legendreCache_;
};

/*!
 * Class modeling the distribution of a property on the surface of a sphere, such as albedo or emissivity. The
 * distribution is periodic. Spatial variations are given by a second-degree zonal spherical harmonics expansion (i.e.
 * the distribution is constant w.r.t. longitude). This model corresponds to the albedo and emissivity for Earth in
 * Knocke (1988).
 *
 * Only the first-degree zonal coefficient varies with time.
 */
class SecondDegreeZonalPeriodicSurfacePropertyDistribution : public SurfacePropertyDistribution
{
public:
    /*!
     * Constructor.
     *
     * @param a0 Zeroeth-degree coefficient
     * @param c0 Constant term of first-degree zonal coefficient
     * @param c1 Cosine coefficient of first-degree zonal coefficient
     * @param c2 Sine coefficient of first-degree zonal coefficient
     * @param a2 Second-degree zonal coefficient
     * @param referenceEpoch Reference epoch for periodicity [seconds]
     * @param period Period of periodicity [days]
     */
    explicit SecondDegreeZonalPeriodicSurfacePropertyDistribution(
            const double a0,
            const double c0,
            const double c1,
            const double c2,
            const double a2,
            const double referenceEpoch,
            const double period) :
            a0(a0),
            c0(c0),
            c1(c1),
            c2(c2),
            a2(a2),
            referenceEpoch(referenceEpoch),
            angularFrequency(2 * PI / period) {}

    double getValue(double latitude, double longitude) override
    {
        return getValue(latitude);
    };

    double getValue(double latitude) const;

    bool isTimeInvariant() override
    {
        return false;
    };

private:
    void updateMembers_(double currentTime) override;

    const double a0;
    double a1{TUDAT_NAN};
    const double c0;
    const double c1;
    const double c2;
    const double a2;
    const double referenceEpoch;
    const double angularFrequency;
};

class CustomSurfacePropertyDistribution : public SurfacePropertyDistribution
{
public:
    CustomSurfacePropertyDistribution(
        const std::function< double( const double, const double, const double ) > customFunction ):
        customFunction_( customFunction ){ }

    double getValue(double latitude, double longitude)
    {
        return customFunction_( latitude, longitude, currentTime_ );
    }

    bool isTimeInvariant( ){ return false; }

private:
    std::function< double( const double, const double, const double ) > customFunction_;
};

} // electromagnetism

} // tudat

#endif //TUDAT_SURFACEPROPERTYDISTRIBUTION_H
