/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

/*!
 * Class modeling the distribution of a property on the surface of a sphere, such as albedo or emissivity. The
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

private:
    double constantValue_;
};

/*!
 * Class modeling the distribution of a property on the surface of a sphere, such as albedo or emissivity. The
 * distribution is constant with respect to time. Spatial variations are given by spherical harmonics.
 */
 // A single distribution per paneled source -> shared between all panels
 // source updateMembers calls distribution update members
 // property is queried every time step -> store in map with lat/lon std::pair as key
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
            maximumOrder_( cosineCoefficients.cols() - 1)
            {
                if((cosineCoefficients.rows() != sineCoefficients.rows()) ||
                    (cosineCoefficients.cols() != sineCoefficients.cols()))
                {
                    throw std::runtime_error( "Error when creating spherical harmonics gravity field; sine and cosine sizes are incompatible" );
                }

                sphericalHarmonicsCache_.resetMaximumDegreeAndOrder(maximumDegree_ + 2, maximumOrder_ + 2);
            }

    double getValue(double latitude, double longitude) override;

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

    basic_mathematics::SphericalHarmonicsCache sphericalHarmonicsCache_{false};
};

} // tudat
} // electromagnetism


#endif //TUDAT_SURFACEPROPERTYDISTRIBUTION_H
