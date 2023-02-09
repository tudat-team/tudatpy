/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/electromagnetism/surfacePropertyDistribution.h"

#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/astro/basic_astro/physicalConstants.h"


namespace tudat
{
namespace electromagnetism
{

void SurfacePropertyDistribution::updateMembers(const double currentTime)
{
    if(currentTime_ != currentTime)
    {
        currentTime_ = currentTime;
        updateMembers_(currentTime);
    }
}

double SphericalHarmonicsSurfacePropertyDistribution::getValue(
        double latitude,
        double longitude)
{
    auto sinOfLat = std::sin(latitude);
    sphericalHarmonicsCache_.update(TUDAT_NAN, sinOfLat, longitude, TUDAT_NAN);

    double value = 0;
    double legendrePolynomial;
    for( int degree = 0; degree <= maximumDegree_; degree++ )
    {
        for( int order = 0; order <= degree; order++ )
        {
            legendrePolynomial = basic_mathematics::computeLegendrePolynomialFromCache(
                    degree, order, legendreCache_);
            value += legendrePolynomial * (
                    cosineCoefficients_(degree, order) * sphericalHarmonicsCache_.getCosineOfMultipleLongitude(order) +
                    sineCoefficients_(degree, order) * sphericalHarmonicsCache_.getSineOfMultipleLongitude(order)
            );
        }
    }

    return value;
}

double SecondDegreeZonalPeriodicSurfacePropertyDistribution::getValue(double latitude) const
{
    const auto sinOfLatitude = sin(latitude);

    // Compute first-degree and second-degree zonal Legendre polynomials
    const double P1 = sinOfLatitude;
    const double P2 = (3 * sinOfLatitude * sinOfLatitude - 1) / 2;

    const double value = a0 + a1 * P1 + a2 * P2;

    return value;
}

void SecondDegreeZonalPeriodicSurfacePropertyDistribution::updateMembers_(const double currentTime)
{
    const double daysSinceReferenceEpoch = (currentTime - referenceEpoch) / physical_constants::JULIAN_DAY;
    a1 = c0 + c1 * cos(angularFrequency * daysSinceReferenceEpoch) + c2 * sin(angularFrequency * daysSinceReferenceEpoch);
}

} // tudat
} // electromagnetism
