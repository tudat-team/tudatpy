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

    auto legendreCache = *sphericalHarmonicsCache_.getLegendreCache();

    double property = 0;
    double legendrePolynomial;
    for( int degree = 0; degree <= maximumDegree_; degree++ )
    {
        for( int order = 0; order <= degree; order++ )
        {
            legendrePolynomial = basic_mathematics::computeLegendrePolynomialFromCache(
                    degree, order, legendreCache);
            property += legendrePolynomial * (
                    cosineCoefficients_(degree, order) * sphericalHarmonicsCache_.getCosineOfMultipleLongitude(order) +
                    sineCoefficients_(degree, order) * sphericalHarmonicsCache_.getSineOfMultipleLongitude(order)
            );
        }
    }

    return property;
}

} // tudat
} // electromagnetism
